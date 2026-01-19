#!/usr/bin/env python3
import os
import shutil
import argparse
import subprocess
import numpy as np

try:
    import tomllib  # Python 3.11+
except ImportError:
    try:
        import toml
    except ImportError:
        print("Error: 'toml' or 'tomllib' (Python 3.11+) is required.")
        print("Please install it using: pip install toml")
        exit(1)

class QESetup:
    def __init__(self, config_file="input.toml", struct_file=None):
        self.work_dir = "scf"
        # 1. 加载配置文件
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"找不到配置文件: {config_file}")
        
        if 'tomllib' in globals():
            with open(config_file, "rb") as f:
                self.config = tomllib.load(f)
        else:
            self.config = toml.load(config_file)
        
        self.struct_file = struct_file
        
        # 2. 默认 QE 参数模板 (Namelists)
        self.qe_params = {
            "CONTROL": {
                "calculation": "'scf'",
                "restart_mode": "'from_scratch'",
                "prefix": "'qe'",
                "outdir": "'./out'",
                "pseudo_dir": "'./pseudo'",
                "tprnfor": ".true.",
                "tstress": ".true.",
            },
            "SYSTEM": {
                "ecutwfc": 50.0,
                "ecutrho": 400.0,
                "occupations": "'smearing'",
                "smearing": "'gaussian'",
                "degauss": 0.01,
            },
            "ELECTRONS": {
                "conv_thr": 1.0e-8,
                "mixing_beta": 0.7,
            }
        }

    def parse_poscar(self, poscar_path):
        """从 POSCAR 中读取结构信息"""
        if not os.path.exists(poscar_path):
            raise FileNotFoundError(f"找不到 {poscar_path} 文件")
        
        with open(poscar_path, "r") as f:
            lines = f.readlines()
            
        # 第2行: 缩放因子
        scale = float(lines[1].strip())
        
        # 第3-5行: 晶格矢量
        lattice = []
        for i in range(2, 5):
            lattice.append([float(x) * scale for x in lines[i].split()])
            
        # 第6行: 元素名称
        elements = lines[5].split()
        if elements[0].isdigit():
            # 简单处理 VASP 4 格式，假设元素名在第1行注释中
            print("警告: POSCAR 可能是 VASP 4 格式，尝试从第1行提取元素。")
            elements = lines[0].split()
            if not elements:
                raise ValueError("无法在 POSCAR 中找到元素名称。")
        
        # 第7行: 原子数量
        counts = [int(x) for x in lines[6].split()]
        
        # 第8行: 坐标类型 (Direct 或 Cartesian)
        coord_line = lines[7].strip()
        if coord_line.lower().startswith('s'): # Selective dynamics
            coord_line = lines[8].strip()
            start_line = 9
        else:
            start_line = 8
            
        # 第9行开始: 坐标
        positions = []
        total_atoms = sum(counts)
        for i in range(total_atoms):
            positions.append(lines[start_line + i].split()[:3])
            
        return {
            "lattice": np.array(lattice),
            "elements": elements,
            "counts": counts,
            "coord_type": coord_line,
            "positions": np.array(positions, dtype=float)
        }

    def get_kpoints(self, lattice, kmesh):
        """计算 K 点网格"""
        if isinstance(kmesh, list) and len(kmesh) == 3:
            return kmesh
        
        # 计算倒格矢长度
        b = 2 * np.pi * np.linalg.inv(lattice).T
        b_lengths = np.linalg.norm(b, axis=1)
        
        # 如果 kmesh 是密度 (spacing)
        # 默认 spacing 0.04 * 2pi / Angstrom
        spacing = kmesh if isinstance(kmesh, float) else 0.04
        kpts = [int(max(1, np.ceil(l / (spacing * 2 * np.pi)))) for l in b_lengths]
        return kpts

    def generate_qe_input(self, struct_info):
        """生成 QE pw.x 输入文件"""
        print(f"正在生成 qe.in...")
        
        # 合并配置文件中的参数
        if "qe_params" in self.config:
            for section, params in self.config["qe_params"].items():
                section_upper = section.upper()
                if section_upper in self.qe_params:
                    self.qe_params[section_upper].update(params)
                else:
                    self.qe_params[section_upper] = params

        qe_input_path = os.path.join(self.work_dir, "qe.in")
        with open(qe_input_path, "w") as f:
            # Namelists
            for section in ["CONTROL", "SYSTEM", "ELECTRONS"]:
                f.write(f"&{section}\n")
                if section == "SYSTEM":
                    f.write(f"  nat  = {struct_info['positions'].shape[0]}\n")
                    f.write(f"  ntyp = {len(struct_info['elements'])}\n")
                
                params = self.qe_params.get(section, {})
                for key, value in params.items():
                    f.write(f"  {key} = {value}\n")
                f.write("/")
                f.write("\n\n")
            
            # ATOMIC_SPECIES
            f.write("ATOMIC_SPECIES\n")
            pseudo_map = self.config.get("pseudo_map", {})
            for el in struct_info["elements"]:
                info = pseudo_map.get(el, {"mass": 1.0, "pseudo": f"{el}.UPF"})
                f.write(f"  {el:3} {info.get('mass', 1.0):8.4f} {info.get('pseudo', f'{el}.UPF')}\n")
            f.write("\n")
            
            # ATOMIC_POSITIONS
            coord_type = struct_info["coord_type"]
            if coord_type.lower().startswith('d'):
                f.write("ATOMIC_POSITIONS crystal\n")
            else:
                f.write("ATOMIC_POSITIONS angstrom\n")
                
            atom_idx = 0
            for el_idx, el in enumerate(struct_info["elements"]):
                for _ in range(struct_info["counts"][el_idx]):
                    pos = struct_info["positions"][atom_idx]
                    f.write(f"  {el:3} {pos[0]:12.8f} {pos[1]:12.8f} {pos[2]:12.8f}\n")
                    atom_idx += 1
            f.write("\n")
            
            # K_POINTS
            kmesh_val = self.config.get("kmesh", 0.04)
            kpts = self.get_kpoints(struct_info["lattice"], kmesh_val)
            f.write("K_POINTS automatic\n")
            f.write(f"  {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0\n\n")
            
            # CELL_PARAMETERS
            f.write("CELL_PARAMETERS angstrom\n")
            for vec in struct_info["lattice"]:
                f.write(f"  {vec[0]:12.8f} {vec[1]:12.8f} {vec[2]:12.8f}\n")

    def copy_pseudos(self, elements):
        """拷贝赝势文件"""
        pseudo_dir_src = self.config.get("qe_pseudo_dir", "")
        if not pseudo_dir_src:
            print("提示: 未在配置文件中设置 'qe_pseudo_dir'，跳过赝势拷贝。")
            return
        
        pseudo_dir_src = os.path.expanduser(pseudo_dir_src)
        target_pseudo_dir = os.path.join(self.work_dir, "pseudo")
        os.makedirs(target_pseudo_dir, exist_ok=True)
        
        pseudo_map = self.config.get("pseudo_map", {})
        for el in elements:
            pseudo_file = pseudo_map.get(el, {}).get("pseudo", f"{el}.UPF")
            src_path = os.path.join(pseudo_dir_src, pseudo_file)
            if os.path.exists(src_path):
                shutil.copy2(src_path, os.path.join(target_pseudo_dir, pseudo_file))
                print(f"已拷贝赝势: {pseudo_file}")
            else:
                print(f"警告: 找不到赝势文件 {src_path}")

    def create_run_script(self):
        """创建运行脚本"""
        qe_config = self.config.get("qe", {})
        pw_path = qe_config.get("executable_path", "mpirun -np 4 pw.x")
        slurm_header = self.config.get("slurm", {}).get("header", "#!/bin/bash")

        run_script_path = os.path.join(self.work_dir, "run_qe.sh")
        with open(run_script_path, "w") as f:
            f.write(slurm_header.strip() + "\n\n")
            f.write(f"echo 'Job started at' `date` \n")
            f.write(f"{pw_path} < qe.in > qe.out\n")
            f.write(f"echo 'Job finished at' `date` \n")
        os.chmod(run_script_path, 0o755)
        print(f"生成运行脚本: {run_script_path}")

    def run(self, run_calc=False):
        try:
            os.makedirs(self.work_dir, exist_ok=True)

            target_poscar = os.path.join(self.work_dir, "POSCAR")
            if self.struct_file:
                shutil.copy2(self.struct_file, target_poscar)
                print(f"已将 {self.struct_file} 拷贝为 {target_poscar}")
            elif os.path.exists("POSCAR"):
                shutil.copy2("POSCAR", target_poscar)
                print(f"已将当前目录下的 POSCAR 拷贝为 {target_poscar}")
            else:
                raise FileNotFoundError("找不到 POSCAR 文件，请使用 -i 指定或在当前目录准备。")

            struct_info = self.parse_poscar(target_poscar)
            self.generate_qe_input(struct_info)
            self.copy_pseudos(struct_info["elements"])
            self.create_run_script()
            
            print(f"\n所有 QE 输入文件已在 {self.work_dir} 目录中准备就绪！")

            if run_calc:
                # 执行逻辑
                pass
                
        except Exception as e:
            print(f"错误: {e}")

def main():
    parser = argparse.ArgumentParser(description="Quantum ESPRESSO SCF Setup Script")
    parser.add_argument("-i", "--input", help="输入结构文件 (POSCAR 格式)")
    parser.add_argument("-c", "--config", default="input.toml", help="配置文件路径")
    parser.add_argument("--run", action="store_true", help="直接执行计算")
    
    args = parser.parse_args()
    setup = QESetup(config_file=args.config, struct_file=args.input)
    setup.run(run_calc=args.run)

if __name__ == "__main__":
    main()
