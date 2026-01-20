#!/usr/bin/env python3
import os
import sys
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
        sys.exit(1)

ATOMIC_MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
    "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.63, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Tc": 98, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29,
    "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Th": 232.04, "Pa": 231.04, "U": 238.03
}

def get_formula(elements, counts):
    formula = ""
    for el, count in zip(elements, counts):
        formula += el
        if count > 1:
            formula += str(count)
    return formula

class QESetup:
    def __init__(self, config_file="input.toml", struct_file=None):
        self.work_dir = "qe_scf"
        self.k_points_config = {}
        # 1. 加载配置文件
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"找不到配置文件: {config_file}")
        
        if 'tomllib' in globals():
            with open(config_file, "rb") as f:
                self.config = tomllib.load(f)
        else:
            self.config = toml.load(config_file)
        
        self.struct_file = struct_file
        self.k_points_config = self.config.get("k_points", {})
        
        # 2. 默认 QE 参数模板 (Namelists)
        self.qe_params = {
            "CONTROL": {
                "calculation": "'scf'",
                "restart_mode": "'from_scratch'",
                "prefix": "'qe'",
                "outdir": "'./out'",
                "pseudo_dir": "'./pseudo'",
                "verbosity": "'high'",
                "tprnfor": ".true.",
                "tstress": ".true.",
                "forc_conv_thr": "1.0e-6",
                "etot_conv_thr": "1.0e-7",
            },
            "SYSTEM": {
                "ibrav": 0,
                "ecutwfc": 50.0,
                "ecutrho": 400.0,
                "occupations": "'smearing'",
                "smearing": "'gaussian'",
                "degauss": 0.01,
            },
            "ELECTRONS": {
                "conv_thr": 1.0e-8,
                "mixing_beta": 0.7,
                "electron_maxstep": 400,
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
        
        pos_array = np.array(positions, dtype=float)
        # 如果是笛卡尔坐标，也需要乘以缩放因子
        if not coord_line.lower().startswith('d'):
            pos_array *= scale
            
        return {
            "lattice": np.array(lattice),
            "elements": elements,
            "counts": counts,
            "coord_type": coord_line,
            "positions": pos_array
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
        print(f"正在生成 scf.in...")
        
        # 1. 自动设置 prefix (由化学配比决定)
        formula = get_formula(struct_info["elements"], struct_info["counts"])
        self.qe_params["CONTROL"]["prefix"] = f"'{formula}'"

        # 2. 合并配置文件中的参数 (覆盖默认值，但 prefix 已由上面逻辑确定)
        if "qe_params" in self.config:
            for section, params in self.config["qe_params"].items():
                section_upper = section.upper()
                if section_upper in self.qe_params:
                    # 如果配置文件里也写了 prefix，这里会覆盖上面的自动设置
                    # 但根据需求，我们应该优先使用自动生成的，或者如果用户没写才用自动生成的。
                    # 用户说“不用在input.toml里面手动设置”，意味着我们应该自动生成。
                    # 如果用户还是写了，我们决定是否覆盖。这里我们选择如果用户写了则遵循用户，没写用自动。
                    # 修正：用户明确说“不用手动设置”，所以我们总是自动生成。
                    self.qe_params[section_upper].update(params)
                else:
                    self.qe_params[section_upper] = params
        
        # 强制使用自动生成的 prefix (如果用户在 config 里也写了，这里再次确保覆盖)
        if "qe_params" in self.config and "CONTROL" in self.config["qe_params"]:
             if "prefix" not in self.config["qe_params"]["CONTROL"]:
                 self.qe_params["CONTROL"]["prefix"] = f"'{formula}'"
        else:
             self.qe_params["CONTROL"]["prefix"] = f"'{formula}'"

        # 3. 处理 pseudo_dir 中的 ~ 符号 (QE 不支持 ~)
        if "CONTROL" in self.qe_params:
            pdir = self.qe_params["CONTROL"].get("pseudo_dir", "").strip("'").strip('"')
            if pdir.startswith("~"):
                pdir = os.path.expanduser(pdir)
                self.qe_params["CONTROL"]["pseudo_dir"] = f"'{pdir}'"

        qe_input_path = os.path.join(self.work_dir, "scf.in")
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
                info = pseudo_map.get(el, {})
                # 自动获取质量
                mass = info.get("mass")
                if mass is None:
                    mass = ATOMIC_MASSES.get(el, 1.0)
                
                # 剔除伪势文件名中的多余空格
                pseudo_file = info.get('pseudo', f'{el}.UPF').strip()
                f.write(f"  {el:3} {mass:8.4f} {pseudo_file}\n")
            f.write("\n")
            
            # ATOMIC_POSITIONS
            coord_type = struct_info["coord_type"]
            if coord_type.lower().startswith('d'):
                f.write("ATOMIC_POSITIONS {crystal}\n")
            else:
                f.write("ATOMIC_POSITIONS {angstrom}\n")
                
            atom_idx = 0
            for el_idx, el in enumerate(struct_info["elements"]):
                for _ in range(struct_info["counts"][el_idx]):
                    pos = struct_info["positions"][atom_idx]
                    f.write(f"  {el:3} {pos[0]:12.8f} {pos[1]:12.8f} {pos[2]:12.8f}\n")
                    atom_idx += 1
            f.write("\n")
            
            # CELL_PARAMETERS
            f.write("CELL_PARAMETERS {angstrom}\n")
            for vec in struct_info["lattice"]:
                f.write(f"  {vec[0]:12.8f} {vec[1]:12.8f} {vec[2]:12.8f}\n")
            f.write("\n")
            
            # K_POINTS：支持显式 kpoints 或按 kmesh 生成自动网格
            kpoints_list = self.k_points_config.get("kpoints")
            kmesh_val = self.k_points_config.get("kmesh", 0.04)

            if isinstance(kpoints_list, list) and len(kpoints_list) == 3:
                kpts = kpoints_list
            else:
                kpts = self.get_kpoints(struct_info["lattice"], kmesh_val)

            f.write("K_POINTS automatic\n")
            f.write(f"  {int(kpts[0])} {int(kpts[1])} {int(kpts[2])} 0 0 0\n")

    def create_run_script(self):
        """创建运行脚本"""
        qe_config = self.config.get("qe", {})
        pw_path = qe_config.get("executable_path", "mpirun -np 4 pw.x")
        slurm_header = self.config.get("slurm", {}).get("header", "#!/bin/bash")

        run_script_path = os.path.join(self.work_dir, "run_qe.sh")
        with open(run_script_path, "w") as f:
            f.write(slurm_header.strip() + "\n\n")
            f.write(f"echo 'Job started at' `date` \n")
            f.write(f"{pw_path} < scf.in > scf.out\n")
            f.write(f"echo 'Job finished at' `date` \n")
        os.chmod(run_script_path, 0o755)
        print(f"生成运行脚本: {run_script_path}")

    def execute(self):
        """进入工作目录并执行 Quantum ESPRESSO 计算"""
        qe_config = self.config.get("qe", {})
        pw_path = qe_config.get("executable_path", "mpirun -np 4 pw.x")
        
        original_dir = os.getcwd()
        try:
            print(f"\n正在进入 {self.work_dir} 目录执行: {pw_path}")
            os.chdir(self.work_dir)
            # 使用 shell=True 因为 pw_path 可能包含 mpirun 等复杂命令
            # QE 通常通过 < scf.in 读取输入
            cmd = f"{pw_path} < scf.in > scf.out"
            subprocess.run(cmd, shell=True, check=True)
            print(f"QE 计算完成，输出已保存至 {self.work_dir}/scf.out")
        except subprocess.CalledProcessError as e:
            print(f"执行 QE 出错: {e}")
        finally:
            os.chdir(original_dir)
            print(f"已返回目录: {original_dir}")

    def setup(self):
        """准备所有输入文件和计算环境"""
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
            self.create_run_script()
            
            # 确保 outdir 目录存在，否则 QE 会报错
            outdir_val = self.qe_params["CONTROL"].get("outdir", "").strip("'").strip('"')
            if outdir_val:
                full_outdir = os.path.join(self.work_dir, outdir_val)
                os.makedirs(full_outdir, exist_ok=True)
                print(f"确保 outdir 存在: {full_outdir}")

            print(f"\n所有 QE 输入文件已在 {self.work_dir} 目录中准备就绪！")
                
        except Exception as e:
            print(f"准备环境出错: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(description="Quantum ESPRESSO SCF Setup Script")
    parser.add_argument("-i", "--input", help="输入结构文件 (POSCAR 格式)")
    parser.add_argument("-c", "--config", default="input.toml", help="配置文件路径")
    parser.add_argument("--run", action="store_true", help="生成文件后直接执行计算")
    
    args = parser.parse_args()
    setup_obj = QESetup(config_file=args.config, struct_file=args.input)
    
    # 1. 环境准备
    setup_obj.setup()
    
    # 2. 如果指定了 --run，则执行计算
    if args.run:
        setup_obj.execute()

if __name__ == "__main__":
    main()
