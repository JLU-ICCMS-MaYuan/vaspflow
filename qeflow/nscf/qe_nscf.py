#!/usr/bin/env python3
import os
import shutil
import argparse
import subprocess
import sys
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

class QENSCFSetup:
    def __init__(self, config_file="input.toml", struct_file=None):
        self.work_dir = "qe_nscf"
        self.scf_dir = "qe_scf"
        # 1. 加载配置文件
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"找不到配置文件: {config_file}")
        
        if 'tomllib' in globals():
            with open(config_file, "rb") as f:
                self.config = tomllib.load(f)
        else:
            self.config = toml.load(config_file)
        
        self.struct_file = struct_file
        nscf_config = self.config.get("nscf", {})
        self.k_automatic = nscf_config.get("k_automatic", True)
        self.wan = nscf_config.get("wan", False)
        
        # 2. 默认 QE 参数模板 (Namelists)
        # 注意：NSCF 默认使用 tetrahedra 且通常需要 high verbosity
        self.qe_params = {
            "CONTROL": {
                "calculation": "'nscf'",
                "restart_mode": "'from_scratch'",
                "prefix": "'qe'",
                "outdir": "'./out'",
                "pseudo_dir": "'./pseudo'",
                "verbosity": "'high'",
                "tprnfor": ".true.",
                "tstress": ".true.",
            },
            "SYSTEM": {
                "ecutwfc": 50.0,
                "ecutrho": 400.0,
                "occupations": "'tetrahedra'",
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
            
        scale = float(lines[1].strip())
        lattice = []
        for i in range(2, 5):
            lattice.append([float(x) * scale for x in lines[i].split()])
            
        elements = lines[5].split()
        if elements[0].isdigit():
            elements = lines[0].split()
            if not elements:
                raise ValueError("无法在 POSCAR 中找到元素名称。")
        
        counts = [int(x) for x in lines[6].split()]
        
        coord_line = lines[7].strip()
        if coord_line.lower().startswith('s'):
            coord_line = lines[8].strip()
            start_line = 9
        else:
            start_line = 8
            
        positions = []
        total_atoms = sum(counts)
        for i in range(total_atoms):
            positions.append(lines[start_line + i].split()[:3])
        
        pos_array = np.array(positions, dtype=float)
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
        
        b = 2 * np.pi * np.linalg.inv(lattice).T
        b_lengths = np.linalg.norm(b, axis=1)
        
        # NSCF 默认使用更密的网格
        spacing = kmesh if isinstance(kmesh, float) else 0.02
        kpts = [int(max(1, np.ceil(l / (spacing * 2 * np.pi)))) for l in b_lengths]
        return kpts

    def get_kmesh_justlike_kmesh_pl(self, kpoints):
        """
        读取self.kpoints_dense参数, 将其传给n1, n2, n3, 再将n1, n2, n3转化为相应的倒空间的均匀网格点坐标
        """
        # 获取输入的 n1, n2, n3
        n1, n2, n3 = kpoints

        # 参数检查：确保 n1, n2, n3 都大于 0
        if n1 <= 0:
            print("n1 must be > 0")
            sys.exit()
        if n2 <= 0:
            print("n2 must be > 0")
            sys.exit()
        if n3 <= 0:
            print("n3 must be > 0")
            sys.exit()

        # 计算总的 k 点数量
        totpts = n1 * n2 * n3

        kpoints_coords = []
        if not self.wan: # 前三列写k点倒空间分数坐标，第四列写其权重
            print("K_POINTS crystal")
            print(totpts)
            for x in range(n1):
                for y in range(n2):
                    for z in range(n3):
                        # 格式化输出 k 点信息
                        # print(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}{1/totpts:14.6e}")
                        kpoints_coords.append(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}{1/totpts:14.6e}")
        else:  # 只写前三列写k点倒空间分数坐标
            for x in range(n1):
                for y in range(n2):
                    for z in range(n3):
                        # 格式化输出 k 点信息（没有权重）
                        # print(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}")
                        kpoints_coords.append(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}")
        return kpoints_coords, totpts

    def generate_qe_input(self, struct_info):
        """生成 QE pw.x 输入文件"""
        print(f"正在生成 nscf.in...")
        
        # 合并配置文件中的参数
        if "qe_params" in self.config:
            for section, params in self.config["qe_params"].items():
                section_upper = section.upper()
                if section_upper in self.qe_params:
                    self.qe_params[section_upper].update(params)
                else:
                    self.qe_params[section_upper] = params

        # 合并 NSCF 特有参数覆盖
        if "nscf_params" in self.config:
            for section, params in self.config["nscf_params"].items():
                section_upper = section.upper()
                if section_upper in self.qe_params:
                    self.qe_params[section_upper].update(params)
                else:
                    self.qe_params[section_upper] = params

        qe_input_path = os.path.join(self.work_dir, "nscf.in")
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
            
            # K_POINTS
            nscf_config = self.config.get("nscf", {})
            kpoints_dense = nscf_config.get("kpoints_dense")
            if kpoints_dense is None:
                kmesh_val = nscf_config.get("kmesh", 0.02)
                kpoints_dense = self.get_kpoints(struct_info["lattice"], kmesh_val)
            if self.k_automatic:
                f.write("K_POINTS automatic\n")
                f.write(f"  {kpoints_dense[0]} {kpoints_dense[1]} {kpoints_dense[2]} 0 0 0\n")
            else:
                kpoints_coords, totpts = self.get_kmesh_justlike_kmesh_pl(kpoints_dense)
                f.write("K_POINTS crystal\n")
                f.write(f"{totpts}\n")
                for kinfo in kpoints_coords:
                    f.write(f" {kinfo}\n")

    def create_run_script(self):
        """创建运行脚本"""
        qe_config = self.config.get("qe", {})
        pw_path = qe_config.get("executable_path", "mpirun -np 4 pw.x")
        slurm_header = self.config.get("slurm", {}).get("header", "#!/bin/bash")

        run_script_path = os.path.join(self.work_dir, "run_qe.sh")
        with open(run_script_path, "w") as f:
            f.write(slurm_header.strip() + "\n\n")
            f.write(f"echo 'Job started at' `date` \n")
            f.write(f"{pw_path} < nscf.in > nscf.out\n")
            f.write(f"echo 'Job finished at' `date` \n")
        os.chmod(run_script_path, 0o755)
        print(f"生成运行脚本: {run_script_path}")

    def run(self, run_calc=False):
        try:
            os.makedirs(self.work_dir, exist_ok=True)

            # 1. 处理结构文件 (优先从 qe_scf 拷贝)
            target_poscar = os.path.join(self.work_dir, "POSCAR")
            scf_poscar = os.path.join(self.scf_dir, "POSCAR")
            
            if self.struct_file:
                shutil.copy2(self.struct_file, target_poscar)
                print(f"已将输入文件 {self.struct_file} 拷贝为 {target_poscar}")
            elif os.path.exists(scf_poscar):
                shutil.copy2(scf_poscar, target_poscar)
                print(f"已从 {self.scf_dir} 目录拷贝 POSCAR")
            elif os.path.exists("POSCAR"):
                shutil.copy2("POSCAR", target_poscar)
                print(f"已将当前目录下的 POSCAR 拷贝为 {target_poscar}")
            else:
                raise FileNotFoundError("找不到结构文件，且无法从 qe_scf 自动获取。")

                        # 2. 拷贝 SCF 生成的数据 (outdir 及其内部的 prefix.save)

                        # 自动获取 prefix 和 outdir，确保 NSCF 与 SCF 一致

                        outdir_val = self.qe_params["CONTROL"]["outdir"].strip("'").strip('"')

                        prefix_val = self.qe_params["CONTROL"]["prefix"].strip("'").strip('"')

                        

                        save_dir_name = f"{prefix_val}.save"

                        src_out_base = os.path.join(self.scf_dir, outdir_val)

                        src_save_path = os.path.join(src_out_base, save_dir_name)

                        

                        dst_out_base = os.path.join(self.work_dir, outdir_val)

                        dst_save_path = os.path.join(dst_out_base, save_dir_name)

                        

                        if os.path.exists(src_save_path):

                            # 创建目标 outdir 层级

                            os.makedirs(dst_out_base, exist_ok=True)

                            if os.path.exists(dst_save_path):

                                shutil.rmtree(dst_save_path)

                            

                            # 拷贝整个 .save 目录以包含 xml, dat, paw 等所有必要文件

                            shutil.copytree(src_save_path, dst_save_path)

                            print(f"已成功从 {src_save_path} 拷贝数据至 {dst_save_path}")

                            print(f"  - 包含: data-file-schema.xml, charge-density.dat, paw.txt 等")

                        else:

                            print(f"警告: 找不到 SCF 的数据目录 {src_save_path}")

                            print(f"      请确保 {self.scf_dir} 中已完成计算且 prefix='{prefix_val}', outdir='{outdir_val}'")

            

            struct_info = self.parse_poscar(target_poscar)
            self.generate_qe_input(struct_info)
            self.create_run_script()
            
            print(f"\n所有 QE NSCF 输入文件及数据已在 {self.work_dir} 目录中准备就绪！")

            if run_calc:
                pass
                
        except Exception as e:
            print(f"错误: {e}")

def main():
    parser = argparse.ArgumentParser(description="Quantum ESPRESSO NSCF Setup Script")
    parser.add_argument("-i", "--input", help="输入结构文件 (POSCAR 格式)")
    parser.add_argument("-c", "--config", default="input.toml", help="配置文件路径")
    parser.add_argument("--run", action="store_true", help="直接执行计算")
    
    args = parser.parse_args()
    setup = QENSCFSetup(config_file=args.config, struct_file=args.input)
    setup.run(run_calc=args.run)

if __name__ == "__main__":
    main()
