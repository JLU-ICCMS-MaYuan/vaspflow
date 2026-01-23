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


def resolve_qe_executable(bin_path, exe_name):
    if not bin_path:
        return exe_name
    bin_path = os.path.expanduser(str(bin_path))
    return os.path.join(bin_path, exe_name)


def get_formula(elements, counts):
    formula = ""
    for el, count in zip(elements, counts):
        formula += el
        if count >= 1:
            formula += str(count)
    return formula


class QEEleDosSetup:
    def __init__(self, config_file="input.toml", struct_file=None):
        self.work_dir = "qe_eledos"
        self.scf_dir = "qe_scf"

        if not os.path.exists(config_file):
            raise FileNotFoundError(f"找不到配置文件: {config_file}")

        if "tomllib" in globals():
            with open(config_file, "rb") as f:
                self.config = tomllib.load(f)
        else:
            self.config = toml.load(config_file)

        self.struct_file = struct_file
        self.qe_params = {
            "CONTROL": {
                "outdir": "'./tmp'",
            }
        }

    def parse_poscar(self, poscar_path):
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
        if coord_line.lower().startswith("s"):
            coord_line = lines[8].strip()
            start_line = 9
        else:
            start_line = 8

        positions = []
        total_atoms = sum(counts)
        for i in range(total_atoms):
            positions.append(lines[start_line + i].split()[:3])

        pos_array = np.array(positions, dtype=float)
        if not coord_line.lower().startswith("d"):
            pos_array *= scale

        return {
            "lattice": np.array(lattice),
            "elements": elements,
            "counts": counts,
            "coord_type": coord_line,
            "positions": pos_array,
        }

    def generate_dos_inputs(self, prefix):
        eledos_cfg = self.config.get("eledos", {})
        outdir_val = self.qe_params["CONTROL"].get("outdir", "'./tmp'").strip("'").strip('"')
        outdir_cfg = eledos_cfg.get("outdir")
        if outdir_cfg:
            outdir_val = str(outdir_cfg).strip("'").strip('"')

        fildos = eledos_cfg.get("fildos", f"{prefix}.tdos")
        filpdos = eledos_cfg.get("filpdos", f"{prefix}.pdos")
        filproj = eledos_cfg.get("filproj", f"{prefix}.proj")

        delta_e = eledos_cfg.get("DeltaE", 0.01)
        emin = eledos_cfg.get("emin", -10)
        emax = eledos_cfg.get("emax", 30)

        ngauss = eledos_cfg.get("ngauss", 0)
        degauss = eledos_cfg.get("degauss", 0.0)

        eletdos_path = os.path.join(self.work_dir, "eletdos.in")
        with open(eletdos_path, "w") as f:
            f.write("&dos\n")
            f.write(f"   prefix = '{prefix}',\n")
            f.write(f"   outdir = '{outdir_val}',\n")
            f.write(f"   fildos = '{fildos}',\n")
            f.write(f"   DeltaE = {delta_e},\n")
            f.write(f"   emin = {emin},\n")
            f.write(f"   emax = {emax}\n")
            f.write("/\n")

        elepdos_path = os.path.join(self.work_dir, "elepdos.in")
        with open(elepdos_path, "w") as f:
            f.write("&projwfc\n")
            f.write(f"   prefix = '{prefix}',\n")
            f.write(f"   outdir = '{outdir_val}',\n")
            f.write(f"   filpdos = '{filpdos}',\n")
            f.write(f"   filproj = '{filproj}',\n")
            f.write(f"   ngauss = {ngauss},\n")
            f.write(f"   DeltaE = {delta_e},\n")
            f.write(f"   degauss = {degauss},\n")
            f.write(f"   emin = {emin},\n")
            f.write(f"   emax = {emax}\n")
            f.write("/\n")

    def create_run_script(self):
        qe_config = self.config.get("qe", {})
        qe_bin = qe_config.get("executable_path", "")
        dos_cfg = self.config.get("dos", {})
        dos_exec = resolve_qe_executable(dos_cfg.get("executable_path", qe_bin), "dos.x")
        dos_args = dos_cfg.get("args", "-pd .true.")
        proj_cfg = self.config.get("projwfc", {})
        proj_exec = resolve_qe_executable(proj_cfg.get("executable_path", qe_bin), "projwfc.x")
        proj_args = proj_cfg.get("args", "-pd .true.")
        slurm_header = self.config.get("slurm", {}).get("header", "#!/bin/bash")

        run_script_path = os.path.join(self.work_dir, "run_qe.sh")
        with open(run_script_path, "w") as f:
            f.write(slurm_header.strip() + "\n\n")
            f.write("echo 'Job started at' `date` \n")
            if dos_args:
                f.write(f"{dos_exec} {dos_args} < eletdos.in > eletdos.out\n")
            else:
                f.write(f"{dos_exec} < eletdos.in > eletdos.out\n")
            if proj_args:
                f.write(f"{proj_exec} {proj_args} < elepdos.in > elepdos.out\n")
            else:
                f.write(f"{proj_exec} < elepdos.in > elepdos.out\n")
            f.write("echo 'Job finished at' `date` \n")
        os.chmod(run_script_path, 0o755)
        print(f"生成运行脚本: {run_script_path}")

    def execute(self):
        qe_config = self.config.get("qe", {})
        qe_bin = qe_config.get("executable_path", "")
        dos_cfg = self.config.get("dos", {})
        dos_exec = resolve_qe_executable(dos_cfg.get("executable_path", qe_bin), "dos.x")
        dos_args = dos_cfg.get("args", "-pd .true.")
        proj_cfg = self.config.get("projwfc", {})
        proj_exec = resolve_qe_executable(proj_cfg.get("executable_path", qe_bin), "projwfc.x")
        proj_args = proj_cfg.get("args", "-pd .true.")

        original_dir = os.getcwd()
        try:
            work_dir_name = os.path.basename(os.path.normpath(self.work_dir))
            os.chdir(self.work_dir)
            dos_cmd = f"{dos_exec} {dos_args}".strip() if dos_args else dos_exec
            dos_full_cmd = f"{dos_cmd} < eletdos.in > eletdos.out"
            print(f"\n正在进入 {work_dir_name} 目录执行: {dos_full_cmd}")
            subprocess.run(dos_full_cmd, shell=True, check=True)
            proj_cmd = f"{proj_exec} {proj_args}".strip() if proj_args else proj_exec
            proj_full_cmd = f"{proj_cmd} < elepdos.in > elepdos.out"
            print(f"正在进入 {work_dir_name} 目录执行: {proj_full_cmd}")
            subprocess.run(proj_full_cmd, shell=True, check=True)
            print(f"QE eledos 计算完成，输出已保存至 {self.work_dir}")
        except subprocess.CalledProcessError as e:
            print(f"执行 QE eledos 出错: {e}")
        finally:
            os.chdir(original_dir)
            print(f"已返回目录: {original_dir}")

    def setup(self):
        try:
            os.makedirs(self.work_dir, exist_ok=True)

            if self.struct_file:
                struct_path = self.struct_file
            elif os.path.exists("POSCAR"):
                struct_path = "POSCAR"
            else:
                raise FileNotFoundError("找不到结构文件。请在当前目录准备 POSCAR 或使用 -i 指定。")

            target_poscar = os.path.join(self.work_dir, "POSCAR")
            if os.path.abspath(struct_path) != os.path.abspath(target_poscar):
                shutil.copy2(struct_path, target_poscar)

            struct_info = self.parse_poscar(struct_path)
            formula = get_formula(struct_info["elements"], struct_info["counts"])

            if "qe_params" in self.config and "CONTROL" in self.config["qe_params"]:
                self.qe_params["CONTROL"].update(self.config["qe_params"]["CONTROL"])

            if "prefix" not in self.qe_params["CONTROL"] or self.qe_params["CONTROL"]["prefix"] in ("''", '""'):
                self.qe_params["CONTROL"]["prefix"] = f"'{formula}'"

            prefix_val = self.qe_params["CONTROL"]["prefix"].strip("'").strip('"')
            outdir_val = self.qe_params["CONTROL"]["outdir"].strip("'").strip('"')

            src_save_path = os.path.join(self.scf_dir, outdir_val, f"{prefix_val}.save")
            dst_save_path = os.path.join(self.work_dir, outdir_val, f"{prefix_val}.save")
            if os.path.exists(src_save_path):
                os.makedirs(os.path.join(self.work_dir, outdir_val), exist_ok=True)
                shutil.copytree(src_save_path, dst_save_path, dirs_exist_ok=True)
                print(f"已从 {src_save_path} 拷贝 SCF 数据到 {dst_save_path}")
            else:
                print(f"警告: 找不到 SCF 数据目录 {src_save_path}，DOS 计算可能会失败。")

            self.generate_dos_inputs(prefix_val)
            self.create_run_script()

            if outdir_val:
                full_outdir = os.path.join(self.work_dir, outdir_val)
                os.makedirs(full_outdir, exist_ok=True)
                print(f"确保 outdir 存在: {full_outdir}")

            print(f"\n所有 QE eledos 输入文件及数据已在 {self.work_dir} 目录中准备就绪！")

        except Exception as e:
            print(f"准备 eledos 环境出错: {e}")
            raise


def main():
    parser = argparse.ArgumentParser(description="Quantum ESPRESSO eledos Setup Script")
    parser.add_argument("-i", "--input", help="输入结构文件 (POSCAR 格式)")
    parser.add_argument("-c", "--config", default="input.toml", help="配置文件路径")
    parser.add_argument("--run", action="store_true", help="生成文件后直接执行计算")

    args = parser.parse_args()
    setup_obj = QEEleDosSetup(config_file=args.config, struct_file=args.input)

    setup_obj.setup()

    if args.run:
        setup_obj.execute()


if __name__ == "__main__":
    main()
