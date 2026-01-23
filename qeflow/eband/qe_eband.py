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
        if count >= 1:
            formula += str(count)
    return formula


def resolve_qe_executable(bin_path, exe_name):
    if not bin_path:
        return exe_name
    bin_path = os.path.expanduser(str(bin_path))
    return os.path.join(bin_path, exe_name)


def parse_kpath(path):
    points = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 4:
                try:
                    kx, ky, kz = map(float, parts[:3])
                except ValueError:
                    continue
                label = parts[3]
                points.append(([kx, ky, kz], label))
            else:
                continue

    if len(points) % 2 != 0:
        raise ValueError(f"KPATH.in 格式异常，点数为奇数: {len(points)}")

    path_points = []
    for i in range(0, len(points), 2):
        (k1, l1), (k2, l2) = points[i], points[i + 1]
        if not path_points or path_points[-1][1] != l1 or path_points[-1][0] != k1:
            path_points.append((k1, l1))
        path_points.append((k2, l2))

    return path_points


class QEEBandSetup:
    def __init__(self, config_file="input.toml", struct_file=None):
        self.work_dir = "qe_eband"
        self.scf_dir = "qe_scf"

        if not os.path.exists(config_file):
            raise FileNotFoundError(f"找不到配置文件: {config_file}")

        if "tomllib" in globals():
            with open(config_file, "rb") as f:
                self.config = tomllib.load(f)
        else:
            self.config = toml.load(config_file)

        self.struct_file = struct_file
        self.k_points_config = self.config.get("k_points", {})
        self.kpath_points = int(self.k_points_config.get("kpath_points", 100))
        self.kpath_file = self.k_points_config.get("kpath_file", "KPATH.in")

        self.qe_params = {
            "CONTROL": {
                "calculation": "'bands'",
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
                "occupations": "'tetrahedra'",
            },
            "ELECTRONS": {
                "conv_thr": 1.0e-8,
                "mixing_beta": 0.7,
                "electron_maxstep": 400,
            },
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

    def run_vaspkit(self):
        print("正在调用 vaspkit 生成能带路径 (303)...")
        try:
            subprocess.run(
                ["vaspkit"],
                input="3\n303\n",
                text=True,
                check=True,
                cwd=self.work_dir,
            )
        except Exception as e:
            raise RuntimeError(f"vaspkit 运行失败: {e}")

    def write_kpath_labels(self, path_points, lattice):
        bohr_to_ang = 0.52917721092
        lattice_bohr = lattice / bohr_to_ang
        alat_bohr = float(np.linalg.norm(lattice_bohr[0]))
        # QE 的 scf.out 中倒格子以 2π/alat 为单位；此处按同一单位自算
        recip_lat = np.linalg.inv(lattice_bohr).T * alat_bohr
        total_dist = 0.0
        last_cart = None
        label_path = os.path.join(self.work_dir, "qe_k_lable.dat")
        with open(label_path, "w") as f:
            f.write("# label dist kx ky kz\n")
            for coords, label in path_points:
                cart = np.dot(coords, recip_lat)
                if last_cart is not None:
                    total_dist += float(np.linalg.norm(cart - last_cart))
                f.write(
                    f"{label:8} {total_dist:12.6f} {coords[0]:12.8f} {coords[1]:12.8f} {coords[2]:12.8f}\n"
                )
                last_cart = cart
        print(f"已生成高对称点坐标: {label_path}")

    def generate_qe_input(self, struct_info):
        print("正在生成 eband.in ...")

        formula = get_formula(struct_info["elements"], struct_info["counts"])
        self.qe_params["CONTROL"]["prefix"] = f"'{formula}'"

        if "qe_params" in self.config:
            for section, params in self.config["qe_params"].items():
                section_upper = section.upper()
                if section_upper in self.qe_params:
                    self.qe_params[section_upper].update(params)
                else:
                    self.qe_params[section_upper] = params

        if "eband_params" in self.config:
            for section, params in self.config["eband_params"].items():
                section_upper = section.upper()
                if section_upper in self.qe_params:
                    self.qe_params[section_upper].update(params)
                else:
                    self.qe_params[section_upper] = params

        if "qe_params" in self.config and "CONTROL" in self.config["qe_params"]:
            if "prefix" not in self.config["qe_params"]["CONTROL"]:
                self.qe_params["CONTROL"]["prefix"] = f"'{formula}'"
        else:
            self.qe_params["CONTROL"]["prefix"] = f"'{formula}'"

        if "CONTROL" in self.qe_params:
            pdir = self.qe_params["CONTROL"].get("pseudo_dir", "").strip("'").strip('"')
            if pdir.startswith("~"):
                pdir = os.path.expanduser(pdir)
                self.qe_params["CONTROL"]["pseudo_dir"] = f"'{pdir}'"

        qe_input_path = os.path.join(self.work_dir, "eband.in")
        with open(qe_input_path, "w") as f:
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

            f.write("ATOMIC_SPECIES\n")
            pseudo_map = self.config.get("pseudo_map", {})
            for el in struct_info["elements"]:
                info = pseudo_map.get(el)
                if not info or "pseudo" not in info:
                    raise KeyError(f"pseudo_map 缺少元素 {el} 的 pseudo 配置")
                mass = info.get("mass")
                if mass is None:
                    mass = ATOMIC_MASSES.get(el, 1.0)

                pseudo_file = str(info["pseudo"]).strip()
                f.write(f"  {el:3} {mass:8.4f} {pseudo_file}\n")
            f.write("\n")

            coord_type = struct_info["coord_type"]
            if coord_type.lower().startswith("d"):
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

            f.write("CELL_PARAMETERS {angstrom}\n")
            for vec in struct_info["lattice"]:
                f.write(f"  {vec[0]:12.8f} {vec[1]:12.8f} {vec[2]:12.8f}\n")
            f.write("\n")

            kpath_path = self.kpath_file
            if not os.path.isabs(kpath_path):
                kpath_path = os.path.join(self.work_dir, kpath_path)

            if not os.path.exists(kpath_path):
                raise FileNotFoundError(f"找不到 KPATH.in: {kpath_path}")

            path_points = parse_kpath(kpath_path)
            if not path_points:
                raise ValueError("KPATH.in 未解析到任何高对称点。")

            self.write_kpath_labels(path_points, struct_info["lattice"])

            f.write("K_POINTS crystal_b\n")
            f.write(f"{len(path_points)}\n")
            for coords, label in path_points:
                f.write(
                    f"  {coords[0]:12.8f} {coords[1]:12.8f} {coords[2]:12.8f} {self.kpath_points:6d} ! {label}\n"
                )

        return formula

    def generate_banddata_inputs(self, prefix):
        outdir_val = self.qe_params["CONTROL"].get("outdir", "'./out'").strip("'").strip('"')
        eband_cfg = self.config.get("eband", {})
        bands_filband = eband_cfg.get("filband", "elebanddata.dat")
        proj_filproj = eband_cfg.get("filproj", "elebandprojdata")
        lsym = eband_cfg.get("lsym", False)

        bands_input_path = os.path.join(self.work_dir, "elebanddata.in")
        with open(bands_input_path, "w") as f:
            f.write("&BANDS\n")
            f.write(f" prefix = '{prefix}',\n")
            f.write(f" outdir = '{outdir_val}',\n")
            f.write(f" filband = '{bands_filband}',\n")
            f.write(" lp = .true.\n")
            f.write("/\n")

        proj_input_path = os.path.join(self.work_dir, "elebandprojdata.in")
        with open(proj_input_path, "w") as f:
            f.write("&projwfc\n")
            f.write(f" prefix = '{prefix}',\n")
            f.write(f" outdir = '{outdir_val}',\n")
            f.write(f" lsym = .{str(lsym).lower()}.,\n")
            f.write(f" filproj = '{proj_filproj}'\n")
            f.write("/\n")

    def create_run_script(self):
        qe_config = self.config.get("qe", {})
        qe_bin = qe_config.get("executable_path", "")
        pw_path = resolve_qe_executable(qe_bin, "pw.x")
        bands_cfg = self.config.get("bands", {})
        bands_exec = resolve_qe_executable(bands_cfg.get("executable_path", qe_bin), "bands.x")
        bands_args = bands_cfg.get("args", "-pd .true.")
        proj_cfg = self.config.get("projwfc", {})
        proj_exec = resolve_qe_executable(proj_cfg.get("executable_path", qe_bin), "projwfc.x")
        proj_args = proj_cfg.get("args", "-pd .true.")
        slurm_header = self.config.get("slurm", {}).get("header", "#!/bin/bash")

        run_script_path = os.path.join(self.work_dir, "run_qe.sh")
        with open(run_script_path, "w") as f:
            f.write(slurm_header.strip() + "\n\n")
            f.write("echo 'Job started at' `date` \n")
            f.write(f"{pw_path} < eband.in > eband.out\n")
            if bands_args:
                f.write(f"{bands_exec} {bands_args} < elebanddata.in > elebanddata.out\n")
            else:
                f.write(f"{bands_exec} < elebanddata.in > elebanddata.out\n")
            if proj_args:
                f.write(f"{proj_exec} {proj_args} < elebandprojdata.in > elebandprojdata.out\n")
            else:
                f.write(f"{proj_exec} < elebandprojdata.in > elebandprojdata.out\n")
            f.write("echo 'Job finished at' `date` \n")
        os.chmod(run_script_path, 0o755)
        print(f"生成运行脚本: {run_script_path}")

    def execute(self):
        qe_config = self.config.get("qe", {})
        qe_bin = qe_config.get("executable_path", "")
        pw_path = resolve_qe_executable(qe_bin, "pw.x")
        bands_cfg = self.config.get("bands", {})
        bands_exec = resolve_qe_executable(bands_cfg.get("executable_path", qe_bin), "bands.x")
        bands_args = bands_cfg.get("args", "-pd .true.")
        proj_cfg = self.config.get("projwfc", {})
        proj_exec = resolve_qe_executable(proj_cfg.get("executable_path", qe_bin), "projwfc.x")
        proj_args = proj_cfg.get("args", "-pd .true.")

        original_dir = os.getcwd()
        try:
            work_dir_name = os.path.basename(os.path.normpath(self.work_dir))
            os.chdir(self.work_dir)
            pw_cmd = f"{pw_path} < eband.in > eband.out"
            print(f"\n正在进入 {work_dir_name} 目录执行: {pw_cmd}")
            subprocess.run(pw_cmd, shell=True, check=True)
            bands_cmd = f"{bands_exec} {bands_args}".strip() if bands_args else bands_exec
            proj_cmd = f"{proj_exec} {proj_args}".strip() if proj_args else proj_exec
            bands_full_cmd = f"{bands_cmd} < elebanddata.in > elebanddata.out"
            print(f"正在进入 {work_dir_name} 目录执行: {bands_full_cmd}")
            subprocess.run(bands_full_cmd, shell=True, check=True)
            proj_full_cmd = f"{proj_cmd} < elebandprojdata.in > elebandprojdata.out"
            print(f"正在进入 {work_dir_name} 目录执行: {proj_full_cmd}")
            subprocess.run(proj_full_cmd, shell=True, check=True)
            print(f"QE eband 计算完成，输出已保存至 {self.work_dir}")
        except subprocess.CalledProcessError as e:
            print(f"执行 QE eband 出错: {e}")
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
            self.qe_params["CONTROL"]["prefix"] = f"'{formula}'"

            if "qe_params" in self.config and "CONTROL" in self.config["qe_params"]:
                self.qe_params["CONTROL"].update(self.config["qe_params"]["CONTROL"])
            if "eband_params" in self.config and "CONTROL" in self.config["eband_params"]:
                self.qe_params["CONTROL"].update(self.config["eband_params"]["CONTROL"])

            outdir_val = self.qe_params["CONTROL"]["outdir"].strip("'").strip('"')
            prefix_val = self.qe_params["CONTROL"]["prefix"].strip("'").strip('"')

            save_dir_name = f"{prefix_val}.save"
            src_save_path = os.path.join(self.scf_dir, outdir_val, save_dir_name)
            dst_save_path = os.path.join(self.work_dir, outdir_val, save_dir_name)

            if os.path.exists(src_save_path):
                os.makedirs(dst_save_path, exist_ok=True)
                needed_files = ["data-file-schema.xml", "charge-density.dat", "paw.txt"]
                for filename in needed_files:
                    src_file = os.path.join(src_save_path, filename)
                    dst_file = os.path.join(dst_save_path, filename)
                    if os.path.exists(src_file):
                        shutil.copy2(src_file, dst_file)
                        print(f"已拷贝 SCF 数据: {filename}")
                print(f"已从 {src_save_path} 准备好 eband 所需的关键数据。")
            else:
                print(f"警告: 找不到 SCF 的数据目录 {src_save_path}，计算可能会失败。")

            self.run_vaspkit()

            prefix = self.generate_qe_input(struct_info)
            self.generate_banddata_inputs(prefix)
            self.create_run_script()

            if outdir_val:
                full_outdir = os.path.join(self.work_dir, outdir_val)
                os.makedirs(full_outdir, exist_ok=True)
                print(f"确保 outdir 存在: {full_outdir}")

            print(f"\n所有 QE eband 输入文件及数据已在 {self.work_dir} 目录中准备就绪！")

        except Exception as e:
            print(f"准备 eband 环境出错: {e}")
            raise


def main():
    parser = argparse.ArgumentParser(description="Quantum ESPRESSO eband Setup Script")
    parser.add_argument("-i", "--input", help="输入结构文件 (POSCAR 格式)")
    parser.add_argument("-c", "--config", default="input.toml", help="配置文件路径")
    parser.add_argument("--run", action="store_true", help="生成文件后直接执行计算")

    args = parser.parse_args()
    setup_obj = QEEBandSetup(config_file=args.config, struct_file=args.input)

    setup_obj.setup()

    if args.run:
        setup_obj.execute()


if __name__ == "__main__":
    main()
