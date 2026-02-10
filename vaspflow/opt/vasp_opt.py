#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess

try:
    import tomllib  # Python 3.11+
except ImportError:
    try:
        import toml
    except ImportError:
        print("Error: 'toml' or 'tomllib' (Python 3.11+) is required.")
        print("Please install it using: pip install toml")
        raise SystemExit(1)


class VaspOptSetup:
    def __init__(self, config_file="input.toml", struct_file=None):
        self.work_dir = "vasp_opt"
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"找不到配置文件: {config_file}")
        if "tomllib" in globals():
            with open(config_file, "rb") as f:
                self.config = tomllib.load(f)
        else:
            self.config = toml.load(config_file)
        self.struct_file = struct_file
        self.vasp_inputpara = self._load_vasp_inputpara()

    def _load_vasp_inputpara(self):
        vasp_inputpara = self.config.get("vasp_inputpara")
        if not vasp_inputpara:
            raise ValueError("配置缺少 [vasp_inputpara] 块")
        return vasp_inputpara

    def _write_incar(self, filepath, lines):
        with open(filepath, "w") as incar:
            for line in lines:
                incar.write(line if line.endswith("\n") else line + "\n")

    def opt_incar1(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_1")
        lines = [
            "ISTART   = 0",
            "ICHARG   = 2",
            f"ISYM     = {self.vasp_inputpara['isym']}",
            "ENCUT    = 300",
            "PREC     = Normal",
            f"SYMPREC  = {self.vasp_inputpara['symprec']}",
            "NCORE    = 4",
            "KSPACING = 0.8",
            "ISMEAR   = 0",
            "SIGMA    = 0.2",
            "NELM     = 200",
            "NELMIN   = 6",
            "EDIFF    = 1e-3",
            "EDIFFG   = -0.2",
            "NSW      = 100",
            "IBRION   = 2",
            "ISIF     = 2",
            "POTIM    = 0.3",
            "LWAVE  = .FALSE.",
            "LCHARG = .FALSE.",
            f"PSTRESS  = {float(self.vasp_inputpara['press']) * 10}",
        ]
        self._write_incar(incar_filepath, lines)
        return incar_filepath

    def opt_incar2(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_2")
        lines = [
            "ISTART   = 0",
            "ICHARG   = 2",
            f"ISYM     = {self.vasp_inputpara['isym']}",
            "ENCUT    = 400",
            "PREC     = Normal",
            f"SYMPREC  = {self.vasp_inputpara['symprec']}",
            "NCORE    = 4",
            "KSPACING = 0.6",
            "ISMEAR   = 0",
            "SIGMA    = 0.1",
            "NELM     = 200",
            "NELMIN   = 6",
            "EDIFF    = 1e-4",
            "EDIFFG   = -0.1",
            "NSW      = 200",
            "IBRION   = 2",
            "ISIF     = 4",
            "POTIM    = 0.1",
            "LWAVE  = .FALSE.",
            "LCHARG = .FALSE.",
            f"PSTRESS  = {float(self.vasp_inputpara['press']) * 10}",
        ]
        self._write_incar(incar_filepath, lines)
        return incar_filepath

    def opt_incar3(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_3")
        lines = [
            "ISTART   = 0",
            "ICHARG   = 2",
            f"ISYM     = {self.vasp_inputpara['isym']}",
            "ENCUT    = 500",
            "PREC     = Accurate",
            f"SYMPREC  = {self.vasp_inputpara['symprec']}",
            "NCORE    = 4",
            "KSPACING = 0.4",
            "ISMEAR   = 0",
            "SIGMA    = 0.05",
            "NELM     = 90",
            "NELMIN   = 6",
            "EDIFF    = 1e-5",
            "EDIFFG   = -0.05",
            "NSW      = 300",
            "IBRION   = 2",
            "ISIF     = 3",
            "POTIM    = 0.05",
            "LWAVE    = .FALSE.",
            "LCHARG   = .FALSE.",
            f"PSTRESS  = {float(self.vasp_inputpara['press']) * 10}",
        ]
        self._write_incar(incar_filepath, lines)
        return incar_filepath

    def opt_incar4(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_4")
        lines = [
            "ISTART   = 0",
            "ICHARG   = 2",
            f"ISYM     = {self.vasp_inputpara['isym']}",
            f"ENCUT    = {self.vasp_inputpara['encut']}",
            "PREC     = Accurate",
            f"SYMPREC  = {self.vasp_inputpara['symprec']}",
        ]
        if self.vasp_inputpara.get("npar") is not None:
            lines.append(f"NPAR   = {self.vasp_inputpara['npar']}")
        if self.vasp_inputpara.get("ncore") is not None:
            lines.append(f"NCORE  = {self.vasp_inputpara['ncore']}")
        if self.vasp_inputpara.get("kpar") is not None:
            lines.append(f"KPAR   = {self.vasp_inputpara['kpar']}")
        lines.extend([
            f"KSPACING = {self.vasp_inputpara['kspacing']}",
            f"ISMEAR   = {self.vasp_inputpara['ismear']}",
            f"SIGMA    = {self.vasp_inputpara['sigma']}",
            f"NELM     = {self.vasp_inputpara['nelm']}",
            "NELMIN   = 6",
            f"EDIFF    = {self.vasp_inputpara['ediff']}",
            f"EDIFFG   = {self.vasp_inputpara['ediffg']}",
            "NSW      = 500",
            f"IBRION   = {self.vasp_inputpara['ibrion']}",
            f"ISIF     = {self.vasp_inputpara['isif']}",
            f"POTIM    = {self.vasp_inputpara['potim']}",
            "LWAVE  = .FALSE.",
            "LCHARG = .FALSE.",
            f"PSTRESS  = {float(self.vasp_inputpara['press']) * 10}",
        ])
        self._write_incar(incar_filepath, lines)
        return incar_filepath

    def _get_vasp_execmd(self):
        vasp_config = self.config.get("vasp", {})
        return vasp_config.get("executable_path", "vasp_std")

    def _get_slurm_header(self):
        slurm_config = self.config.get("slurm", {})
        return slurm_config.get("header", "")

    def _write_job_header(self, submit_file):
        header = self._get_slurm_header()
        if header:
            submit_file.write(header.strip() + "\n\n")
        else:
            submit_file.write("#!/bin/bash\n")

    def fouropt(self, submit_dirpath):
        jobname = "fouropt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        vasp_execmd = self._get_vasp_execmd()
        vaspstd_path = "vasp_std"
        with open(submit_script_filepath, "w") as submit:
            self._write_job_header(submit)
            submit.write("for i in {1..4}; do                  \n")
            submit.write("cp INCAR_$i INCAR                    \n")
            submit.write("killall -9 vasp_std                  \n")
            submit.write("sleep 3                              \n")
            submit.write(f"{vasp_execmd} {vaspstd_path} > vasp.log_$i 2>&1             \n")
            submit.write("cp CONTCAR CONTCAR_$i                \n")
            submit.write("cp OUTCAR  OUTCAR_$i                 \n")
            submit.write("cp CONTCAR POSCAR                    \n")
            submit.write("done                                 \n")
        os.chmod(submit_script_filepath, 0o755)
        return jobname

    def oneopt(self, submit_dirpath):
        jobname = "oneopt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        vasp_execmd = self._get_vasp_execmd()
        vaspstd_path = "vasp_std"
        with open(submit_script_filepath, "w") as submit:
            self._write_job_header(submit)
            submit.write(f"{vasp_execmd} {vaspstd_path} > vasp.log 2>&1               \n")
            submit.write("cp CONTCAR POSCAR\n")
        os.chmod(submit_script_filepath, 0o755)
        return jobname

    def fopt(self, submit_dirpath):
        jobname = "fopt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        vasp_execmd = self._get_vasp_execmd()
        vaspstd_path = "vasp_std"
        with open(submit_script_filepath, "w") as submit:
            self._write_job_header(submit)
            submit.write("num=0                                         \n")
            submit.write("while true;do                                 \n")
            submit.write("        let num+=1                            \n")
            submit.write("        echo \"run fine vasp opt-$num\"         \n")
            submit.write("        killall -9 vasp_std                                \n")
            submit.write("        sleep 3                                            \n")
            submit.write(f"        timeout 14400s {vasp_execmd} {vaspstd_path} > vasp.log 2>&1               \n")
            submit.write("        cp -f CONTCAR CONTCAR-fine &&  cp -f CONTCAR POSCAR\n")
            submit.write("        rows=`sed -n '/F\=/p' OSZICAR | wc -l`             \n")
            submit.write("        echo \"rows-$rows\"                                  \n")
            submit.write("        echo $num >> count_opt_times                       \n")
            submit.write("        if [ \"$rows\" -eq \"1\" ];then                        \n")
            submit.write("                break                                      \n")
            submit.write("        fi                                                 \n")
            submit.write("        if [ \"$num\" -gt \"50\" ]; then                       \n")
            submit.write("                break                                      \n")
            submit.write("        fi                                                 \n")
            submit.write("done                                                       \n")
        os.chmod(submit_script_filepath, 0o755)
        return jobname

    def _prepare_poscar(self):
        os.makedirs(self.work_dir, exist_ok=True)
        target_poscar = os.path.join(self.work_dir, "POSCAR")
        if self.struct_file:
            if not os.path.exists(self.struct_file):
                raise FileNotFoundError(f"找不到结构文件: {self.struct_file}")
            shutil.copy2(self.struct_file, target_poscar)
            print(f"已将 {self.struct_file} 拷贝为 {target_poscar}")
        elif os.path.exists("POSCAR") and os.path.abspath("POSCAR") != os.path.abspath(target_poscar):
            shutil.copy2("POSCAR", target_poscar)
            print(f"已将当前目录下的 POSCAR 拷贝为 {target_poscar}")

    def _copy_aux_files(self):
        for name in ["POTCAR", "KPOINTS"]:
            src = name
            dst = os.path.join(self.work_dir, name)
            if os.path.exists(src) and os.path.abspath(src) != os.path.abspath(dst):
                shutil.copy2(src, dst)

    def run(self, mode):
        self._prepare_poscar()
        self._copy_aux_files()
        if mode == "rv4":
            self.opt_incar1(self.work_dir)
            self.opt_incar2(self.work_dir)
            self.opt_incar3(self.work_dir)
            self.opt_incar4(self.work_dir)
            script = self.fouropt(self.work_dir)
        elif mode == "rv1":
            self.opt_incar4(self.work_dir)
            script = self.oneopt(self.work_dir)
        elif mode == "rvf":
            self.opt_incar4(self.work_dir)
            script = self.fopt(self.work_dir)
        else:
            raise ValueError("mode must be one of rv4/rv1/rvf")
        print(f"已生成优化输入于 {self.work_dir}，提交脚本: {script}")


def main():
    parser = argparse.ArgumentParser(description="VASP structure optimization (rv4/rv1/rvf)")
    parser.add_argument("-i", "--input", help="输入结构文件 (将被拷贝为 POSCAR)")
    parser.add_argument("-c", "--config", default="input.toml", help="配置文件路径 (默认为 input.toml)")
    parser.add_argument("--mode", choices=["rv4", "rv1", "rvf"], default="rv4", help="优化模式")
    args = parser.parse_args()

    setup = VaspOptSetup(config_file=args.config, struct_file=args.input)
    setup.run(mode=args.mode)


if __name__ == "__main__":
    main()
