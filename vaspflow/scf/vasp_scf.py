#!/usr/bin/env python3
import os
import shutil
import argparse
import math
import subprocess

try:
    import tomllib  # Python 3.11+
except ImportError:
    try:
        import toml
    except ImportError:
        print("Error: 'toml' or 'tomllib' (Python 3.11+) is required.")
        print("Please install it using: pip install toml")
        exit(1)

class VaspSetup:
    def __init__(self, config_file="input.toml", struct_file=None):
        self.work_dir = "vasp_scf"
        # 1. 加载配置文件
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"找不到配置文件: {config_file}")
        
        if 'tomllib' in globals():
            with open(config_file, "rb") as f:
                self.config = tomllib.load(f)
        else:
            self.config = toml.load(config_file)
        
        self.struct_file = struct_file
        
        # 2. 默认 INCAR 模板
        self.incar_template = {
            "SYSTEM" : "VASP_Job",
            "ISTART" : 0,
            "ICHARG" : 2,
            "PREC"   : "Accurate",
            "IBRION" : -1,
            "LWAVE"  : ".TRUE.",
            "LCHARG" : ".TRUE.",
            "NCORE"  : 4
        }

    def get_elements_from_poscar(self):
        """从 POSCAR 中读取元素名称和顺序"""
        poscar_path = os.path.join(self.work_dir, "POSCAR")
        if not os.path.exists(poscar_path):
            raise FileNotFoundError(f"找不到 {poscar_path} 文件")
        
        with open(poscar_path, "r") as f:
            lines = f.readlines()
            # 这里的逻辑假设 POSCAR 符合 VASP 5 标准（第6行是元素名称）
            # 第1行: 注释
            # 第2行: 缩放因子
            # 第3-5行: 晶格矢量
            # 第6行: 元素名称
            elements = lines[5].split()
            # 如果第6行全是数字，说明是旧版 POSCAR，元素名可能在第1行或者缺失
            if elements[0].isdigit():
                print("警告: 检测到可能不包含元素名称的 POSCAR (VASP 4 格式)。")
                print("请确保 POSCAR 第6行包含元素名称。")
            return elements

    def generate_incar(self):
        """生成 INCAR 文件"""
        # 合并模板和配置文件中的参数
        final_params = self.incar_template.copy()
            
        # 1. 合并 [incar_params] 中的参数 (最高优先级)
        if "incar_params" in self.config:
            final_params.update(self.config["incar_params"])
        
        print(f"正在生成 INCAR...")
        incar_path = os.path.join(self.work_dir, "INCAR")
        with open(incar_path, "w") as f:
            for key, value in final_params.items():
                # 处理字符串布尔值或特殊格式
                f.write(f"{key:12} = {value}\n")

    def generate_potcar(self, elements):
        """根据 POSCAR 顺序从指定目录合并 POTCAR"""
        print(f"检测到元素顺序: {' '.join(elements)}")
        potcar_dir = self.config.get("potcar_dir", "")
        if not potcar_dir:
            raise ValueError("未在配置文件中找到 'potcar_dir' (赝势目录路径)")
        
        potcar_dir = os.path.expanduser(potcar_dir)
        potcar_content = []
        for el in elements:
            path = os.path.join(potcar_dir, el)
            if not os.path.exists(path):
                raise FileNotFoundError(f"在 {potcar_dir} 中找不到元素 {el} 的 POTCAR 文件")
            
            with open(path, "r") as f:
                potcar_content.append(f.read())
        
        potcar_path = os.path.join(self.work_dir, "POTCAR")
        with open(potcar_path, "w") as f:
            f.writelines(potcar_content)
        print(f"POTCAR 已从 {potcar_dir} 合并完成。")

    def generate_kpoints(self):
        """使用 vaspkit 生成 KPOINTS 文件"""
        kmesh = self.config.get("kmesh")
        if not kmesh:
            print("未在配置文件中找到 'kmesh'，跳过生成 KPOINTS。")
            return

        print(f"检测到 kmesh = {kmesh}，正在使用 vaspkit 生成 KPOINTS...")
        try:
            process = subprocess.Popen(
                ["vaspkit"],
                cwd=self.work_dir,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            # vaspkit 102 -> 2 -> kmesh
            stdout, stderr = process.communicate(input=f"102\n2\n{kmesh}\n")
            
            kpoints_file = os.path.join(self.work_dir, "KPOINTS")
            if os.path.exists(kpoints_file):
                print(f"KPOINTS 已通过 vaspkit 生成完成。")
            else:
                print(f"错误: vaspkit 未生成 KPOINTS。")
                print(f"vaspkit 输出: {stdout}")
        except Exception as e:
            print(f"调用 vaspkit 出错: {e}")

    def create_run_script(self):
        """创建一个提交脚本 (支持 Slurm)"""
        vasp_config = self.config.get("vasp", {})
        vasp_path = vasp_config.get("executable_path", "vasp_std")
        
        slurm_config = self.config.get("slurm", {})
        slurm_header = slurm_config.get("header", "")

        run_script_path = os.path.join(self.work_dir, "run_vasp.sh")
        with open(run_script_path, "w") as f:
            if slurm_header:
                f.write(slurm_header.strip() + "\n\n")
            else:
                # 默认简单的 bash 头部
                f.write("#!/bin/bash\n")
            
            f.write(f"echo 'Job started at' `date` \n")
            f.write(f"{vasp_path}\n")
            f.write(f"echo 'Job finished at' `date` \n")
        
        os.chmod(run_script_path, 0o755)
        print(f"生成运行脚本: {run_script_path}")

    def execute_vasp(self):
        """进入工作目录并执行 VASP 命令"""
        vasp_config = self.config.get("vasp", {})
        vasp_path = vasp_config.get("executable_path", "vasp_std")
        
        original_dir = os.getcwd()
        try:
            print(f"\n正在进入 {self.work_dir} 目录执行: {vasp_path}")
            os.chdir(self.work_dir)
            # 使用 shell=True 因为 executable_path 可能包含 mpirun 等复杂命令
            subprocess.run(vasp_path, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"执行 VASP 出错: {e}")
        finally:
            os.chdir(original_dir)
            print(f"已返回目录: {original_dir}")

    def run(self, run_calc=False):
        try:
            # 确保工作目录存在
            os.makedirs(self.work_dir, exist_ok=True)

            # 0. 如果指定了输入结构文件，先拷贝为 POSCAR
            target_poscar = os.path.join(self.work_dir, "POSCAR")
            if self.struct_file:
                if not os.path.exists(self.struct_file):
                    raise FileNotFoundError(f"找不到结构文件: {self.struct_file}")
                shutil.copy2(self.struct_file, target_poscar)
                print(f"已将 {self.struct_file} 拷贝为 {target_poscar}")
            elif os.path.exists("POSCAR") and os.path.abspath("POSCAR") != os.path.abspath(target_poscar):
                shutil.copy2("POSCAR", target_poscar)
                print(f"已将当前目录下的 POSCAR 拷贝为 {target_poscar}")

            elements = self.get_elements_from_poscar()
            self.generate_incar()
            self.generate_potcar(elements)
            self.generate_kpoints()
            self.create_run_script()
            print(f"\n所有 VASP 输入文件已在 {self.work_dir} 目录中准备就绪！")

            if run_calc:
                self.execute_vasp()
                
        except Exception as e:
            print(f"错误: {e}")

def main():
    parser = argparse.ArgumentParser(description="VASP SCF Setup Script")
    parser.add_argument("-i", "--input", help="输入结构文件 (将被拷贝为 POSCAR)")
    parser.add_argument("-c", "--config", default="input.toml", help="配置文件路径 (默认为 input.toml)")
    parser.add_argument("--run", action="store_true", help="生成文件后直接执行计算")
    
    args = parser.parse_args()
    
    setup = VaspSetup(config_file=args.config, struct_file=args.input)
    setup.run(run_calc=args.run)

if __name__ == "__main__":
    main()
