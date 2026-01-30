#!/usr/bin/env python3
"""
postw90.x 相关流程：生成运行脚本并可选本地执行。
"""

import argparse
import os
import subprocess
from typing import Any, Dict, Tuple

from .wannier_init import DEFAULT_SLURM_HEADER, load_config


def resolve_win_path(input_path: str) -> Tuple[str, str]:
    win_path = input_path
    if os.path.isdir(win_path):
        raise ValueError("输入应为 *.win 文件路径，而不是目录。")
    if not win_path.endswith(".win"):
        win_path = f"{win_path}.win"
    if not os.path.exists(win_path):
        raise FileNotFoundError(f"未找到 {win_path}，请检查 -i 参数。")
    prefix = os.path.splitext(os.path.basename(win_path))[0]
    return win_path, prefix


def create_run_script(work_dir: str, prefix: str, cfg: Dict[str, Any]) -> None:
    post_cfg = cfg.get("postw90", {})
    post_exec_path = post_cfg.get("executable_path", "postw90.x")
    script_path = os.path.join(work_dir, "run_postw90.sh")
    with open(script_path, "w", encoding="utf-8") as f:
        f.write(DEFAULT_SLURM_HEADER.strip() + "\n\n")
        f.write("set -e\n")
        f.write('cd "$(dirname "$0")"\n\n')
        f.write('echo "Run postw90 at $(date)"\n')
        f.write(f'{post_exec_path} "{prefix}" > POSTW90.log 2>&1\n')
    os.chmod(script_path, 0o755)
    print(f"已生成提交脚本: {script_path}")


def run_postw90(work_dir: str, prefix: str, cfg: Dict[str, Any]) -> None:
    post_cfg = cfg.get("postw90", {})
    post_exec_path = post_cfg.get("executable_path", "postw90.x")
    cmd_post = f'{post_exec_path} "{prefix}" > POSTW90.log 2>&1'
    print(f"执行命令: {cmd_post}")
    subprocess.run(cmd_post, shell=True, check=True, cwd=work_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description="生成/执行 postw90.x 流程")
    parser.add_argument(
        "-c",
        "--config",
        default="inputwannier.toml",
        help="配置文件路径（TOML/JSON）",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="输入 .win 文件路径（或不带扩展名的前缀）",
    )
    parser.add_argument("--run", action="store_true", help="生成脚本后直接执行 postw90.x")
    args = parser.parse_args()

    cfg = load_config(args.config)
    win_path, prefix = resolve_win_path(args.input)
    work_dir = os.path.dirname(os.path.abspath(win_path)) or "."

    create_run_script(work_dir, prefix, cfg)
    if args.run:
        run_postw90(work_dir, prefix, cfg)


if __name__ == "__main__":
    main()
