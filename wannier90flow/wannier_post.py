#!/usr/bin/env python3
"""
postw90.x 相关流程：生成运行脚本并可选本地执行。
"""

import argparse
import os
import subprocess
from typing import Any, Dict, List, Tuple

from .wannier_init import DEFAULT_SLURM_HEADER
try:
    import tomllib  # Python 3.11+
except ImportError:
    tomllib = None
    try:
        import toml  # type: ignore
    except ImportError:
        toml = None


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


def load_post_config(path: str) -> Dict[str, Any]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到配置文件: {path}")

    with open(path, "r", encoding="utf-8") as f:
        raw_lines = f.read().splitlines()

    dos_projects: List[str] = []
    dos_projects_literal = None
    out_lines: List[str] = []
    current_section = ""
    dos_section_start = None
    dos_section_end = None

    for idx, line in enumerate(raw_lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            if current_section == "[dos]":
                dos_section_end = len(out_lines)
            current_section = stripped
            if current_section == "[dos]":
                dos_section_start = len(out_lines)
            out_lines.append(line)
            continue

        if current_section == "[dos]":
            if stripped.startswith("dos_project"):
                parts = stripped.split("=", 1)
                if len(parts) == 2:
                    value = parts[1].split("#", 1)[0].strip()
                    if value:
                        if value.startswith("["):
                            dos_projects_literal = value
                        elif value[0] in ("'", '"') and value[-1] == value[0]:
                            dos_projects.append(value)
                        else:
                            dos_projects.append(f'"{value}"')
                continue
        out_lines.append(line)

    if current_section == "[dos]":
        dos_section_end = len(out_lines)

    if dos_projects_literal:
        insert_at = dos_section_end if dos_section_end is not None else len(out_lines)
        out_lines.insert(insert_at, f"dos_project = {dos_projects_literal}")
    elif dos_projects:
        insert_at = dos_section_end if dos_section_end is not None else len(out_lines)
        merged = f"dos_project = [{', '.join(dos_projects)}]"
        out_lines.insert(insert_at, merged)

    merged_text = "\n".join(out_lines) + "\n"
    if tomllib:
        return tomllib.loads(merged_text)
    if toml:
        return toml.loads(merged_text)
    raise ImportError("缺少 tomllib/toml 模块，无法读取配置。")


def format_win_value(key: str, value: Any) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, (int, float)):
        return str(value)
    if key == "dos_project" and isinstance(value, (list, tuple)):
        return ", ".join(str(x) for x in value)
    if key == "dos_kmesh" and isinstance(value, (list, tuple)):
        return " ".join(str(int(x)) for x in value)
    return str(value)


def update_win_content(content: str, updates: Dict[str, Any]) -> str:
    keys = list(updates.keys())
    lines = []
    for line in content.splitlines():
        stripped = line.lstrip()
        if any(stripped.startswith(f"{key} ") or stripped.startswith(f"{key}=") for key in keys):
            continue
        lines.append(line)
    for key in keys:
        value = format_win_value(key, updates[key])
        lines.append(f"{key} = {value}")
    return "\n".join(lines) + "\n"


def sanitize_project_tag(value: str) -> str:
    tag = value.strip().strip('"').strip("'")
    for ch in [":", ",", " "]:
        tag = tag.replace(ch, "_")
    while "__" in tag:
        tag = tag.replace("__", "_")
    return tag.strip("_")


def resolve_dos_output(work_dir: str, prefix: str) -> str:
    cand_underscore = os.path.join(work_dir, f"{prefix}_dos.dat")
    cand_hyphen = os.path.join(work_dir, f"{prefix}-dos.dat")
    if os.path.exists(cand_underscore):
        return cand_underscore
    if os.path.exists(cand_hyphen):
        return cand_hyphen
    return ""


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

    cfg = load_post_config(args.config)
    win_path, prefix = resolve_win_path(args.input)
    work_dir = os.path.dirname(os.path.abspath(win_path)) or "."

    create_run_script(work_dir, prefix, cfg)
    if args.run:
        with open(win_path, "r", encoding="utf-8") as f:
            original_content = f.read()
        dos_cfg = cfg.get("dos", {})
        dos_projects = dos_cfg.get("dos_project")
        projects: List[Any] = []
        if isinstance(dos_projects, (list, tuple)):
            projects = list(dos_projects)
        elif dos_projects is not None:
            projects = [dos_projects]
        else:
            projects = []

        common_updates = {}
        for key, value in dos_cfg.items():
            if key == "dos_project":
                continue
            common_updates[f"dos_{key}" if not key.startswith("dos") else key] = value
        if "dos" not in common_updates:
            common_updates["dos"] = True

        try:
            if projects:
                for proj in projects:
                    updates = dict(common_updates)
                    updates["dos_project"] = proj
                    content = update_win_content(original_content, updates)
                    with open(win_path, "w", encoding="utf-8") as f:
                        f.write(content)
                    run_postw90(work_dir, prefix, cfg)
                    output_path = resolve_dos_output(work_dir, prefix)
                    if output_path:
                        tag = sanitize_project_tag(str(proj))
                        target = os.path.join(work_dir, f"{prefix}_dos_{tag}.dat")
                        os.replace(output_path, target)
                        print(f"已保存投影 DOS: {target}")
            else:
                content = update_win_content(original_content, common_updates)
                with open(win_path, "w", encoding="utf-8") as f:
                    f.write(content)
                run_postw90(work_dir, prefix, cfg)
        finally:
            with open(win_path, "w", encoding="utf-8") as f:
                f.write(original_content)


if __name__ == "__main__":
    main()
