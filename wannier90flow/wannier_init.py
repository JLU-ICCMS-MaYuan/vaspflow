#!/usr/bin/env python3
"""
准备 Wannier90 初始输入（*.win）。

功能：
- 读取 POSCAR 获取晶格与原子坐标（保留 Direct/Cartesian 体系）。
- 读取配置 (TOML/JSON) 设置投影、能窗、k 点网格等；kpoint_path 优先从 vaspkit 303 生成的 KPATH.in 读取。
- 生成包含 unit_cell_cart、atoms_*、projections、mp_grid、kpoints、可选 kpoint_path 的 .win 文件。

示例：
python wannier90flow/wannier_init.py \
  --input POSCAR \
  --config wannier90.toml
"""

import argparse
import json
import os
import subprocess
from typing import Any, Dict, List, Tuple

import numpy as np

DEFAULT_SLURM_HEADER = """#!/bin/bash
#SBATCH --job-name=wannier90
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
"""

try:
    import tomllib  # Python 3.11+
except ImportError:
    tomllib = None
    try:
        import toml  # type: ignore
    except ImportError:
        toml = None


def load_config(path: str) -> Dict[str, Any]:
    """读取 TOML 或 JSON 配置。"""
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到配置文件: {path}")

    if path.endswith(".json"):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    if tomllib:
        with open(path, "rb") as f:
            return tomllib.load(f)
    if toml:
        return toml.load(path)
    raise ImportError("缺少 tomllib/toml 模块，且配置不是 JSON。")


def parse_poscar(path: str) -> Dict[str, Any]:
    """解析 POSCAR，保留 Direct/Cartesian 坐标类型。"""
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到 POSCAR: {path}")

    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.rstrip() for ln in f if ln.strip() != ""]

    if len(lines) < 8:
        raise ValueError("POSCAR 行数不足，无法解析。")

    comment = lines[0]
    scale = float(lines[1])

    lattice = []
    for i in range(2, 5):
        lattice.append([float(x) * scale for x in lines[i].split()])

    elements = lines[5].split()
    counts = [int(x) for x in lines[6].split()]

    coord_line = lines[7].strip().lower()
    start_line = 8
    if coord_line.startswith("s"):  # Selective dynamics
        coord_line = lines[8].strip().lower()
        start_line = 9

    total_atoms = sum(counts)
    positions: List[List[float]] = []
    for idx in range(total_atoms):
        positions.append([float(x) for x in lines[start_line + idx].split()[:3]])

    coords = np.array(positions, dtype=float)
    coord_type = "direct" if coord_line.startswith("d") else "cartesian"
    if coord_type == "cartesian":
        coords *= scale

    return {
        "comment": comment,
        "lattice": np.array(lattice),
        "elements": elements,
        "counts": counts,
        "coord_type": coord_type,
        "positions": coords,
    }


def generate_kpoints_from_grid(grid: List[int]) -> List[Tuple[float, float, float]]:
    """生成 Gamma-centered 均匀网格的 k 点列表。"""
    if len(grid) != 3:
        raise ValueError("mp_grid 必须包含 3 个整数。")
    n1, n2, n3 = [int(x) for x in grid]
    if min(n1, n2, n3) <= 0:
        raise ValueError("mp_grid 中的每个值必须 > 0。")

    kpts: List[Tuple[float, float, float]] = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                kpts.append((i / n1, j / n2, k / n3))
    return kpts


def parse_kpath(path: str) -> List[Dict[str, Any]]:
    """
    解析 vaspkit 303 生成的 KPATH.in，提取能带路径段。
    期望行格式：Label1 k1x k1y k1z Label2 k2x k2y k2z
    """
    segments: List[Dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 8:
                # 兼容含有前缀数字的行，如 "4\n" 或 header
                continue
            from_label, k1x, k1y, k1z, to_label, k2x, k2y, k2z = parts
            segments.append(
                {
                    "from": from_label,
                    "from_k": [float(k1x), float(k1y), float(k1z)],
                    "to": to_label,
                    "to_k": [float(k2x), float(k2y), float(k2z)],
                }
            )
    return segments


def format_vector(vec: List[float]) -> str:
    return f"{vec[0]:16.10f} {vec[1]:16.10f} {vec[2]:16.10f}"


def format_kpoint(pt: Tuple[float, float, float]) -> str:
    return f"{pt[0]:12.8f} {pt[1]:12.8f} {pt[2]:12.8f}"


def create_run_script(work_dir: str, prefix: str, cfg: Dict[str, Any]) -> None:
    """生成提交脚本，包含 -pp 预处理与主计算占位。"""
    w90_cfg = cfg.get("wannier90", {})
    exec_path = w90_cfg.get("executable_path", "wannier90.x")
    pw2_cfg = cfg.get("pw2wannier90", {})
    pw2_exec_path = pw2_cfg.get("executable_path", "pw2wannier90.x")
    pw2_input_path, pw2_input_name = resolve_pw2_input_path(work_dir, prefix, cfg)
    ensure_pw2wan_input(pw2_input_path, prefix, cfg)
    slurm_header = DEFAULT_SLURM_HEADER

    script_path = os.path.join(work_dir, "run_wannier90.sh")
    with open(script_path, "w", encoding="utf-8") as f:
        f.write(slurm_header.strip() + "\n\n")
        f.write("set -e\n")
        f.write('cd "$(dirname "$0")"\n\n')
        f.write('echo "Prep NNKP at $(date)"\n')
        f.write(f'{exec_path} -pp "{prefix}"\n')
        f.write('\n')
        f.write('echo "Run pw2wannier90 at $(date)"\n')
        f.write(f'if [ ! -f "{pw2_input_name}" ]; then\n')
        f.write(f'  echo "缺少 pw2wannier90 输入文件: {pw2_input_name}"\n')
        f.write('  exit 1\n')
        f.write('fi\n')
        f.write(f'{pw2_exec_path} < "{pw2_input_name}" > PW2WAN.log 2>&1\n')
        f.write('\n')
        f.write('echo "Run Wannier90 (requires *.amn/*.mmn ready)"\n')
        f.write(f'{exec_path} "{prefix}" > WANNIER90.log 2>&1\n')
    os.chmod(script_path, 0o755)
    print(f"已生成提交脚本: {script_path}")


def format_pw2_value(value: Any) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, (int, float)):
        return str(value)
    return f"'{value}'"


def resolve_pw2_input_path(
    work_dir: str, prefix: str, cfg: Dict[str, Any]
) -> Tuple[str, str]:
    pw2_cfg = cfg.get("pw2wannier90", {})
    input_file = pw2_cfg.get("input_file")
    candidates: List[str] = []
    if input_file:
        candidates.append(input_file)
    else:
        candidates.extend([f"{prefix}.pw2wan", f"{prefix}.pw2wan.in"])

    for cand in candidates:
        path = cand if os.path.isabs(cand) else os.path.join(work_dir, cand)
        if os.path.exists(path):
            return path, os.path.basename(path)

    default_name = input_file or f"{prefix}.pw2wan"
    path = default_name if os.path.isabs(default_name) else os.path.join(work_dir, default_name)
    return path, os.path.basename(path)


def ensure_pw2wan_input(path: str, prefix: str, cfg: Dict[str, Any]) -> None:
    """确保 pw2wannier90 输入文件存在，若不存在则生成默认文件。"""
    if os.path.exists(path):
        print(f"pw2wannier90 输入文件已存在: {path}")
        return
    pw2_cfg = cfg.get("pw2wannier90", {})
    outdir = pw2_cfg.get("outdir", "./")
    prefix_val = pw2_cfg.get("prefix", prefix)
    seedname = pw2_cfg.get("seedname", prefix)

    lines = [
        "&inputpp",
        f"  outdir = {format_pw2_value(outdir)}",
        f"  prefix = {format_pw2_value(prefix_val)}",
        f"  seedname = {format_pw2_value(seedname)}",
    ]
    if "write_amn" in pw2_cfg:
        lines.append(f"  write_amn = {format_pw2_value(pw2_cfg['write_amn'])}")
    else:
        lines.append("  write_amn = .true.")
    if "write_mmn" in pw2_cfg:
        lines.append(f"  write_mmn = {format_pw2_value(pw2_cfg['write_mmn'])}")
    else:
        lines.append("  write_mmn = .true.")
    if "write_unk" in pw2_cfg:
        lines.append(f"  write_unk = {format_pw2_value(pw2_cfg['write_unk'])}")
    if "spin_component" in pw2_cfg:
        lines.append(f"  spin_component = {format_pw2_value(pw2_cfg['spin_component'])}")
    if "wan_mode" in pw2_cfg:
        lines.append(f"  wan_mode = {format_pw2_value(pw2_cfg['wan_mode'])}")
    lines.append("/\n")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"已生成 pw2wannier90 输入文件: {path}")


def run_wannier90_pipeline(work_dir: str, prefix: str, cfg: Dict[str, Any]) -> None:
    w90_cfg = cfg.get("wannier90", {})
    exec_path = w90_cfg.get("executable_path", "wannier90.x")
    pw2_cfg = cfg.get("pw2wannier90", {})
    pw2_exec_path = pw2_cfg.get("executable_path", "pw2wannier90.x")
    pw2_input_path, pw2_input_name = resolve_pw2_input_path(work_dir, prefix, cfg)
    ensure_pw2wan_input(pw2_input_path, prefix, cfg)

    print("开始执行 Wannier90 流程（-pp -> pw2wannier90 -> wannier90）...")
    cmd_pp = f'{exec_path} -pp "{prefix}"'
    print(f"执行命令: {cmd_pp}")
    subprocess.run(cmd_pp, shell=True, check=True, cwd=work_dir)
    cmd_pw2 = f'{pw2_exec_path} < "{pw2_input_name}" > PW2WAN.log 2>&1'
    print(f"执行命令: {cmd_pw2}")
    subprocess.run(cmd_pw2, shell=True, check=True, cwd=work_dir)
    cmd_w90 = f'{exec_path} "{prefix}.win" > WANNIER90.log 2>&1'
    print(f"执行命令: {cmd_w90}")
    subprocess.run(cmd_w90, shell=True, check=True, cwd=work_dir)


def write_win(
    output: str,
    struct: Dict[str, Any],
    cfg: Dict[str, Any],
    kpoints: List[Tuple[float, float, float]],
    band_path: List[Dict[str, Any]],
) -> None:
    """写出 .win 文件。"""
    win_cfg = cfg.get("win", {})
    system_name = win_cfg.get("system_name") or cfg.get("system_name") or struct["comment"] or "wannier_system"
    projections = cfg.get("projections", {}).get("list") or cfg.get("projections") or []
    kpt_cfg = cfg.get("k_points", {})

    # 必选参数检查
    if not projections:
        raise ValueError("配置缺少投影 (projections)。")

    num_wann = win_cfg.get("num_wann", cfg.get("num_wann"))
    if num_wann is None:
        raise ValueError("配置缺少 num_wann。")

    with open(output, "w", encoding="utf-8") as f:
        f.write(f"num_wann = {int(num_wann)}\n")
        num_bands = win_cfg.get("num_bands", cfg.get("num_bands"))
        if num_bands is not None:
            f.write(f"num_bands = {int(num_bands)}\n")
        exclude_bands = win_cfg.get("exclude_bands", cfg.get("exclude_bands"))
        if exclude_bands is not None:
            f.write(f"exclude_bands = {exclude_bands}\n")
        iprint = win_cfg.get("iprint", cfg.get("iprint"))
        if iprint is not None:
            f.write(f"iprint = {int(iprint)}\n")
        num_iter = win_cfg.get("num_iter", cfg.get("num_iter"))
        if num_iter is not None:
            f.write(f"num_iter = {int(num_iter)}\n")
        dis_num_iter = win_cfg.get("dis_num_iter", cfg.get("dis_num_iter"))
        if dis_num_iter is not None:
            f.write(f"dis_num_iter = {int(dis_num_iter)}\n")
        dis_froz_min = win_cfg.get("dis_froz_min", cfg.get("dis_froz_min"))
        if dis_froz_min is not None:
            f.write(f"dis_froz_min = {float(dis_froz_min):12.6f}\n")
        dis_froz_max = win_cfg.get("dis_froz_max", cfg.get("dis_froz_max"))
        if dis_froz_max is not None:
            f.write(f"dis_froz_max = {float(dis_froz_max):12.6f}\n")
        dis_win_min = win_cfg.get("dis_win_min", cfg.get("dis_win_min"))
        if dis_win_min is not None:
            f.write(f"dis_win_min = {float(dis_win_min):12.6f}\n")
        dis_win_max = win_cfg.get("dis_win_max", cfg.get("dis_win_max"))
        if dis_win_max is not None:
            f.write(f"dis_win_max = {float(dis_win_max):12.6f}\n")
        write_hr = win_cfg.get("write_hr", cfg.get("write_hr"))
        if write_hr is not None:
            f.write(f"write_hr = {str(write_hr).lower()}\n")
        write_bvec = win_cfg.get("write_bvec", cfg.get("write_bvec"))
        if write_bvec is not None:
            f.write(f"write_bvec = {str(write_bvec).lower()}\n")
        bands_plot = win_cfg.get("bands_plot", cfg.get("bands_plot"))
        bands_plot_format = win_cfg.get("bands_plot_format", cfg.get("bands_plot_format"))
        has_band_path = len(band_path) > 0
        if bands_plot is not None:
            if bands_plot and not has_band_path:
                print("提示：已请求 bands_plot，但未找到 kpoint_path，将自动关闭 bands_plot。")
                bands_plot = False
            f.write(f"bands_plot = {str(bands_plot).lower()}\n")
        if bands_plot_format is not None and bands_plot:
            f.write(f"bands_plot_format = {bands_plot_format}\n")
        wvfn_formatted = win_cfg.get("wvfn_formatted", cfg.get("wvfn_formatted"))
        if wvfn_formatted is not None:
            f.write(f"wvfn_formatted = {str(wvfn_formatted).lower()}\n")
        wannier_plot = win_cfg.get("wannier_plot", cfg.get("wannier_plot"))
        if wannier_plot is not None:
            f.write(f"wannier_plot = {str(wannier_plot).lower()}\n")
        wannier_plot_format = win_cfg.get("wannier_plot_format", cfg.get("wannier_plot_format"))
        if wannier_plot_format is not None and wannier_plot:
            f.write(f"wannier_plot_format = {wannier_plot_format}\n")
        fermi_energy = win_cfg.get("fermi_energy", cfg.get("fermi_energy"))
        if fermi_energy is not None:
            f.write(f"fermi_energy = {float(fermi_energy):12.6f}\n")
        fermi_surface_plot = win_cfg.get("fermi_surface_plot", cfg.get("fermi_surface_plot"))
        if fermi_surface_plot is not None:
            f.write(f"fermi_surface_plot = {str(fermi_surface_plot).lower()}\n")

        f.write(f"mp_grid = {' '.join(str(int(x)) for x in kpt_cfg.get('mp_grid', []))}\n")
        f.write("\n")

        f.write("begin projections\n")
        for proj in projections:
            f.write(f" {proj}\n")
        f.write("end projections\n\n")

        f.write("begin unit_cell_cart\n")
        f.write("Ang\n")
        for vec in struct["lattice"]:
            f.write(f"{format_vector(vec.tolist())}\n")
        f.write("end unit_cell_cart\n\n")

        coord_type = struct["coord_type"]
        if coord_type == "direct":
            f.write("begin atoms_frac\n")
        else:
            f.write("begin atoms_cart\n")
            f.write("Ang\n")

        atom_idx = 0
        for el, count in zip(struct["elements"], struct["counts"]):
            for _ in range(count):
                pos = struct["positions"][atom_idx]
                f.write(f"{el:3} {format_kpoint(tuple(pos))}\n")
                atom_idx += 1
        if coord_type == "direct":
            f.write("end atoms_frac\n\n")
        else:
            f.write("end atoms_cart\n\n")

        f.write("begin kpoints\n")
        for pt in kpoints:
            f.write(f"{format_kpoint(pt)}\n")
        f.write("end kpoints\n\n")

        if band_path:
            f.write("bands_plot = .true.\n")
            f.write("begin kpoint_path\n")
            for seg in band_path:
                from_label = seg.get("from", "G")
                to_label = seg.get("to", "G")
                k1 = seg.get("from_k", [0.0, 0.0, 0.0])
                k2 = seg.get("to_k", [0.0, 0.0, 0.0])
                f.write(
                    f"{from_label:2} {k1[0]:8.4f} {k1[1]:8.4f} {k1[2]:8.4f} "
                    f"{to_label:2} {k2[0]:8.4f} {k2[1]:8.4f} {k2[2]:8.4f}\n"
                )
            f.write("end kpoint_path\n\n")

        f.write(f"! system_name: {system_name}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="准备 Wannier90 .win 输入文件")
    parser.add_argument("-i", "--input", default="POSCAR", help="结构文件 (VASP POSCAR)")
    parser.add_argument("-c", "--config", required=True, help="配置文件 (TOML/JSON)")
    parser.add_argument("--run", action="store_true", help="生成文件后直接执行 Wannier90 流程")
    args = parser.parse_args()

    work_dir = "wannier90"
    os.makedirs(work_dir, exist_ok=True)

    struct = parse_poscar(args.input)
    cfg = load_config(args.config)
    kpt_cfg = cfg.get("k_points", {})
    win_cfg = cfg.get("win", {})

    explicit_kpts = kpt_cfg.get("kpoints")
    if explicit_kpts:
        kpoints = [tuple(map(float, kp)) for kp in explicit_kpts]
    elif "mp_grid" in kpt_cfg:
        kpoints = generate_kpoints_from_grid(kpt_cfg["mp_grid"])
    else:
        raise ValueError("k_points 配置缺少 mp_grid 或 kpoints。")

    bands_plot_cfg = win_cfg.get("bands_plot", cfg.get("bands_plot"))
    kpath_file_cfg = kpt_cfg.get("kpath_file", "KPATH.in")
    kpath_file = (
        kpath_file_cfg
        if os.path.isabs(kpath_file_cfg)
        else os.path.join(work_dir, kpath_file_cfg)
    )
    band_segments: List[Dict[str, Any]] = []
    if bands_plot_cfg and not os.path.exists(kpath_file):
        print("检测到 bands_plot=true，开始运行 vaspkit 303 生成 KPATH.in ...")
        try:
            subprocess.run(
                ["vaspkit"],
                input="3\n303\n",
                text=True,
                check=True,
                cwd=work_dir,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except Exception as exc:
            raise RuntimeError("vaspkit 运行失败，无法生成 KPATH.in") from exc
    if os.path.exists(kpath_file):
        band_segments = parse_kpath(kpath_file)
        print(f"从 {kpath_file} 读取到 {len(band_segments)} 段能带路径。")
    else:
        print(f"未找到 {kpath_file}，将不写入 kpoint_path。")

    win_cfg = cfg.get("win", {})
    basename = win_cfg.get("system_name") or cfg.get("system_name") or struct["comment"] or "wannier90"
    output = os.path.join(work_dir, f"{basename}.win")

    write_win(output, struct, cfg, kpoints, band_segments)
    create_run_script(work_dir, basename, cfg)
    print(f"已生成 {output}，包含 {len(kpoints)} 个 k 点。")
    if args.run:
        run_wannier90_pipeline(work_dir, basename, cfg)


if __name__ == "__main__":
    main()
