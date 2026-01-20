#!/usr/bin/env python3
"""
准备 Wannier90 初始输入（*.win）。

功能：
- 读取 POSCAR 获取晶格与原子坐标（保留 Direct/Cartesian 体系）。
- 读取配置 (TOML/JSON) 设置投影、能窗、k 点网格等；kpoint_path 优先从 vaspkit 303 生成的 KPATH.in 读取。
- 生成包含 unit_cell_cart、atoms_*、projections、mp_grid、kpoints、可选 kpoint_path 的 .win 文件。

示例：
python wannier90flow/prepare_wannier90.py \
  --poscar POSCAR \
  --config wannier90.toml
"""

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Tuple

import numpy as np

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


def write_win(
    output: str,
    struct: Dict[str, Any],
    cfg: Dict[str, Any],
    kpoints: List[Tuple[float, float, float]],
    band_path: List[Dict[str, Any]],
) -> None:
    """写出 .win 文件。"""
    system_name = cfg.get("system_name") or struct["comment"] or "wannier_system"
    projections = cfg.get("projections", {}).get("list") or cfg.get("projections") or []
    kpt_cfg = cfg.get("k_points", {})

    # 必选参数检查
    if not projections:
        raise ValueError("配置缺少投影 (projections)。")

    num_wann = cfg.get("num_wann")
    if num_wann is None:
        raise ValueError("配置缺少 num_wann。")

    with open(output, "w", encoding="utf-8") as f:
        f.write(f"num_wann = {int(num_wann)}\n")
        if "num_bands" in cfg:
            f.write(f"num_bands = {int(cfg['num_bands'])}\n")
        if "exclude_bands" in cfg:
            f.write(f"exclude_bands = {cfg['exclude_bands']}\n")
        if "iprint" in cfg:
            f.write(f"iprint = {int(cfg['iprint'])}\n")
        if "num_iter" in cfg:
            f.write(f"num_iter = {int(cfg['num_iter'])}\n")
        if "dis_num_iter" in cfg:
            f.write(f"dis_num_iter = {int(cfg['dis_num_iter'])}\n")
        if "dis_froz_min" in cfg:
            f.write(f"dis_froz_min = {float(cfg['dis_froz_min']):12.6f}\n")
        if "dis_froz_max" in cfg:
            f.write(f"dis_froz_max = {float(cfg['dis_froz_max']):12.6f}\n")
        if "dis_win_min" in cfg:
            f.write(f"dis_win_min = {float(cfg['dis_win_min']):12.6f}\n")
        if "dis_win_max" in cfg:
            f.write(f"dis_win_max = {float(cfg['dis_win_max']):12.6f}\n")
        if "write_hr" in cfg:
            f.write(f"write_hr = {str(cfg['write_hr']).lower()}\n")
        if "write_bvec" in cfg:
            f.write(f"write_bvec = {str(cfg['write_bvec']).lower()}\n")
        if "bands_plot" in cfg:
            f.write(f"bands_plot = {str(cfg['bands_plot']).lower()}\n")
        if "bands_plot_format" in cfg:
            f.write(f"bands_plot_format = {cfg['bands_plot_format']}\n")

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
    parser.add_argument("-p", "--poscar", default="POSCAR", help="结构文件 (VASP POSCAR)")
    parser.add_argument("-c", "--config", required=True, help="配置文件 (TOML/JSON)")
    args = parser.parse_args()

    struct = parse_poscar(args.poscar)
    cfg = load_config(args.config)
    kpt_cfg = cfg.get("k_points", {})

    explicit_kpts = kpt_cfg.get("kpoints")
    if explicit_kpts:
        kpoints = [tuple(map(float, kp)) for kp in explicit_kpts]
    elif "mp_grid" in kpt_cfg:
        kpoints = generate_kpoints_from_grid(kpt_cfg["mp_grid"])
    else:
        raise ValueError("k_points 配置缺少 mp_grid 或 kpoints。")

    kpath_file = kpt_cfg.get("kpath_file", "KPATH.in")
    band_segments: List[Dict[str, Any]] = []
    if os.path.exists(kpath_file):
        band_segments = parse_kpath(kpath_file)
        print(f"从 {kpath_file} 读取到 {len(band_segments)} 段能带路径。")
    else:
        print(f"未找到 {kpath_file}，将不写入 kpoint_path（可先用 vaspkit 303 生成）。")

    basename = cfg.get("system_name") or struct["comment"] or "wannier90"
    output = f"{basename}.win"

    write_win(output, struct, cfg, kpoints, band_segments)
    print(f"已生成 {output}，包含 {len(kpoints)} 个 k 点。")


if __name__ == "__main__":
    main()
