#!/usr/bin/env python3
import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np

def get_recip_base(cell):
    """从实空间基矢计算倒空间基矢 (2*pi*b)"""
    v1, v2, v3 = cell
    vol = np.dot(v1, np.cross(v2, v3))
    b1 = 2 * np.pi * np.cross(v2, v3) / vol
    b2 = 2 * np.pi * np.cross(v3, v1) / vol
    b3 = 2 * np.pi * np.cross(v1, v2) / vol
    return np.array([b1, b2, b3])

def pick_file(directory, patterns, desc):
    matches = []
    for pattern in patterns:
        matches.extend(glob.glob(os.path.join(directory, pattern)))
    matches = sorted(set(matches))
    if not matches:
        raise FileNotFoundError(f"找不到 {desc}，目录: {directory}")
    if len(matches) > 1:
        print(f"Warning: 找到多个{desc}，默认使用 {matches[0]}")
    return matches[0]


def resolve_inputs(qe_dir, w90_dir):
    qe_dir = os.path.abspath(qe_dir)
    w90_dir = os.path.abspath(w90_dir)

    qe_band_file = pick_file(
        qe_dir,
        ["*_band.dat"],
        "QE 能带数据",
    )
    w90_band_file = pick_file(
        w90_dir,
        ["*_band.dat"],
        "Wannier90 能带数据",
    )

    label_patterns = ["*_band.labelinfo.dat"]
    label_path = None
    try:
        label_path = pick_file(w90_dir, label_patterns, "高对称点信息")
    except FileNotFoundError:
        label_path = pick_file(qe_dir, label_patterns, "高对称点信息")

    return qe_band_file, w90_band_file, label_path


def parse_qe_plot_file(qe_band_file, b_basis, fermi_qe):
    qe_k_dist = []
    qe_bands = []

    with open(qe_band_file, "r") as f:
        lines = f.readlines()

    header = lines[0].strip()
    nbnd = None
    nks = None
    match_nbnd = re.search(r"nbnd\s*=\s*(\d+)", header)
    match_nks = re.search(r"nks\s*=\s*(\d+)", header)
    if match_nbnd and match_nks:
        nbnd = int(match_nbnd.group(1))
        nks = int(match_nks.group(1))
    else:
        parts = header.split()
        if len(parts) >= 6:
            try:
                nbnd = int(parts[2].replace(",", ""))
                nks = int(parts[5].replace(",", ""))
            except ValueError:
                nbnd = None
                nks = None

    if nbnd is None or nks is None:
        raise ValueError(f"无法解析 QE 能带文件头部: {lines[0].strip()}")

    cur_dist = 0.0
    prev_k_cart = None
    idx = 1
    for _ in range(nks):
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
        if idx >= len(lines):
            break
        line_parts = lines[idx].split()
        k_frac = np.array([float(x) for x in line_parts[:3]])
        k_cart = k_frac @ b_basis
        if prev_k_cart is not None:
            cur_dist += np.linalg.norm(k_cart - prev_k_cart)
        qe_k_dist.append(cur_dist)
        prev_k_cart = k_cart

        idx += 1
        ebands = []
        while len(ebands) < nbnd and idx < len(lines):
            ebands.extend([float(x) for x in lines[idx].split()])
            idx += 1
        qe_bands.append(ebands)

    qe_bands = np.array(qe_bands).T - fermi_qe
    return qe_k_dist, qe_bands


def parse_qe_gnu_file(qe_band_file, fermi_qe):
    x_vals = []
    y_vals = []
    with open(qe_band_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                if x_vals:
                    yield np.array(x_vals), np.array(y_vals) - fermi_qe
                x_vals = []
                y_vals = []
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            x_vals.append(float(parts[0]))
            y_vals.append(float(parts[1]))
    if x_vals:
        yield np.array(x_vals), np.array(y_vals) - fermi_qe


def parse_w90_band_file(w90_band_file, fermi_w90):
    segments = []
    x_vals = []
    y_vals = []
    with open(w90_band_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                if x_vals:
                    segments.append((np.array(x_vals), np.array(y_vals) - fermi_w90))
                x_vals = []
                y_vals = []
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            x_vals.append(float(parts[0]))
            y_vals.append(float(parts[1]))
    if x_vals:
        segments.append((np.array(x_vals), np.array(y_vals) - fermi_w90))
    return segments


def plot_comparison(qe_band_file, w90_band_file, w90_label_file, fermi_energy, output_img):
    # --- 1. 参数设置 ---
    # 晶胞参数 (来自 eband.in 或 Y1H6.win)
    cell = np.array([
        [-1.6848999334,  1.6848999334,  1.6848999334],
        [ 1.6848999334, -1.6848999334,  1.6848999334],
        [ 1.6848999334,  1.6848999334, -1.6848999334]
    ])
    
    # 费米能级 (来自 Y1H6.win)
    fermi_w90 = fermi_energy
    fermi_qe = fermi_energy 

    if not os.path.exists(qe_band_file):
        print(f"Error: 找不到 QE 数据文件 {qe_band_file}")
        return
    if not os.path.exists(w90_band_file):
        print(f"Error: 找不到 Wannier90 数据文件 {w90_band_file}")
        return

    # --- 2. 处理 QE 数据 ---
    b_basis = get_recip_base(cell)
    qe_is_gnu = qe_band_file.endswith(".gnu")
    if qe_is_gnu:
        qe_segments = list(parse_qe_gnu_file(qe_band_file, fermi_qe))
    else:
        qe_segments = parse_w90_band_file(qe_band_file, fermi_qe)

    # --- 3. 读取 Wannier90 数据 ---
    w90_segments = parse_w90_band_file(w90_band_file, fermi_w90)
    if not w90_segments:
        raise ValueError(f"Wannier90 能带文件为空或无法解析: {w90_band_file}")

    # --- 4. 绘图 ---
    plt.figure(figsize=(10, 7))

    # 绘制 Wannier90 能带 (蓝色线，底层)
    for idx, (w90_k, w90_e) in enumerate(w90_segments):
        label = 'Wannier90' if idx == 0 else ""
        plt.plot(w90_k, w90_e, 'b-', linewidth=1.5, alpha=0.8, label=label, zorder=2)

    # 绘制 QE 能带 (红色虚线，上层)
    for idx, (x_vals, y_vals) in enumerate(qe_segments):
        label = 'DFT (QE)' if idx == 0 else ""
        plt.plot(x_vals, y_vals, 'r--', linewidth=1.2, alpha=0.8, label=label, zorder=3)

    # 绘制高对称点垂直线和标签
    if w90_label_file and os.path.exists(w90_label_file):
        with open(w90_label_file, 'r') as f:
            tick_coords = []
            tick_labels = []
            for line in f:
                parts = line.split()
                if not parts: continue
                label = parts[0].replace('GAMMA', r'$\Gamma$')
                dist = float(parts[2])
                plt.axvline(x=dist, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
                tick_coords.append(dist)
                tick_labels.append(label)
            plt.xticks(tick_coords, tick_labels)

    plt.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    plt.ylabel('Energy - $E_f$ (eV)')
    plt.title('Band Structure Comparison: QE (Dots) vs Wannier90 (Lines)')
    
    # 自动去重 legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right')

    plt.xlim(0, max(seg[0].max() for seg in qe_segments))
    w90_min = min(seg[1].min() for seg in w90_segments)
    w90_max = max(seg[1].max() for seg in w90_segments)
    plt.ylim(w90_min, w90_max * 1.2)
    plt.grid(True, axis='y', linestyle=':', alpha=0.4)
    
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Comparison plot saved to {output_img}")

def main():
    parser = argparse.ArgumentParser(description='Compare QE and Wannier90 band structures.')
    parser.add_argument('--qe', default='.', help='QE 数据目录')
    parser.add_argument('--w90', default='.', help='Wannier90 数据目录')
    parser.add_argument('--fermi', type=float, default=23.3313, help='Fermi energy to align (eV)')
    parser.add_argument('--out', default='band_comparison.png', help='Output image filename')
    
    args = parser.parse_args()
    qe_band_file, w90_band_file, label_file = resolve_inputs(args.qe, args.w90)
    plot_comparison(qe_band_file, w90_band_file, label_file, args.fermi, args.out)

if __name__ == "__main__":
    main()
