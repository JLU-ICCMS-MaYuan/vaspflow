#!/usr/bin/env python3
import os
import argparse
import re
import numpy as np


def is_float(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def guess_prefix_from_filename(path):
    name = os.path.basename(path)
    for suffix in ["_band.projwfc_up", "_band.projwfc_down", "_band.projwfc", "_band.dat"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if "." in name:
        return name.split(".")[0]
    return name


def orbital_from_lm(l_val, m_val):
    if l_val == 0:
        return "s"
    if l_val == 1:
        p_map = {1: "pz", 2: "px", 3: "py"}
        return p_map.get(m_val, f"p{m_val}")
    if l_val == 2:
        d_map = {1: "dz2", 2: "dxz", 3: "dyz", 4: "dx2-y2", 5: "dxy"}
        return d_map.get(m_val, f"d{m_val}")
    if l_val == 3:
        f_map = {
            1: "fz3",
            2: "fxz2",
            3: "fyz2",
            4: "fxyz",
            5: "fx2y2",
            6: "fy3x2-y2",
            7: "fzx2-y2",
        }
        return f_map.get(m_val, f"f{m_val}")
    return f"l{l_val}m{m_val}"


def read_cell_from_eband(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到 eband.in: {path}")
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        if line.strip().upper().startswith("CELL_PARAMETERS"):
            cell = []
            for j in range(1, 4):
                if idx + j >= len(lines):
                    raise ValueError("CELL_PARAMETERS 段不完整。")
                parts = lines[idx + j].split()
                if len(parts) < 3:
                    raise ValueError("CELL_PARAMETERS 行格式错误。")
                cell.append([float(parts[0]), float(parts[1]), float(parts[2])])
            return np.array(cell, dtype=float)
    raise ValueError("eband.in 中未找到 CELL_PARAMETERS 段。")


def read_kpath_points_from_eband(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到 eband.in: {path}")
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    kpoint_line = None
    for idx, line in enumerate(lines):
        if line.strip().upper().startswith("K_POINTS"):
            kpoint_line = idx
            break
    if kpoint_line is None:
        raise ValueError("eband.in 中未找到 K_POINTS 段。")
    if "CRYSTAL_B" not in lines[kpoint_line].upper():
        raise ValueError("K_POINTS 段不是 crystal_b 格式。")

    idx = kpoint_line + 1
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    if idx >= len(lines):
        raise ValueError("K_POINTS 段缺少点数。")
    try:
        num_points = int(lines[idx].split()[0])
    except ValueError as exc:
        raise ValueError("K_POINTS 点数解析失败。") from exc
    idx += 1

    points = []
    labels = []
    for _ in range(num_points):
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
        if idx >= len(lines):
            raise ValueError("K_POINTS 段数据不足。")
        line = lines[idx]
        idx += 1
        if "!" in line:
            coord_part, label_part = line.split("!", 1)
            label = label_part.strip()
        else:
            coord_part = line
            label = ""
        parts = coord_part.split()
        if len(parts) < 3:
            raise ValueError(f"K_POINTS 行格式错误: {line}")
        points.append([float(parts[0]), float(parts[1]), float(parts[2])])
        labels.append(label)
    return np.array(points, dtype=float), labels


def read_kpoints_weights_from_eband(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到 eband.in: {path}")
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    kpoint_line = None
    for idx, line in enumerate(lines):
        if line.strip().upper().startswith("K_POINTS"):
            kpoint_line = idx
            break
    if kpoint_line is None:
        raise ValueError("eband.in 中未找到 K_POINTS 段。")
    if "CRYSTAL_B" not in lines[kpoint_line].upper():
        raise ValueError("K_POINTS 段不是 crystal_b 格式。")

    idx = kpoint_line + 1
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    if idx >= len(lines):
        raise ValueError("K_POINTS 段缺少点数。")
    try:
        num_points = int(lines[idx].split()[0])
    except ValueError as exc:
        raise ValueError("K_POINTS 点数解析失败。") from exc
    idx += 1

    weights = []
    for _ in range(num_points):
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
        if idx >= len(lines):
            raise ValueError("K_POINTS 段数据不足。")
        line = lines[idx]
        idx += 1
        coord_part = line.split("!", 1)[0]
        parts = coord_part.split()
        if len(parts) < 4:
            raise ValueError(f"K_POINTS 行格式错误: {line}")
        try:
            weight = int(float(parts[3]))
        except ValueError as exc:
            raise ValueError(f"K_POINTS 权重解析失败: {line}") from exc
        weights.append(weight)
    return weights


def get_recip_lattice(cell):
    vol = np.dot(cell[0], np.cross(cell[1], cell[2]))
    b1 = 2 * np.pi * np.cross(cell[1], cell[2]) / vol
    b2 = 2 * np.pi * np.cross(cell[2], cell[0]) / vol
    b3 = 2 * np.pi * np.cross(cell[0], cell[1]) / vol
    return np.array([b1, b2, b3], dtype=float)


def compute_kpath_dist(coords, cell, weights=None):
    recip = get_recip_lattice(cell)
    metric = recip @ recip.T
    dist = [0.0]
    jump_points = set()
    if weights:
        expected_1 = sum(weights) + 1
        expected_2 = sum(weights)
        if len(coords) in (expected_1, expected_2):
            indices = [1]
            for i in range(1, len(weights)):
                indices.append(indices[-1] + int(weights[i - 1]))
            jump_indices = set(indices[i] for i in range(1, len(indices)) if weights[i - 1] == 1)
            jump_points = set(jump_indices)
        else:
            weights = None

    if not weights:
        steps = []
        for i in range(1, len(coords)):
            delta = np.array(coords[i]) - np.array(coords[i - 1])
            steps.append(float(np.sqrt(delta @ metric @ delta)))
        nonzero_steps = [s for s in steps if s > 1.0e-10]
        median_step = np.median(nonzero_steps) if nonzero_steps else 0.0
        jump_threshold = median_step * 5.0 if median_step > 0 else float("inf")

        for i in range(1, len(coords)):
            delta = np.array(coords[i]) - np.array(coords[i - 1])
            step = float(np.sqrt(delta @ metric @ delta))
            if step > jump_threshold:
                dist.append(dist[-1])
                continue
            dist.append(dist[-1] + step)
    else:
        for i in range(1, len(coords)):
            if i + 1 in jump_points:
                dist.append(dist[-1])
                continue
            delta = np.array(coords[i]) - np.array(coords[i - 1])
            dist.append(dist[-1] + float(np.sqrt(delta @ metric @ delta)))
    return dist


def normalize_band_prefix(prefix):
    if prefix.endswith("_band_proj"):
        return prefix[: -len("_band_proj")]
    if prefix.endswith("_proj"):
        return prefix[: -len("_proj")]
    return prefix




def build_channels(labels):
    group_map = {}
    for label in labels:
        parsed = parse_label(label)
        if not parsed:
            continue
        elem, atom_index, shell_num, base_letter, _suffix = parsed
        if base_letter not in ("p", "d", "f"):
            continue
        group = f"{elem}{atom_index}-{shell_num}{base_letter}"
        group_map.setdefault(group, []).append(label)

    channels = order_channels(labels, group_map)
    return channels, group_map


def build_sum_channels(labels):
    sum_map = {}
    ordered = []
    for label in labels:
        parsed = parse_label(label)
        if not parsed:
            continue
        elem, _atom_index, shell_num, base_letter, suffix = parsed
        suffix = suffix or ""
        sum_label = f"{elem}-{shell_num}{base_letter}{suffix}"
        if sum_label not in sum_map:
            sum_map[sum_label] = []
            ordered.append(sum_label)
        sum_map[sum_label].append(label)

    group_map = {}
    for label in ordered:
        parsed = parse_label(label)
        if not parsed:
            continue
        elem, _atom_index, shell_num, base_letter, _suffix = parsed
        if base_letter not in ("p", "d", "f"):
            continue
        group = f"{elem}-{shell_num}{base_letter}"
        group_map.setdefault(group, []).append(label)

    channels = order_channels(ordered, group_map)
    return channels, sum_map, group_map


def parse_label(label):
    match = re.match(r"^([A-Z][a-z]?)(\d+)?-(\d+)([spdf])(.+)?$", label)
    if not match:
        return None
    elem, atom_index, shell_num, base_letter, suffix = match.groups()
    atom_index = int(atom_index) if atom_index else 0
    return elem, atom_index, int(shell_num), base_letter, (suffix or "")


def order_channels(labels, group_map):
    suffix_orders = {
        "s": [""],
        "p": ["x", "y", "z"],
        "d": ["z2", "xz", "yz", "x2-y2", "xy"],
        "f": ["z3", "xz2", "yz2", "xyz", "x2y2", "y3x2-y2", "zx2-y2"],
    }
    base_order = ["s", "p", "d", "f"]

    grouped = {}
    for label in labels:
        parsed = parse_label(label)
        if not parsed:
            grouped.setdefault(label, []).append(label)
            continue
        elem, atom_index, shell_num, base_letter, suffix = parsed
        grouped.setdefault((elem, atom_index, shell_num, base_letter), []).append((suffix, label))

    ordered = []
    for key in sorted(
        [k for k in grouped if isinstance(k, tuple)],
        key=lambda x: (x[0], x[1], x[2], base_order.index(x[3])),
    ):
        elem, atom_index, shell_num, base_letter = key
        suffix_labels = {suffix: lbl for suffix, lbl in grouped[key]}
        for suffix in suffix_orders.get(base_letter, []):
            label = suffix_labels.get(suffix)
            if label:
                ordered.append(label)
        if atom_index == 0:
            group_label = f"{elem}-{shell_num}{base_letter}"
        else:
            group_label = f"{elem}{atom_index}-{shell_num}{base_letter}"
        if group_label in group_map:
            ordered.append(f"{group_label}__sum")

    for key in grouped:
        if not isinstance(key, tuple):
            ordered.extend(grouped[key])
    return ordered


def resolve_weight(channel, group_map, proj_index, proj_weights, k_idx, band_idx):
    if channel.endswith("__sum"):
        base = channel[:-5]
        weight = 0.0
        for label in group_map.get(base, []):
            weight += proj_weights[proj_index[label], k_idx, band_idx]
        return weight
    return proj_weights[proj_index[channel], k_idx, band_idx]


def parse_elebanddata(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到能带数据文件: {path}")

    with open(path, "r", encoding="utf-8") as f:
        lines = [line.rstrip() for line in f if line.strip() != ""]

    header = lines[0]
    match = re.search(r"nbnd\s*=\s*(\d+).+nks\s*=\s*(\d+)", header)
    if not match:
        raise ValueError("无法从 elebanddata.dat 解析 nbnd/nks。")
    nbnd = int(match.group(1))
    nks = int(match.group(2))

    coords = []
    energies = []
    idx = 1
    while idx < len(lines) and len(coords) < nks:
        parts = lines[idx].split()
        if len(parts) >= 3 and all(is_float(p) for p in parts[:3]):
            coords.append([float(parts[0]), float(parts[1]), float(parts[2])])
            idx += 1
            band_vals = []
            while idx < len(lines) and len(band_vals) < nbnd:
                vals = [float(v) for v in lines[idx].split() if is_float(v)]
                band_vals.extend(vals)
                idx += 1
            energies.append(band_vals[:nbnd])
        else:
            idx += 1

    if len(coords) != nks or len(energies) != nks:
        raise ValueError("elebanddata.dat 解析的 k 点数量与头部不一致。")

    return np.array(coords), np.array(energies), nbnd, nks


def parse_projwfc_file(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = [line.rstrip() for line in f]

    header_idx = None
    nproj = nk = nbnd = None
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 3 and all(p.isdigit() for p in parts):
            if i + 1 < len(lines) and re.match(r"^[FT](\s+[FT])*$", lines[i + 1].strip(), re.IGNORECASE):
                nproj, nk, nbnd = map(int, parts)
                header_idx = i + 2
                break
    if header_idx is None:
        return None

    labels = []
    proj_weights = np.zeros((nproj, nk, nbnd))

    idx = header_idx
    for proj_idx in range(nproj):
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
        if idx >= len(lines):
            raise ValueError("projwfc 文件结构不完整，缺少投影头部。")
        tokens = lines[idx].split()
        idx += 1
        if len(tokens) < 6:
            raise ValueError(f"投影头部解析失败: {lines[idx-1]}")

        atom_index = tokens[1]
        element = tokens[2]
        shell_token = None
        for token in tokens:
            if re.search(r"[SPDFspdf]$", token) and re.search(r"\d", token):
                shell_token = token
                break
        trailing_ints = [int(t) for t in tokens if t.isdigit()]
        if len(trailing_ints) < 2:
            raise ValueError(f"无法解析 l/m: {lines[idx-1]}")
        l_val = trailing_ints[-2]
        m_val = trailing_ints[-1]

        orbital = orbital_from_lm(l_val, m_val)
        if shell_token:
            shell_label = shell_token.strip().lower()
            base_letter = shell_label[-1]
            if orbital.startswith(base_letter):
                suffix = orbital[1:]
            else:
                suffix = orbital
            label = f"{element}{atom_index}-{shell_label}{suffix}"
        else:
            label = f"{element}{atom_index}-{orbital}"
        labels.append(label)

        total_lines = nk * nbnd
        for _ in range(total_lines):
            if idx >= len(lines):
                raise ValueError("projwfc 文件数据长度不足。")
            parts = lines[idx].split()
            idx += 1
            if len(parts) < 3:
                continue
            k_idx = int(parts[0]) - 1
            b_idx = int(parts[1]) - 1
            proj_weights[proj_idx, k_idx, b_idx] = float(parts[2])

    return labels, proj_weights, nk, nbnd


def main():
    parser = argparse.ArgumentParser(description="Process QE projwfc band projections.")
    parser.add_argument("-i", "--input", required=True, help="projwfc 投影文件（如 *.projwfc_up）")
    parser.add_argument("-b", "--bandfile", default="elebanddata.dat", help="bands.x 输出文件")
    parser.add_argument("--cell", default="eband.in", help="包含 CELL_PARAMETERS 的 eband.in")
    parser.add_argument("-o", "--output", help="输出文件名，默认 <prefix>_eleband_proj.dat")
    parser.add_argument("--prefix", help="输出文件前缀（覆盖自动推断）")
    args = parser.parse_args()

    input_path = args.input
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"找不到输入文件: {input_path}")

    proj_result = parse_projwfc_file(input_path)
    if proj_result is None:
        raise ValueError("未识别 projwfc 投影文件格式。")

    proj_labels, proj_weights, nk, nbnd = proj_result
    coords, energies, nbnd_band, nks_band = parse_elebanddata(args.bandfile)

    if nk != nks_band or nbnd != nbnd_band:
        raise ValueError(f"投影文件与能带文件尺寸不一致: projwfc(nk={nk}, nbnd={nbnd}) vs bands(nk={nks_band}, nbnd={nbnd_band})")

    cell = read_cell_from_eband(args.cell)
    kpath_points, kpath_labels = read_kpath_points_from_eband(args.cell)
    weights = read_kpoints_weights_from_eband(args.cell)
    kpath_dist = compute_kpath_dist(coords, cell, weights)
    channels, group_map = build_channels(proj_labels)
    proj_index = {label: idx for idx, label in enumerate(proj_labels)}

    prefix = args.prefix or guess_prefix_from_filename(input_path)
    output_path = args.output or f"{prefix}_eleband_proj.dat"
    sum_output_path = f"{prefix}_eleband_proj_sum.dat"
    band_prefix = normalize_band_prefix(prefix)
    band_output_path = f"{band_prefix}_band.dat"

    col_width = max(12, max(len(name) for name in channels) if channels else 12)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"# source = {input_path}\n")
        f.write(f"# bands = {args.bandfile}\n")
        f.write(f"# cell = {args.cell}\n")
        header_parts = [f"{'kpath_dist':>12}", f"{'Energy(eV)':>14}"]
        header_parts.extend(f"{name.replace('__sum',''):>{col_width}}" for name in channels)
        f.write(" ".join(header_parts) + "\n")
        for band_idx in range(nbnd):
            for k_idx in range(nk):
                energy = energies[k_idx, band_idx]
                weights = []
                for channel in channels:
                    weight = resolve_weight(channel, group_map, proj_index, proj_weights, k_idx, band_idx)
                    weights.append(f"{weight:{col_width}.6f}")
                row = f"{kpath_dist[k_idx]:12.6f} {energy:14.6f} " + " ".join(weights)
                f.write(row + "\n")
            f.write(f"#band={band_idx + 1}\n\n")

    sum_channels, sum_map, sum_group_map = build_sum_channels(proj_labels)
    sum_col_width = max(12, max(len(name) for name in sum_channels) if sum_channels else 12)

    with open(sum_output_path, "w", encoding="utf-8") as f:
        f.write(f"# source = {input_path}\n")
        f.write(f"# bands = {args.bandfile}\n")
        f.write(f"# cell = {args.cell}\n")
        f.write("# summed by element and shell (e.g., Y-4s)\n")
        header_parts = [f"{'kpath_dist':>12}", f"{'Energy(eV)':>14}"]
        header_parts.extend(f"{name.replace('__sum',''):>{sum_col_width}}" for name in sum_channels)
        f.write(" ".join(header_parts) + "\n")
        for band_idx in range(nbnd):
            for k_idx in range(nk):
                energy = energies[k_idx, band_idx]
                weights = []
                for channel in sum_channels:
                    if channel.endswith("__sum"):
                        base = channel[:-5]
                        weight = 0.0
                        for label in sum_group_map.get(base, []):
                            for base_label in sum_map[label]:
                                weight += proj_weights[proj_index[base_label], k_idx, band_idx]
                    else:
                        weight = 0.0
                        if channel in sum_map:
                            for base_label in sum_map[channel]:
                                weight += proj_weights[proj_index[base_label], k_idx, band_idx]
                    weights.append(f"{weight:{sum_col_width}.6f}")
                row = f"{kpath_dist[k_idx]:12.6f} {energy:14.6f} " + " ".join(weights)
                f.write(row + "\n")
            f.write(f"#band={band_idx + 1}\n\n")

    with open(band_output_path, "w", encoding="utf-8") as f:
        f.write(f"# source = {args.bandfile}\n")
        f.write(f"# cell = {args.cell}\n")
        f.write("# kpath_dist energy\n")
        for band_idx in range(nbnd):
            for k_idx in range(nk):
                energy = energies[k_idx, band_idx]
                f.write(f"{kpath_dist[k_idx]:12.6f} {energy:14.6f}\n")
            f.write(f"#band={band_idx + 1}\n\n")

    print(f"已生成投影能带数据: {output_path}")
    print(f"已生成元素汇总投影能带数据: {sum_output_path}")
    print(f"已生成能带数据: {band_output_path}")


if __name__ == "__main__":
    main()
