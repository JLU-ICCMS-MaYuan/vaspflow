#!/usr/bin/env python3
import argparse
import os
import re

import numpy as np


def is_float(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def find_first_file(directory, patterns):
    for pattern in patterns:
        matches = [f for f in os.listdir(directory) if re.fullmatch(pattern, f)]
        if matches:
            return os.path.join(directory, sorted(matches)[0])
    return None


def parse_tdos(path):
    energies = []
    tdos = []
    itdos = []
    fermi = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip().startswith("#"):
                match = re.search(r"EFermi\s*=\s*([-\d.]+)", line)
                if match:
                    fermi = float(match.group(1))
                continue
            parts = line.split()
            if len(parts) >= 3 and is_float(parts[0]):
                energies.append(float(parts[0]))
                tdos.append(float(parts[1]))
                itdos.append(float(parts[2]))
    if not energies:
        raise ValueError(f"未从 {path} 解析到 TDOS 数据。")
    return np.array(energies), np.array(tdos), np.array(itdos), fermi


def parse_projwfc_shells(path):
    if not path or not os.path.exists(path):
        return {}
    mapping = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            match = re.match(
                r"^\s*\d+\s+(\d+)\s+([A-Za-z]{1,2})\s+(\d+)([SPDF])\s+",
                line,
            )
            if not match:
                continue
            atom_index = int(match.group(1))
            shell_num = match.group(3)
            letter = match.group(4).lower()
            shell_label = f"{shell_num}{letter}"
            mapping.setdefault(atom_index, [])
            if shell_label not in mapping[atom_index]:
                mapping[atom_index].append(shell_label)
    return mapping


def parse_pdos_filename(path):
    name = os.path.basename(path)
    match = re.search(
        r"pdos_atm#(\d+)\(([A-Za-z]{1,2})\)_wfc#(\d+)\(([spdf])\)",
        name,
    )
    if not match:
        return None
    atom_index = int(match.group(1))
    element = match.group(2)
    wfc_index = int(match.group(3))
    l_letter = match.group(4).lower()
    return atom_index, element, wfc_index, l_letter


def orbital_components(l_letter):
    if l_letter == "s":
        return ["s"]
    if l_letter == "p":
        return ["pz", "px", "py"]
    if l_letter == "d":
        return ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]
    if l_letter == "f":
        return ["fz3", "fxz2", "fyz2", "fxyz", "fx2y2", "fy3x2-y2", "fzx2-y2"]
    return [l_letter]


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


def resolve_channel(channel, group_map, data, idx):
    if channel.endswith("__sum"):
        base = channel[:-5]
        return sum(data[label][idx] for label in group_map.get(base, []))
    return data[channel][idx]


def cumulative_integral(energy, values):
    integral = np.zeros_like(values)
    for i in range(1, len(values)):
        de = energy[i] - energy[i - 1]
        integral[i] = integral[i - 1] + 0.5 * (values[i] + values[i - 1]) * de
    return integral


def write_table(
    path,
    energy,
    channels,
    group_map,
    data,
    tdos_values,
    fermi,
    total_label,
):
    labels = ["Energy"] + [c.replace("__sum", "") for c in channels] + [total_label]
    width = max(12, max(len(label) for label in labels) + 2)
    fmt = f"{{:>{width}.6e}}"

    with open(path, "w", encoding="utf-8") as f:
        if fermi is not None:
            f.write(f"# EFermi = {fermi:.6f} eV\n")
        f.write("#" + "".join(f"{label:>{width}}" for label in labels) + "\n")
        for i in range(len(energy)):
            parts = [fmt.format(energy[i])]
            for channel in channels:
                parts.append(fmt.format(resolve_channel(channel, group_map, data, i)))
            parts.append(fmt.format(tdos_values[i]))
            f.write("".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser(description="处理 QE DOS/PDOS 数据输出表格。")
    parser.add_argument(
        "-i",
        "--input",
        default=".",
        help="TDOS 文件路径或包含 DOS 文件的目录。",
    )
    parser.add_argument(
        "--proj",
        default=None,
        help="projwfc 文件路径（用于解析壳层编号）。",
    )
    args = parser.parse_args()

    input_path = os.path.abspath(args.input)
    if os.path.isdir(input_path):
        work_dir = input_path
        tdos_path = find_first_file(work_dir, [r".+\.tdos"])
    else:
        work_dir = os.path.dirname(input_path)
        tdos_path = input_path if input_path.endswith(".tdos") else find_first_file(work_dir, [r".+\.tdos"])

    if not tdos_path or not os.path.exists(tdos_path):
        raise FileNotFoundError("未找到 *.tdos 文件。")

    proj_path = args.proj
    if not proj_path:
        proj_path = find_first_file(
            work_dir,
            [r".+\.projwfc_up", r".+\.projwfc_down", r".+\.projwfc", r".+\.proj"],
        )

    energy, tdos, itdos, fermi = parse_tdos(tdos_path)
    shell_map = parse_projwfc_shells(proj_path)

    pdos_files = sorted(
        f for f in os.listdir(work_dir) if re.search(r"pdos_atm#\d+\([^)]+\)_wfc#\d+\([spdf]\)", f)
    )
    if not pdos_files:
        raise FileNotFoundError("未找到 pdos_atm 文件。")

    data = {}
    for filename in pdos_files:
        file_path = os.path.join(work_dir, filename)
        parsed = parse_pdos_filename(file_path)
        if not parsed:
            continue
        atom_index, element, wfc_index, l_letter = parsed
        shell_label = None
        if atom_index in shell_map and wfc_index <= len(shell_map[atom_index]):
            shell_label = shell_map[atom_index][wfc_index - 1]
        if not shell_label:
            shell_label = f"{wfc_index}{l_letter}"

        shell_match = re.match(r"(\d+)([spdf])", shell_label)
        shell_num = shell_match.group(1) if shell_match else str(wfc_index)
        base_letter = shell_match.group(2) if shell_match else l_letter

        values = []
        energies = []
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 3 and is_float(parts[0]):
                    energies.append(float(parts[0]))
                    values.append([float(x) for x in parts[2:]])
        if not energies:
            continue

        energies = np.array(energies)
        if len(energies) != len(energy) or np.max(np.abs(energies - energy)) > 1e-6:
            raise ValueError(f"能量网格与 TDOS 不一致: {filename}")

        values = np.array(values)
        components = orbital_components(l_letter)
        if values.shape[1] < len(components):
            raise ValueError(f"{filename} 的 pdos 列数不足。")

        for idx, comp in enumerate(components):
            label = f"{element}{atom_index}-{shell_num}{comp}"
            data[label] = values[:, idx]

    if not data:
        raise ValueError("未解析到 PDOS 数据。")

    labels = list(data.keys())
    channels, group_map = build_channels(labels)
    sum_channels, sum_map, sum_group_map = build_sum_channels(labels)

    sum_data = {}
    for sum_label, members in sum_map.items():
        sum_data[sum_label] = np.zeros_like(energy)
        for label in members:
            sum_data[sum_label] += data[label]

    idos_data = {}
    for label, vals in data.items():
        idos_data[label] = cumulative_integral(energy, vals)

    idos_sum_data = {}
    for label, vals in sum_data.items():
        idos_sum_data[label] = cumulative_integral(energy, vals)

    base = os.path.basename(tdos_path)
    prefix = base.split(".")[0]

    dos_path = os.path.join(work_dir, f"{prefix}_dos_proj.dat")
    dos_sum_path = os.path.join(work_dir, f"{prefix}_dos_proj_sum.dat")
    idos_path = os.path.join(work_dir, f"{prefix}_idos_proj.dat")
    idos_sum_path = os.path.join(work_dir, f"{prefix}_idos_proj_sum.dat")

    write_table(dos_path, energy, channels, group_map, data, tdos, fermi, "TDOS")
    write_table(dos_sum_path, energy, sum_channels, sum_group_map, sum_data, tdos, fermi, "TDOS")
    write_table(idos_path, energy, channels, group_map, idos_data, itdos, fermi, "ITDOS")
    write_table(idos_sum_path, energy, sum_channels, sum_group_map, idos_sum_data, itdos, fermi, "ITDOS")

    print(f"已生成: {dos_path}")
    print(f"已生成: {dos_sum_path}")
    print(f"已生成: {idos_path}")
    print(f"已生成: {idos_sum_path}")


if __name__ == "__main__":
    main()
