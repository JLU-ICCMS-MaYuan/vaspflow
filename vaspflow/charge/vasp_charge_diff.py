#!/usr/bin/env python3
import argparse
import os


def _read_vasp_charge(path):
    with open(path, "r") as f:
        lines = f.readlines()

    if len(lines) < 10:
        raise ValueError(f"文件内容过短，无法解析: {path}")

    def _is_int_tokens(tokens):
        if not tokens:
            return False
        for token in tokens:
            try:
                int(token)
            except ValueError:
                return False
        return True

    line5_tokens = lines[5].split()
    if _is_int_tokens(line5_tokens):
        counts_line_index = 5
        coord_line_index = 6
    else:
        counts_line_index = 6
        coord_line_index = 7

    atom_count_line = lines[counts_line_index].strip()
    coord_line = lines[coord_line_index].strip()
    if coord_line.lower().startswith("selective"):
        coord_line_index += 1
        coord_line = lines[coord_line_index].strip()

    coords_start = coord_line_index + 1

    try:
        atom_counts = [int(x) for x in atom_count_line.split()]
    except ValueError as e:
        raise ValueError(f"无法解析原子数行: {atom_count_line}") from e

    total_atoms = sum(atom_counts)
    coord_end = coords_start + total_atoms
    if coord_end >= len(lines):
        raise ValueError(f"坐标区长度不足: {path}")

    coords_block = lines[coords_start:coord_end]

    grid_line_index = coord_end
    grid_line = lines[grid_line_index].strip()
    while grid_line == "" and grid_line_index < len(lines):
        grid_line_index += 1
        grid_line = lines[grid_line_index].strip()

    grid_parts = grid_line.split()
    if len(grid_parts) != 3:
        raise ValueError(f"网格行格式异常: {grid_line}")

    try:
        ngx, ngy, ngz = (int(x) for x in grid_parts)
    except ValueError as e:
        raise ValueError(f"网格行无法解析为整数: {grid_line}") from e

    data_start = grid_line_index + 1
    data_lines = []
    for line in lines[data_start:]:
        if line.strip() == "":
            continue
        data_lines.extend(line.split())

    expected = ngx * ngy * ngz
    if len(data_lines) < expected:
        raise ValueError(f"电荷密度数据不足: {len(data_lines)} < {expected}")

    data = [float(x) for x in data_lines[:expected]]

    header_block = lines[:coords_start]

    return {
        "header": header_block,
        "coords": coords_block,
        "grid": (ngx, ngy, ngz),
        "data": data,
    }


def _write_vasp_charge(path, header, coords, grid, data):
    ngx, ngy, ngz = grid
    with open(path, "w") as f:
        for line in header:
            f.write(line if line.endswith("\n") else line + "\n")
        for line in coords:
            f.write(line if line.endswith("\n") else line + "\n")
        f.write("\n")
        f.write(f"  {ngx} {ngy} {ngz}\n")
        for i, value in enumerate(data, start=1):
            f.write(f"{value:12.4f}")
            if i % 10 == 0:
                f.write("\n")
        if len(data) % 10 != 0:
            f.write("\n")


def main():
    parser = argparse.ArgumentParser(description="VASP 电荷密度差分 (A - B)")
    parser.add_argument("a", help="电荷密度 A 文件路径 (PARCHG/CHGCAR)")
    parser.add_argument("b", help="电荷密度 B 文件路径 (PARCHG/CHGCAR)")
    parser.add_argument("-o", "--output", default="A_minus_B.vasp", help="输出文件名")
    args = parser.parse_args()

    a_data = _read_vasp_charge(args.a)
    b_data = _read_vasp_charge(args.b)

    if a_data["header"] != b_data["header"]:
        raise ValueError("A 与 B 的头信息不一致，无法相减")
    if a_data["coords"] != b_data["coords"]:
        raise ValueError("A 与 B 的原子坐标不一致，无法相减")
    if a_data["grid"] != b_data["grid"]:
        raise ValueError("A 与 B 的网格不一致，无法相减")

    diff = [av - bv for av, bv in zip(a_data["data"], b_data["data"])]
    _write_vasp_charge(args.output, a_data["header"], a_data["coords"], a_data["grid"], diff)

    print(f"已生成差分电荷密度: {args.output}")


if __name__ == "__main__":
    main()
