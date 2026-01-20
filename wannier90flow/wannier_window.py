#!/usr/bin/env python3
import argparse
import xml.dom.minidom


def get_energies(xml_name):
    # VASP EIGENVAL format
    if xml_name == "EIGENVAL":
        with open("EIGENVAL", "r") as eigenval:
            raw_content = eigenval.readlines()
        content = raw_content[7:]
        entries = []
        for i, line in enumerate(content):
            if len(line.split()) == 4:
                entries.append(i)
        eng_full = []
        for i in entries:
            eng_ki = []
            for line in content[i+1:]:
                if len(line.split()) == 3:
                    eng_ki.append(float(line.split()[1]))
                else:
                    break
            eng_full.append(eng_ki)
    # Quantum ESPRESSO and FLEUR xml format
    else:
        Har2eV = 13.60569253 * 2
        dom = xml.dom.minidom.parse(xml_name)
        root = dom.documentElement
        eng_full = []
        if root.nodeName == "fleurOutput":
            eigenvalues = root.getElementsByTagName("eigenvalues")[-1]
            eks = eigenvalues.getElementsByTagName("eigenvaluesAt")
            eng_full = [[float(f) * Har2eV for f in ek.childNodes[0].data.split()] for ek in eks]
        elif root.nodeName == "qes:espresso":
            eigenvalues = root.getElementsByTagName("eigenvalues")
            eng_full = [[float(f) * Har2eV for f in ek.childNodes[0].data.split()] for ek in eigenvalues]
        else:
            raise RuntimeError("Unknown xml output")
    return eng_full


def build_parser():
    parser = argparse.ArgumentParser(
        description="解析能量窗口：-n 指定能带编号；-e 指定能量范围。"
    )
    parser.add_argument(
        "xml",
        help="能量文件：EIGENVAL 或 QE/FLEUR XML（含能带信息）",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-n",
        nargs="+",
        type=int,
        help="能带序号（1 或 2 个整数）。1 个：返回该带范围；2 个：返回 band1 最低与 band2 最高。",
    )
    group.add_argument(
        "-e",
        nargs="+",
        type=float,
        help="能量范围（1 或 2 个浮点）。2 个：统计区间内带数 min/max；1 个：视作 band_index。",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    eng_full = get_energies(args.xml)

    if args.n is not None:
        if len(args.n) not in (1, 2):
            parser.error("-n 需要 1 或 2 个整数参数")
        for idx in args.n:
            if idx < 1:
                parser.error("band_index 必须从 1 开始")
        if len(args.n) == 1:
            band_index = args.n[0]
            eng_selected = [eng[band_index - 1] for eng in eng_full]
            emin = min(eng_selected)
            emax = max(eng_selected)
            print(f"emin = {emin:.6f}")
            print(f"emax = {emax:.6f}")
        else:
            band_low, band_high = args.n
            eng_low = [eng[band_low - 1] for eng in eng_full]
            eng_high = [eng[band_high - 1] for eng in eng_full]
            emin = min(eng_low)
            emax = max(eng_high)
            print(f"emin = {emin:.6f}")
            print(f"emax = {emax:.6f}")
    else:
        if len(args.e) not in (1, 2):
            parser.error("-e 需要 1 或 2 个数值参数")
        if len(args.e) == 1:
            band_index = int(args.e[0])
            if band_index < 1:
                parser.error("band_index 必须从 1 开始")
            eng_selected = [eng[band_index - 1] for eng in eng_full]
            emin = min(eng_selected)
            emax = max(eng_selected)
            print(f"emin = {emin:.6f}")
            print(f"emax = {emax:.6f}")
        else:
            emin, emax = args.e
            nbnd = []
            for ik, ek in enumerate(eng_full):
                num_bands = 0
                for eng in ek:
                    if emin <= eng <= emax:
                        num_bands += 1
                print(f"ik = {ik + 1}, nbnd = {num_bands}")
                nbnd.append(num_bands)
            print(min(nbnd), max(nbnd))


if __name__ == "__main__":
    main()
