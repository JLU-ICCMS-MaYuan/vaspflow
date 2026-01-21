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
        description="解析能量窗口：-n 指定能带编号；-e 指定能量范围；--auto 自动估计能窗。"
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
    group.add_argument(
        "--auto",
        nargs=3,
        type=int,
        metavar=("nbnd1", "num_wann", "nbnd3"),
        help="自动估计能窗：nbnd1 起始能带，num_wann 目标带数，nbnd3 用于 dis_win_max 的偏移带数。",
    )
    parser.add_argument(
        "--auto-delta",
        type=float,
        default=0.02,
        help="自动能窗边界偏移（默认 0.02 eV）。",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    eng_full = get_energies(args.xml)

    def plot_band_ranges(
        start_band: int,
        end_band: int,
        dis_froz_min: float,
        dis_froz_max: float,
        dis_win_min: float,
        dis_win_max: float,
        output_path: str,
    ) -> None:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("未安装 matplotlib，跳过能带范围绘图。")
            return

        band_indices = list(range(start_band, end_band + 1))
        band_mins = [min(eng[i - 1] for eng in eng_full) for i in band_indices]
        band_maxs = [max(eng[i - 1] for eng in eng_full) for i in band_indices]

        fig, ax = plt.subplots(figsize=(8, 6))
        for idx, emin, emax in zip(band_indices, band_mins, band_maxs):
            ax.bar(idx, emax - emin, bottom=emin, width=0.6, color="red", alpha=0.3, edgecolor="red")
        ax.axhline(dis_froz_min, color="red", linestyle="--", linewidth=1.2, label="dis_froz_min")
        ax.axhline(dis_froz_max, color="red", linestyle="--", linewidth=1.2, label="dis_froz_max")
        ax.axhline(dis_win_min, color="black", linestyle="--", linewidth=1.2, label="dis_win_min")
        ax.axhline(dis_win_max, color="black", linestyle="--", linewidth=1.2, label="dis_win_max")
        ax.set_xlabel("Band Index")
        ax.set_ylabel("Energy")
        ax.set_title("Energy spectrum of bands")
        ax.set_xticks(band_indices)
        ax.legend(frameon=False, loc="upper left")
        fig.tight_layout()
        fig.savefig(output_path, dpi=300)
        plt.close(fig)

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
    elif args.e is not None:
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
    else:
        nbnd1, num_wann, nbnd3 = args.auto
        if nbnd1 < 1 or num_wann < 1 or nbnd3 < 0:
            parser.error("--auto 参数必须满足：nbnd1>=1, num_wann>=1, nbnd3>=0")

        if nbnd1 > len(eng_full[0]):
            parser.error("nbnd1 超出能带数量")

        delta = args.auto_delta
        dis_froz_min = min(eng[nbnd1 - 1] for eng in eng_full) - delta
        dis_win_min = dis_froz_min

        nbnd2 = nbnd1 + num_wann
        if nbnd2 > len(eng_full[0]):
            parser.error("nbnd1 + num_wann 超出能带数量")

        def max_band_count(emin, emax):
            max_count = 0
            for ek in eng_full:
                count = sum(1 for eng in ek if emin <= eng <= emax)
                if count > max_count:
                    max_count = count
            return max_count

        while True:
            if nbnd2 < nbnd1:
                parser.error("自适应失败：nbnd2 已小于 nbnd1")
            dis_froz_max = max(eng[nbnd2 - 1] for eng in eng_full) + delta
            max_count = max_band_count(dis_froz_min, dis_froz_max)
            print(f"auto: nbnd2={nbnd2}, dis_froz_max={dis_froz_max:.6f}, max_nbnd={max_count}")
            if max_count <= num_wann:
                break
            nbnd2 -= 1

        nbnd_win = nbnd2 + nbnd3
        if nbnd_win > len(eng_full[0]):
            parser.error("nbnd2 + nbnd3 超出能带数量")
        dis_win_max = max(eng[nbnd_win - 1] for eng in eng_full) + delta

        print(f"dis_froz_min = {dis_froz_min:.6f}")
        print(f"dis_froz_max = {dis_froz_max:.6f}")
        print(f"dis_win_min = {dis_win_min:.6f}")
        print(f"dis_win_max = {dis_win_max:.6f}")
        plot_end = nbnd_win + 8
        if plot_end > len(eng_full[0]):
            plot_end = len(eng_full[0])
        plot_band_ranges(
            nbnd1,
            plot_end,
            dis_froz_min,
            dis_froz_max,
            dis_win_min,
            dis_win_max,
            "band_ranges.png",
        )


if __name__ == "__main__":
    main()
