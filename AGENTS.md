# DFTFlow AGENTS.md

## 需求处理记录

### [2026-01-20] QE K 点配置统一至 [k_points]
- **需求**：
    1. 删除 `[nscf]` 中的 K 点配置，统一在 `[k_points]` 管理 NSCF/SCF 的 K 点设置。
    2. `qe_scf.py` 需支持两种方式：`kmesh` 自动生成网格，或 `kpoints` 显式网格（SCF 直接写入 `K_POINTS automatic`）。
- **方案**：
    1. 更新 `qeflow/nscf/qe_nscf.py`，仅从 `[k_points]` 读取 `kpoints` 或 `kmesh`（`wan` 同步迁移），未配置时使用默认间距生成均匀网格。
    2. 更新 `qeflow/scf/qe_scf.py`，支持 `[k_points].kpoints` 优先，其次 `kmesh` 生成自动网格，不写显式坐标。
    3. 调整 `qeflow/input.toml`，将示例字段改为 `[k_points].kpoints`/`kmesh`，保留 `[nscf_params]` 供 NSCF 参数覆盖。
- **状态**：已完成。

### [2026-01-20] Wannier90 输入准备脚本
- **需求**：
    1. 提供脚本在 `wannier90flow` 下自动生成 `.win` 初始输入。
    2. 规划 Wannier90 使用步骤，简化流程。
- **方案**：
    1. 新增 `wannier90flow/wannier_init.py`：读取 POSCAR 与 TOML/JSON 配置，生成包含 `unit_cell_cart`、`atoms_*`、`projections`、`mp_grid`、显式 `kpoints` 及可选 `kpoint_path` 的 `.win` 文件（`kpoint_path` 通过 vaspkit 303 生成的 `KPATH.in` 自动读取）。
    2. 添加示例配置 (`wannier90flow/sample_wannier90.toml`，另有 `wannier90flow/input.toml` 可直接使用)，复用 Ce1Sc2H24 参数，演示投影、能窗与网格设置。
    3. 在 `pyproject.toml` 注册脚本入口 `wannier90_prepare`（指向 `wannier90flow.wannier_init:main`）。
- **状态**：已完成。

### [2026-01-20] Wannier90 工具参数化
- **需求**：`wannier90flow/wannier_window.py` 使用 argparse 管理输入参数。
- **方案**：改为 argparse 解析 xml 文件、模式（e/n）及对应 band_index 或能量区间。
- **状态**：已完成。

### [2026-01-20] Wannier90 能量窗口参数统一为子命令
- **需求**：`wannier90flow/wannier_window.py` 所有输入参数统一由 argparse 解析，并按模式拆分参数。
- **方案**：使用 `argparse` 子命令 `e`/`n` 分别接收 `band_index` 与 `emin/emax`，保留原有调用结构。
- **状态**：已完成。

### [2026-01-20] Wannier90 脚本入口补齐
- **需求**：在 `pyproject.toml` 中补齐 wannier 相关脚本入口。
- **方案**：在 `[project.scripts]` 新增 `wannier_window` 指向 `wannier90flow.wannier_window:main`，保留现有 `wannier_init`。
- **状态**：已完成。

### [2026-01-20] Wannier90 可执行路径与 --run 支持
- **需求**：在测试配置中补充 Wannier90 与 pw2wannier90 的可执行文件路径，并支持 `--run` 执行流程。
- **方案**：
    1. 在 `tests/inputwannier.toml` 增加 `[wannier90]` 与 `[pw2wannier90]` 的 `executable_path`。
    2. 在 `wannier90flow/wannier_init.py` 增加 `--run`，顺序执行 `wannier90 -pp`、`pw2wannier90`、`wannier90`。
    3. 运行脚本中加入 `pw2wannier90` 输入文件检查。
- **状态**：已完成。

### [2026-01-20] Wannier90 配置归并与 Slurm 脚本准备
- **需求**：将 .win 相关参数归并到 `[win]`，并增加 Slurm 脚本生成。
- **方案**：
    1. `wannier90flow/input.toml` 与 `tests/inputwannier.toml` 的 .win 参数移动至 `[win]`，读取逻辑兼容旧顶层参数。
    2. 在 `wannier90flow/wannier_init.py` 生成 `run_wannier90.slurm`，调用 `run_wannier90.sh` 执行流程。
- **状态**：已完成。

### [2026-01-20] Wannier90 bands_plot 与 kpoint_path 联动
- **需求**：避免在未提供 `kpoint_path` 时仍写入 `bands_plot = .true.` 导致 wannier90 报错。
- **方案**：
    1. `wannier90flow/wannier_init.py` 中当 `bands_plot` 为 true 且未检测到 `kpoint_path` 时自动关闭，并提示。
    2. `tests/inputwannier.toml` 默认关闭 `bands_plot`，避免无路径时报错。
- **状态**：已完成。

### [2026-01-20] pw2wannier90 输入文件自动生成
- **需求**：`--run` 时缺少 `seedname.pw2wan` 造成流程中断。
- **方案**：
    1. `wannier90flow/wannier_init.py` 自动解析/生成 `pw2wannier90` 输入文件（含 `outdir/prefix/seedname`），并兼容既有 `input_file`。
    2. `tests/inputwannier.toml` 补充 `pw2wannier90` 的 `input_file/prefix/outdir/seedname`。
- **状态**：已完成。

### [2026-01-20] Wannier90 运行脚本合并
- **需求**：合并 `run_wannier90.sh` 与 `run_wannier90.slurm`，仅保留一个带 Slurm 头的脚本。
- **方案**：`wannier90flow/wannier_init.py` 仅生成 `run_wannier90.sh`，内置默认 Slurm 头，不再生成 `.slurm` 文件。
- **状态**：已完成。

### [2026-01-20] Wannier90 运行脚本简化
- **需求**：`run_wannier90.sh` 按模板简化，直接执行命令不再定义变量。
- **方案**：在 `wannier90flow/wannier_init.py` 中直接写入 `wannier90.x` 与 `pw2wannier90.x` 命令行，并保留日志重定向。
- **状态**：已完成。

### [2026-01-20] Wannier90 本地运行打印命令
- **需求**：本地 `--run` 时打印每一步实际执行的命令，便于排查。
- **方案**：在 `wannier90flow/wannier_init.py` 中为三步命令构造并打印，再执行。
- **状态**：已完成。

### [2026-01-20] Wannier90 能量窗口参数增强
- **需求**：`wannier90flow/wannier_window.py` 支持 `-n`/`-e` 复合参数：能带范围与能量区间统计。
- **方案**：
    1. `-n` 接受 1 或 2 个能带编号，输出对应能量范围。
    2. `-e` 接受 1 或 2 个数值；2 个按能量区间统计 min/max 带数，1 个按旧行为视作 band_index。
- **状态**：已完成。

### [2026-01-20] Wannier90 自动估计能窗
- **需求**：新增 `--auto nbnd1 num_wann nbnd3` 自动确定 `dis_froz_min/max` 与 `dis_win_min/max`。
- **方案**：
    1. `dis_froz_min = dis_win_min = band(nbnd1) 最低值`。
    2. 自适应下调 `nbnd2 = nbnd1 + num_wann` 的上界，直到区间最大带数等于 `num_wann`，否则报错。
    3. `dis_win_max = band(nbnd2 + nbnd3) 最大值`。
- **状态**：已完成。

### [2026-01-21] Wannier90 自动能窗迭代打印
- **需求**：`--auto` 过程中打印每次 `nbnd2` 与对应最大带数，便于排查。
- **方案**：在自适应循环中输出 `nbnd2/dis_froz_max/max_nbnd`。
- **状态**：已完成。

### [2026-01-21] Wannier90 自动能窗边界偏移
- **需求**：`--auto` 计算的 `dis_froz_min/dis_win_min` 使用最小值减 0.02，`dis_froz_max/dis_win_max` 使用最大值加 0.02。
- **方案**：在 `wannier90flow/wannier_window.py` 中对自动能窗边界增加 ±0.02 偏移。
- **状态**：已完成。

### [2026-01-21] Wannier90 自动能窗偏移可配置
- **需求**：支持 `--auto-delta` 配置自动能窗边界偏移，默认 0.02。
- **方案**：新增 `--auto-delta` 参数并应用到 `dis_froz_min/max` 与 `dis_win_min/max`。
- **状态**：已完成。

### [2026-01-21] Wannier90 自动能窗绘图
- **需求**：`--auto` 完成后绘制能带能量范围图，并保存为 `band_ranges.png`。
- **方案**：在 `wannier90flow/wannier_window.py` 中使用 matplotlib 绘制 `nbnd1..(nbnd_win+8)` 的能带能量区间柱状图，并标出 `dis_froz_min/max` 与 `dis_win_min/max` 虚线，若缺少 matplotlib 则提示跳过。
- **状态**：已完成。

### [2026-01-21] Wannier90flow 使用方法文档
- **需求**：补充 `wannier_init` 与 `wannier_window` 的使用方法说明。
- **方案**：新增 `wannier90flow/README.md`，包含功能、示例与参数说明。
- **状态**：已完成。

### [2026-01-21] Wannier90 额外输出参数支持
- **需求**：增加 `wvfn_formatted`、`wannier_plot`、`wannier_plot_format` 并在测试配置中说明用途。
- **方案**：
    1. `wannier90flow/wannier_init.py` 写入 `wannier_plot` 与 `wannier_plot_format` 到 `.win`，并支持 `wvfn_formatted` 写入 `pw2wan` 输入。
    2. `tests/inputwannier.toml` 增加参数并写明功能说明。
- **状态**：已完成。

### [2026-01-21] wvfn_formatted 写入 .win
- **需求**：`wvfn_formatted` 应写入 `*.win` 而非 `pw2wan`。
- **方案**：从 `pw2wan` 生成逻辑移除该字段，改为写入 `*.win`，并在 `tests/inputwannier.toml` 的 `[win]` 中说明。
- **状态**：已完成。

### [2026-01-21] 移除误添加的子模块
- **需求**：撤销误引入的 `wannier90flow/external/wannier90` 子模块配置。
- **方案**：删除 `.gitmodules` 并移除子模块记录，清理 README 引用。
- **状态**：已完成。

### [2026-01-19] 功能扩展与规范化：Quantum ESPRESSO 支持及目录命名规范
- **需求**：
    1. 在 `qeflow` 目录下增加 `qe_scf.py`。
    2. 统一任务目录命名规范：`qeflow` 脚本使用 `qe_xxx`（如 `qe_scf`），`vaspflow` 脚本使用 `vasp_xxx`（如 `vasp_scf`）。
- **方案**：
    1. 在 `qeflow/scf/` 目录下创建 `qe_scf.py`。
    2. 更新 `vaspflow` 下所有脚本（scf, eband, eledos, cohp, fermisurface）的 `work_dir` 为 `vasp_xxx`。
    3. 更新 `qeflow/scf/qe_scf.py` 的 `work_dir` 为 `qe_scf`。
    4. 更新 `vaspflow/input.toml` 中的 `chgcar_path` 和 `wavecar_path` 引用路径。
    5. 在 `pyproject.toml` 中注册 `qe_scf` 脚本。
- **状态**：已完成。
