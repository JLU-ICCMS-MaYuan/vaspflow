# DFTFlow AGENTS.md

## 需求处理记录

### [2026-02-09] vasp_charge_diff 输出格式调整
- **需求**：
    1. `_write_vasp_charge` 以每行 10 个数写出电荷密度。
    2. 数值使用定点小数格式并对齐输出。
- **方案**：
    1. 将每行输出数量由 5 改为 10。
    2. 使用 `{value:12.4f}` 的格式写出电荷密度。
- **状态**：已完成。

### [2026-02-09] vasp_process_locpot 参数合并与多文件电荷密度
- **需求**：
    1. 删除 `--poscar/--locpot` 开关，默认使用当前目录 `POSCAR/LOCPOT`。
    2. 合并 `--parchg/--chgcar` 为 `--chg`，允许多个文件输入。
    3. 自动判断按 CHGCAR/PARCHG 读取与缺失处理，图中 label 使用文件路径名。
- **方案**：
    1. `--chg` 支持多文件并自动探测缺省 `PARCHG/CHGCAR`。
    2. 读取并校验网格一致性后统一按密度输出曲线。
    3. 绘图与数据表头使用文件路径作为标注。
- **状态**：已完成。

### [2026-02-09] vasp_process_locpot 未指定 --chg 时不绘制电荷密度
- **需求**：未指定 `--chg` 时不再自动探测 `PARCHG/CHGCAR`，图中仅绘制 LOCPOT。
- **方案**：仅在提供 `--chg` 时读取电荷密度并绘制曲线，否则跳过。
- **状态**：已完成。

### [2026-02-09] vasp_process_locpot 未指定 --chg 报错修复
- **需求**：未指定 `--chg` 时不应遍历空对象导致报错。
- **方案**：仅在 `--chg` 非空时迭代读取电荷密度文件。
- **状态**：已完成。

### [2026-02-09] vasp_opt 结构优化脚本
- **需求**：
    1. 新增结构优化脚本，支持 rv4/rv1/rvf 三种模式。
    2. 生成对应 INCAR 与 Slurm 脚本，工作目录为 `vasp_opt`。
    3. 注册脚本入口。
- **方案**：
    1. 新增 `vaspflow/opt/vasp_opt.py` 与 `__init__.py`。
    2. 按模板生成 INCAR_1~4 与 `fouropt.sh/oneopt.sh/fopt.sh`。
    3. 在 `pyproject.toml` 中注册 `vasp_opt` 入口。
- **状态**：已完成。

### [2026-02-10] vasp_opt 示例配置补充
- **需求**：补充 `tests/inputvasp.toml` 的 `[vasp_inputpara]` 示例配置。
- **方案**：新增结构优化所需字段示例值。
- **状态**：已完成。

### [2026-02-10] vasp_opt 配置合并至 incar_params
- **需求**：删除 `[vasp_inputpara]`，结构优化所需 INCAR 参数全部并入 `[incar_params]`。
- **方案**：`vasp_opt` 仅读取 `[incar_params]` 生成 rv4 最后一步与 rv1/rvf 的 INCAR。
- **状态**：已完成。

### [2026-02-10] vasp_opt KSPACING 由 kmesh 兜底
- **需求**：`KSPACING` 未配置时使用顶层 `kmesh`。
- **方案**：`incar_params` 缺少 `KSPACING` 时读取 `kmesh` 写入 INCAR。
- **状态**：已完成。

### [2026-02-09] vasp_process_locpot 多点剖面分析
- **需求**：
    1. 支持多点路径抽样 LOCPOT/PARCHG/CHGCAR 的 1D 剖面。
    2. 支持分数/笛卡尔坐标与路径采样步长或点数。
    3. 输出数据文件、曲线图与几何路径图。
    4. 无点输入时保留 LOCPOT min/max 简化模式。
- **方案**：
    1. 复用 POSCAR 头解析与体数据读取逻辑。
    2. 实现路径重采样与三线性插值采样。
    3. 生成 `*.dat`、`*.png` 与 `*_geom.png` 输出。
    4. 无点输入时仅输出 LOCPOT min/max。
- **状态**：已完成。

### [2026-02-07] vaspflow skill（locpot 工作流）
- **需求**：
    1. 设计通用 vaspflow 技能，适配新增 `vasp_locpot` 子功能。
    2. 目录结构与脚本骨架对齐 `vaspflow/scf`。
    3. 明确配置字段、输出文件与命令入口规范。
- **方案**：
    1. 在 `skills/vaspflow/` 创建 `SKILL.md`。
    2. 定义目录模板、脚本方法清单、locpot 关键参数建议。
    3. 说明 input.toml 配置字段与输出规范。
- **状态**：已完成。

### [2026-02-07] vaspflow 概览 PPT
- **需求**：
    1. 基于 `vaspflow` 代码生成面向新入组学生的中文概览 PPT。
    2. 覆盖整体结构、模块职责、配置文件与运行方式。
    3. 输出放在 `vaspflow/` 目录下。
- **方案**：
    1. 读取 `vaspflow` 目录结构与各模块脚本概要。
    2. 按模块与通用流程组织 10+ 张 slides。
    3. 生成 `vaspflow_overview.pptx`。
- **状态**：已完成。

### [2026-02-07] vasp_locpot 脚本与配置
- **需求**：
    1. 新增 `vasp_locpot.py`，支持拷贝 SCF 的 `CHGCAR/WAVECAR`。
    2. 默认写入 `LVHAR = .TRUE.`，网格参数从配置读取。
    3. `input.toml` 增加 `[locpot_params]`，`pyproject.toml` 增加脚本入口。
- **方案**：
    1. 新增 `vaspflow/locpot/vasp_locpot.py` 与 `__init__.py`。
    2. `locpot_params` 覆盖 INCAR 网格参数（无配置则不写入）。
    3. 注册 `vasp_locpot` 脚本入口。
- **状态**：已完成。

### [2026-02-07] vasp_charge 能量范围电荷密度
- **需求**：
    1. 新增能量范围电荷密度脚本 `vasp_charge.py`。
    2. 默认写入 `LORBIT=11`、`LPARD=.TRUE.`。
    3. `EINT` 与网格参数由 `[charge_params]` 提供。
    4. 拷贝 SCF 的 `CHGCAR/WAVECAR`。
- **方案**：
    1. 新增 `vaspflow/charge/vasp_charge.py` 与 `__init__.py`。
    2. `input.toml` 增加 `[charge_params]`（含 `EINT` 与网格参数）。
    3. `pyproject.toml` 注册 `vasp_charge` 入口。
- **状态**：已完成。

### [2026-02-07] vasp_charge EINT 数组写入格式
- **需求**：`EINT = [8, 16.6861]` 读入后写为 `EINT = 8 16.6861`。
- **方案**：在 `vasp_charge.py` 中对 `incar_params.EINT` 做数组转空格字符串的处理，并移除 `charge_params` 读取。
- **状态**：已完成。

### [2026-02-07] vasp_charge 目录名按 EINT
- **需求**：根据 `EINT = [emin, emax]` 自动设置工作目录为 `vasp_charge___emin___emax`。
- **方案**：在 `vasp_charge.py` 初始化阶段读取 `incar_params.EINT`，若为数组则按三位小数格式更新 `work_dir`。
- **状态**：已完成。

### [2026-02-07] vasp_process_locpot 脚本入口
- **需求**：在 `pyproject.toml` 中新增 `vasp_process_locpot` 入口。
- **方案**：将入口指向 `vaspflow.locpot.vasp_process_locpot:main`。
- **状态**：已完成。

### [2026-02-07] vasp_charge_diff 电荷密度差分
- **需求**：
    1. 新增 `vasp_charge_diff`，计算电荷密度差 A - B。
    2. 同时支持 `PARCHG` 与 `CHGCAR`。
    3. 头信息/坐标/网格必须一致，否则报错。
    4. 输出文件名默认 `A_minus_B.vasp`。
- **方案**：
    1. 新增 `vaspflow/charge/vasp_charge_diff.py`。
    2. 读取两份电荷密度并逐点相减。
    3. 在 `pyproject.toml` 注册脚本入口。
- **补充**：兼容未包含元素行的 VASP 头信息，自动识别原子数行与坐标起始行。
- **补充**：写出数值时加空格分隔，避免连写。
- **状态**：已完成。

### [2026-01-30] postw90 多投影自动执行
- **需求**：
    1. `inputw90post.toml` 中出现多个 `dos_project` 时，循环执行多次 `postw90.x`。
    2. 每次输出按 `seedname_dos_1_5.dat` 形式重命名，避免覆盖。
    3. 支持 `dos_project = ["1:5", "6:11"]` 列表写法。
    4. 支持 `dos_project = { Y_d = "1:5", H_p = "6:11" }` 字典写法。
- **方案**：
    1. `wannier_post.py` 读取 `inputw90post.toml`，合并重复 `dos_project` 为列表。
    1. 允许 `dos_project` 直接为列表，按列表顺序逐次执行。
    2. 逐个投影写入临时 `.win`，执行 `postw90.x`，重命名输出文件。
    3. 运行结束后恢复原始 `.win` 内容。
    4. 新增 `tests/inputw90post.toml` 示例并更新文档。
    5. 字典写法使用键名作为输出文件后缀。
- **状态**：已完成。

### [2026-01-30] wannier_post 参数调整
- **需求**：
    1. `wannier_post.py` 不再支持 `-w/-p`，改为 `-i` 指定输入文件。
    2. 使用 `-i` 自动确定工作目录与前缀。
    3. 更新使用文档。
- **方案**：
    1. 解析 `-i` 为 `.win` 文件路径（或前缀），从文件名获取 prefix，从路径获取 work_dir。
    2. `run_postw90` 直接使用推断的目录与前缀执行。
    3. 更新 `wannier90flow/README.md` 示例。
- **状态**：已完成。

### [2026-01-29] Wannier90 postw90 流程拆分
- **需求**：
    1. 将 `postw90.x` 流程从 `wannier_init.py` 迁移到独立脚本。
    2. `wannier_init.py` 移除 `postw90.x` 的执行与脚本生成，避免重复。
    3. 新增 `wannier_post.py` 支持生成脚本与 `--run` 本地执行。
    4. 更新脚本入口与使用文档。
- **方案**：
    1. 新增 `wannier90flow/wannier_post.py`，读取 `[postw90].executable_path` 并生成 `run_postw90.sh`，`--run` 执行 `postw90.x`。
    2. `wannier90flow/wannier_init.py` 移除 `postw90.x` 相关逻辑与提示。
    3. `pyproject.toml` 增加 `wannier_post` 脚本入口。
    4. 更新 `wannier90flow/README.md` 与示例配置说明。
- **状态**：已完成。

### [2026-01-29] Wannier90 DOS 参数写入 .win
- **需求**：
    1. 在 `wannier90flow/wannier_init.py` 增加 `dos`、`dos_kmesh` 与 `dos_project` 参数读取。
    2. `dos`、`dos_kmesh`、`dos_project` 写入 `*.win`，并更新示例配置。
    3. `--run` 模式下补充执行 `postw90.x`。
- **方案**：
    1. 从 `[win]` 读取 `dos/dos_kmesh/dos_project`，写入 `.win` 文件。
    2. `dos_kmesh` 支持空格分隔字符串或数组输入，输出为 `dos_kmesh = 25 25 25` 形式。
    3. `dos_project` 支持字符串或数组输入，输出为 `dos_project = 1:5` 形式。
    4. `--run` 流程追加 `postw90.x`，路径从 `[postw90].executable_path` 读取。
    5. 更新 `tests/inputwannier.toml` 示例。
- **状态**：已完成。

### [2026-01-29] Wannier90 能带对比与 QE eband 输出命名调整
- **需求**：
    1. `wannier90_plot.py` 改为读取 `--qe/--w90` 目录并自动发现能带与高对称信息文件，默认当前目录。
    2. `qe_eband.py` 输出文件名调整：`filband` 为 `{prefix}_band`、`filproj` 为 `{prefix}_band_proj`，高对称点文件名为 `{prefix}__band.labelinfo.dat`。
    3. `wannier90_plot.py` 的 QE 与 Wannier90 文件匹配规则与 `.gnu` 能带格式解析。
    4. Wannier90 能带文件需按空行分段绘制，避免不同能带连线。
    5. `wannier90_plot.py` 删除 `--label` 参数，完全依赖目录自动匹配。
- **方案**：
    1. 新增目录内文件自动匹配逻辑与可选 `--label` 覆盖；标签文件优先从 Wannier90 目录查找，不存在则回退 QE 目录。
    2. `qe_eband.py` 生成输入文件时固定按前缀命名，并在写高对称点文件时使用新文件名。
    3. QE 目录匹配 `*_band.gnu`，Wannier90 目录匹配 `*_band.dat`，并兼容 `&plot`/`.gnu` 两种格式的 QE 能带解析。
    4. 解析 Wannier90 的 `*_band.dat` 时按空行拆分为单条能带并逐条绘制。
    5. 移除 `--label` 分支，标签文件仅通过目录匹配自动发现。
    6. QE/Wannier90 能带文件统一读取 `*_band.dat`，并跳过注释行。
    7. 绘图时 QE 使用虚线且置于上层，W90 置于下层；y 轴范围取 W90 能量最小值与最大值的 1.2 倍。
- **状态**：已完成。

### [2026-01-29] QE 投影能带 kpath 距离与 Wannier90 对齐
- **需求**：
    1. `qe_process_eband.py` 使用 `eband.in` 的 `CELL_PARAMETERS` 计算 kpath 距离，与 Wannier90 度量一致。
    2. 新增 `--cell` 参数（默认当前目录 `eband.in`），弃用旧的分数坐标欧氏距离算法。
    3. 输出 `f"{prefix}_band.labelinfo.dat"` 与两列能带 `f"{prefix}_band.dat"`。
- **方案**：
    1. 读取 `eband.in` 的晶胞参数，构造倒格子度量矩阵并计算 `kpath_dist`。
    2. 生成 `*_band.labelinfo.dat` 仅记录 `eband.in` 高对称点坐标与距离。
    3. 额外输出两列能带数据文件，按能带分段写入，并统一去除 `_band_proj` 前缀用于 `*_band.dat` 与 `*_band.labelinfo.dat` 命名。
- **补充**：读取 `K_POINTS crystal_b` 权重，跳转点（权重=1）不计入距离，以对齐 Wannier90 的 `kpath_dist`。
- **补充**：若 `K_POINTS` 权重与 bands 文件点数不匹配，则改用“步长突变阈值”自动识别跳转点。
- **补充**：KList 由 `K_POINTS crystal_b` 插值生成，能量仅从 `*_band` 读取，并校验生成点数与 nks 一致。
- **状态**：已完成。

### [2026-01-29] QE eband 高对称路径不连续处理
- **需求**：
    1. 解析 `KPATH.in` 识别连续路径段并输出每段路径。
    2. 生成 `*_band.labelinfo.dat` 时，若路径不连续则用上一段终点替换下一段起点，距离连续累积。
- **方案**：
    1. 按两行一段解析 `KPATH.in`，以“上一段终点坐标 == 下一段起点坐标(1e-6)”判定连续。
    2. 生成路径段日志，并在标签输出时替换不连续段起点坐标/标签。
    3. 写 `K_POINTS crystal_b` 时，对段末点设置权重为 1 以形成跳转点，且最后一个高对称点权重固定为 1。
    4. 额外输出 `{prefix}_band.labelinfo.dat`，包含索引列与与 Wannier90 一致的累计距离。
    5. `{prefix}_band.labelinfo.dat` 使用 Wannier90 的 2π 倒格子度量，`*_gun_band.labelinfo.dat` 保持 QE 2π/alat。
    6. 不连续段起点距离不累加（跳转段不计入距离），与 Wannier90 的 labelinfo 对齐。
- **状态**：已完成。

### [2026-01-29] QE 高对称路径标签文件归并
- **需求**：
    1. 将 `qe_process_eband.py` 中的 `*_band.labelinfo.dat` 输出移除。
    2. `qe_eband.py` 继续输出标签文件，但重命名为 `*_gun_band.labelinfo.dat`。
- **方案**：
    1. 删除 `qe_process_eband.py` 的 labelinfo 写入逻辑。
    2. `qe_eband.py` 将标签文件名改为 `_gun_band.labelinfo.dat`。
- **状态**：已完成。

### [2026-01-22] wannier_init 复制输入结构到工作目录
- **需求**：`-i` 指定的结构文件需拷贝到 `work_dir` 下并命名为 `POSCAR`。
- **方案**：在 `wannier90flow/wannier_init.py` 中解析结构后，复制输入文件到 `wannier90/POSCAR`。
- **状态**：已完成。

### [2026-01-22] 修复 vaspkit KPATH.in 解析
- **需求**：`parse_kpath` 需要正确解析 vaspkit 303 生成的 `KPATH.in`（两行一段）。
- **方案**：按“坐标+标签”两行一段收集点，再成对生成路径段。
- **状态**：已完成。

### [2026-01-22] QE 赝势名仅从 pseudo_map 读取
- **需求**：`ATOMIC_SPECIES` 伪势文件名应直接使用 `pseudo_map` 中对应元素的 `pseudo` 字段。
- **方案**：`qeflow/scf/qe_scf.py` 与 `qeflow/nscf/qe_nscf.py` 移除默认回退，缺失配置时报错。
- **状态**：已完成。

### [2026-01-22] 新增 QE eband 输入生成脚本
- **需求**：仿照 `qe_nscf.py` 新增 `qe_eband.py`，高对称路径通过 vaspkit 生成。
- **方案**：新增 `qeflow/eband/qe_eband.py`，调用 vaspkit 生成 `KPATH.in`，解析后写入 `K_POINTS crystal_b`；生成 `eband.in/elebanddata.in/elebandprojdata.in` 与运行脚本，并补充脚本入口。
- **状态**：已完成。

### [2026-01-22] QE 可执行路径统一为 bin 目录
- **需求**：`[qe].executable_path` 只保存到 `qe-7.4.1/bin`，各脚本自动拼接 `pw.x`、`bands.x`、`projwfc.x`。
- **方案**：`qe_scf.py`/`qe_nscf.py` 使用 `executable_path/pw.x`；`qe_eband.py` 使用 `executable_path/pw.x/bands.x/projwfc.x`，并更新示例配置。
- **状态**：已完成。

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
- **方案**：在 `wannier90flow/wannier_window.py` 中使用 matplotlib 绘制 `1..(nbnd_win+8)` 的能带能量区间柱状图，并标出 `dis_froz_min/max` 与 `dis_win_min/max` 虚线，若缺少 matplotlib 则提示跳过。
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

### [2026-01-21] Wannier90 fermi 参数支持
- **需求**：支持 `fermi_energy` 与 `fermi_surface_plot` 写入 `*.win`。
- **方案**：在 `wannier90flow/wannier_init.py` 的 `.win` 输出中增加对应字段。
- **状态**：已完成。

### [2026-01-21] Wannier90 UNK 文件生成提示
- **需求**：在测试配置中标注 `write_unk` 用途，避免 `wannier_plot` 报错缺少 UNK。
- **方案**：在 `tests/inputwannier.toml` 的 `[pw2wannier90]` 中增加 `write_unk = true` 并说明用途。
- **状态**：已完成。

### [2026-01-21] bands_plot 自动生成 KPATH
- **需求**：当 `bands_plot = true` 且缺少 `KPATH.in` 时自动运行 vaspkit 303；`bands_plot = false` 时不执行。
- **方案**：在 `wannier90flow/wannier_init.py` 中检测 `bands_plot` 并调用 vaspkit 生成 `KPATH.in`，失败则报错退出。
- **状态**：已完成。

### [2026-01-23] 新增：QE eledos 输入生成脚本
- **需求**：新增 `qe_eledos.py`，基于 `qe_scf` 的自洽结果生成电子 DOS 所需输入并支持运行。
- **方案**：
    1. 新增 `qeflow/eledos/qe_eledos.py`：生成 `eletdos.in/elepdos.in`，拷贝 `qe_scf` 的 `prefix.save` 数据，支持 `--run` 执行 `dos.x/projwfc.x`。
    2. `pyproject.toml` 注册 `qe_eledos` 入口，`qeflow/input.toml` 补充 `[eledos]` 示例参数。
- **状态**：已完成。

### [2026-01-23] 新增：QE 投影能带数据处理脚本
- **需求**：新增 `qe_process_eband.py`，将 `projwfc` 投影能带文件处理为三列数据（横坐标/能量/轨道贡献），并输出所有轨道及 p/d/f 汇总。
- **方案**：
    1. 新增 `qeflow/eband/qe_process_eband.py`：解析 `*.projwfc_up/down`，按能带顺序输出各轨道与 p/d/f 汇总，生成 `<prefix>_eleband_proj.dat`。
    2. `pyproject.toml` 注册 `qe_process_eband` 入口。
- **状态**：已完成。

### [2026-01-23] 新增：QE 投影 DOS 数据处理脚本
- **需求**：新增 `qe_process_dos.py`，处理 `*.tdos` 与 `pdos_atm` 文件，输出 DOS/IDOS 表格并包含 TDOS/ITDOS。
- **方案**：
    1. 新增 `qeflow/eledos/qe_process_dos.py`：解析 `*.tdos`、`pdos_atm` 与 `projwfc`，生成 `<prefix>_dos_proj*.dat` 与 `<prefix>_idos_proj*.dat`。
    2. 输出包含同元素轨道汇总与未汇总两种版本，TDOS/ITDOS 置于最后一列。
- **状态**：已完成。

### [2026-01-23] 修正：qe_process_eledos 脚本入口
- **需求**：`qe_process_eledos` 入口需可用。
- **方案**：`pyproject.toml` 中将 `qe_process_eledos` 指向 `qeflow.eledos.qe_process_dos:main`。
- **状态**：已完成。

### [2026-01-22] vaspkit 303 静默运行
- **需求**：执行 vaspkit 303 时不打印输出到屏幕。
- **方案**：调用 vaspkit 时重定向 stdout/stderr 到 DEVNULL。
- **状态**：已完成。

### [2026-01-22] bands_plot 写入位置统一
- **需求**：`bands_plot` 仅由统一配置块写入，不在 `kpoint_path` 段重复写入。
- **方案**：移除 `kpoint_path` 写入中的 `bands_plot = .true.`，保留统一逻辑。
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
