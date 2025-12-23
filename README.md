# VASP 计算工具使用指南

## 配置来源（唯一）
- 仅接受 TOML：工作目录 `job_templates.local.toml`，或通过 `--config` 指定；仓库内 `config/job_templates.toml` 仅作模板不被读取。
- `[defaults]` 提供 `vasp_std`/`vasp_gam`、`mpi_procs`；`[templates.*]` 为作业脚本头。
- `[settings]` 收录 `modules`（必填）、`kspacing`/`encut`、`job_system`、`max_workers`、`structure_ext`、`dos_type`、`submit`、`pressure` 等通用参数。
- `[phonon]`、`[md]`、`[tasks]` 为专用子段，可覆盖对应参数。
- `[potcar]` 必填：元素符号 -> POTCAR 绝对路径，程序按 POSCAR 元素顺序拼接。

## 命令与执行
- 入口命令：`vasp -i <结构文件或目录> [-p 压强列表] [--config 路径]`。未指定 `--config` 时默认寻找当前目录 `job_templates.local.toml`，缺失即报错。
- 模块顺序完全由 `[settings].modules` 决定（如 `["relax","scf","dos"]`）；`-p` 覆盖配置中的 `pressure`，默认为 `[0]`。
- 目录输入自动批量；并行度来自 `[tasks].max_workers` 或 `[settings].max_workers`（<=1 串行）。
- `submit=true` 时自动提交生成的 `run_*.sh`，否则仅生成输入与脚本。

## 作业脚本与可执行
- VASP 可执行路径与队列头只读自当前配置文件；不再支持环境变量或用户级配置。
- `[templates]` 仅读取首个定义的队列头（如 `templates.slurm`），其余忽略；若未定义任何模板将直接报错。未显式写 `job_system` 时默认使用首个模板的队列；若指定的队列无对应模板也会报错。
- MPI 启动默认用 `[defaults].mpi_cmd/mpi_procs`：可写数字（`8`，生成 `mpirun -np 8`）或完整命令（如 `"mpirun -np 48"`、`"srun"`），`--mpi-procs` 可覆盖。

## POTCAR 处理
- 不再搜索 `potcar_dir` 或 `potcar_lib`。程序仅使用 `[potcar]` 映射，按 POSCAR 中的元素顺序拼接成 POTCAR。
- 缺失元素或路径不存在将直接报错，请先在配置中补全；无需交互选择。

## 目录与产物（概览）
- 单结构：`<stem>/<pressure>/01_relax -> 02_scf -> ...`（压强目录名为压力数值），批量目录下按文件名分组。
- 每个步骤输出 `pipeline_checkpoint.json`、`pipeline_report.txt` 与 `run_*.sh`，完成后写 `finished` 标记。
- 重复运行会基于产物检查决定跳过/重跑；prepare_only 模式仅生成输入不提交。

## 声子结构开关
- `[phonon].structure` 决定声子计算使用的结构：`primitive`（默认，原胞）、`conventional`（标准晶胞）或 `relaxed`（CONTCAR）。若选定结构不存在将报错。
