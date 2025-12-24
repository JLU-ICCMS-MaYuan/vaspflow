# Repository Guidelines

## 项目结构与模块
- `cli.py`：主 CLI 入口，可 `python cli.py ...` 或安装后使用 `vasp`；重写多模块串写为 combo，统一调度入口。
- `pipelines/`：核心管线实现（relax、电子性质、phonon、md、批量调度），`utils.py` 负责结构校验与通用工具，结果写入工作目录下的压力子目录。
- `workflows/`：面向用户的工作流封装，编排单任务、批量与后处理，便于脚本化调用。
- `io/`：POSCAR/OUTCAR/vasprun.xml 等读写工具；`config/`：队列模板 `job_templates.toml`（仅作为示例，运行时不读取）；`scheduler/`：调度辅助。
- `analysis/`：绘图辅助；`utils/`：作业工具；`examples/`：放置示例输入；`VASP_error_collecting.md`：故障记录，复现场景前先查阅。

## 开发、构建与运行
- Python 3.12+，`python -m venv .venv && source .venv/bin/activate`，开发时 `export PYTHONPATH=$(pwd):$PYTHONPATH`，无额外构建，可用 `pip install -e .` 做可编辑安装。
- CLI 仅读取 TOML：默认查找工作目录 `job_templates.local.toml`，或通过 `--config` 指定；仓库内 `config/job_templates.toml` 仅作模板且不可直接使用。
- 基本用法：`vasp -i POSCAR --config job_templates.local.toml -p 0 5`；批量：`vasp -i ./structures --config job_templates.local.toml -p 0 5 10`。
- 并行与提交：并发来自 `[tasks].max_workers` 或 `[settings].max_workers`，`submit=true` 时自动提交 `run_*.sh`，否则仅生成输入。
- 调试建议：先在 prepare 模式检查目录结构与 POTCAR，确认无误后再打开 `submit`；日志在各工作子目录 `pipeline.log`。 

## 代码风格与命名
- Python/PEP8，4 空格缩进，使用类型标注；函数与模块用 `snake_case`，类用 `PascalCase`。
- 使用 `logging` 取代 `print`，关键路径写清日志上下文；公共函数加 docstring，错误路径返回明确消息。
- 保持 CLI 参数与 pipelines/workflows 同步，避免硬编码集群特定路径；新增功能优先纯函数，少用全局状态。
- 搜索首选 `rg`/`fd.git`/`sg`，提交前清理 `__pycache__`、临时脚本和大文件。
- 模块内按“常量 → 导入 → 类型定义 → 函数 → 类”组织；新文件默认 ASCII，必要时才引入非 ASCII。

## 测试指南
- 当前无自动化测试；新增功能需自测：准备模式与提交模式各跑一条代表性命令，并保留关键日志。
- 推荐补充 `pytest`，文件命名 `test_*.py`，覆盖解析器、队列调度、并发分支及失败路径，新增逻辑优先 80% 以上覆盖。
- 引入新依赖时附最小可复现命令或验证日志，必要时加假数据避免泄露真实产物。
- 组件或工具函数新增后，给出示例命令、预期输出片段和可能的失败提示，方便回归测试。

## 提交与 PR 规范
- 提交信息格式：`[YYYY-MM-DD HH:MM:SS] <需求类型>：<修改内容>`，完成后立即提交，勿包含 POTCAR/真实账号等敏感数据。
- 变更说明需包含背景、主要修改、影响面、验证方式（命令与关键输出）；队列/资源配置请写明假名。
- 分支/PR 标题简短，包含描述、关联 issue、必要截图或日志；若修改子目录，记得同步该目录下的 AGENTS 记录并注明验证方式。
- 代码评审关注：默认参数、路径拼接、安全边界（文件覆盖、队列命令注入）、并发锁与异常处理；发现风险请在 PR 中明确。

## 配置与安全提示
- 配置来源唯一：`--config` 指定的 TOML 或当前目录 `job_templates.local.toml`；不再读取环境变量、用户级配置或仓库模板。
- `[templates]` 仅使用首个定义的队列头；如果既未定义模板又设置了 job_system，将报错，请至少提供一个 `[templates.<queue>]`。
- `[defaults].mpi_cmd/mpi_procs` 支持数字或完整命令（如 `"mpirun -np 48"` / `"srun"`），可被 `--mpi-procs` 覆盖；缺省则回退 `mpirun -np 8`。
- `[potcar]` 必填：元素 -> POTCAR 绝对路径；缺失元素或路径不存在会直接报错，不再支持 potcar_dir/potcar_lib 搜索。
- `[phonon].structure` 新增开关：primitive（默认）/conventional/relaxed，缺少对应结构会报错。
- 压强目录命名为压力数值（不再带 `_GPa` 后缀），便于与上游目录对齐。
- 调整队列模板时勿提交真实账号路径；批量任务先用少量结构自检，关注 `pipeline.log` 与 `pipeline_report.txt`。
- 产出目录包含 `pipeline_checkpoint.json`/`pipeline_report.txt`/`finished` 标记，可用于排查；异常时先查 `VASP_error_collecting.md`，再最小化复现场景。 

## 协作与记录
- 开工前先做需求优化：澄清 do / don't、输入输出、目标队列与资源；确认后再执行，沟通全程使用中文。
- 每次改动代码或配置，同时更新对应子目录的 AGENTS 记录，描述需求、修改方向、依赖冲突与验证方式；无子目录则在根 AGENTS 补充。
- 文档与 README 按层级清晰排版，示例命令使用反引号或代码块；提交前再次 `git status`，避免遗漏或误提交二进制产物。 

## 需求与修改记录
- 2025-12-15：新增 `pyproject.toml` 与 `setup.cfg` 支持 `pip install -e .`，映射包名 `vasp` 并打包默认模板；`parser_vasp` 在缺少 structuregenerator 时给出可选依赖提示。
- 2025-12-15：添加 `.gitignore` 并清理已跟踪的 `__pycache__/*.pyc` 等构建/缓存产物，避免再次提交；本地工作目录可保留 `test/` 示例结构自测。
- 2025-12-15：新增 `vasp_cli_entry.py` 与入口脚本映射 `vasp/vaspflow -> vasp_cli_entry:main`，在运行时强制优先加载当前项目的 `vasp` 包，避免与其他同名包冲突导致 CLI 无法导入 pipelines。
- 2025-12-15：完善 `config_example.json` 使用说明（已废弃，现已移除并合并到 TOML 模板）。
- 2025-12-16：配置与赝势来源统一为 TOML（仅读取 `job_templates.local.toml` 或 `--config` 指定文件），移除环境变量/用户级兼容；POTCAR 仅按 `[potcar]` 映射拼接；CLI 仍只保留 `-i`/`-p`/`--config`，README 更新为新用法。
- 2025-12-16：队列模板仅采纳首个 `[templates.<queue>]`，缺失模板即报错；声子计算新增结构选择开关 `[phonon].structure`（默认原胞）。 
- 2025-12-23：新增 skill `skills/public/vaspflow-test-runner`（含 run_vaspflow.sh、collect_logs.py、report-template），用于仓库内 VASP 跑测与日志汇总；已打包 `skills/vaspflow-test-runner.skill`，未实际运行 VASP（需按测试输入/配置手动验证）。
- 2025-12-24：结构优化改为四阶段 INCAR+串行脚本，`submit=false` 时仅准备输入不提交，prepare 模式下缺失 CHGCAR/POTCAR 不再阻断后续输入准备（需实际运行后补齐）。
- 2025-12-24：RelaxPipeline 默认断点文件更名为 `relax_checkpoint.json`，兼容读取旧的 `pipeline_checkpoint.json`。
- 2025-12-24：新增 `pipelines/defaults.py` 统一管理最低优先级默认值（提交/并发、kspacing/encut、MD/phonon、队列脚本头、四阶段 relax 等），CLI 与各 Pipeline 均改为引用该集中默认。
