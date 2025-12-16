# Repository Guidelines

## 项目结构与模块
- `cli.py`：主 CLI 入口，可 `python cli.py ...` 或安装后使用 `vasp`；重写多模块串写为 combo，统一调度入口。
- `pipelines/`：核心管线实现（relax、电子性质、phonon、md、批量调度），`utils.py` 负责结构校验与通用工具，结果写入工作目录下的压力子目录。
- `workflows/`：面向用户的工作流封装，编排单任务、批量与后处理，便于脚本化调用。
- `io/`：POSCAR/OUTCAR/vasprun.xml 等读写工具；`config/`：队列模板 `job_templates.toml` 与示例 `config_example.json`；`scheduler/`：调度辅助。
- `analysis/`：绘图辅助；`utils/`：作业工具；`examples/`：放置示例输入；`VASP_error_collecting.md`：故障记录，复现场景前先查阅。

## 开发、构建与运行
- Python 3.12+，`python -m venv .venv && source .venv/bin/activate`，开发时 `export PYTHONPATH=$(pwd):$PYTHONPATH`，无额外构建。
- 查看帮助：`python cli.py --help`。仅准备输入：`python cli.py relax -i POSCAR --kspacing 0.2 --encut 520`（不提交队列），便于本地快速验证。
- Slurm 提交示例：`python cli.py relax -i POSCAR -j slurm --mpi-procs "srun -n 32" --submit`；批量并行：`python cli.py relax -i ./structures --tasks 4 --submit`。
- 读取 JSON 配置：`python cli.py dos -i POSCAR --json config_example.json --submit`；优先使用 CLI 参数覆盖配置。
- 本地调试可先运行 prepare-only，再在同一目录用 `--submit` 复用产物；如需可编辑安装，使用 `pip install -e .` 方式。

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
- 优先级：CLI 参数最高，其次 `config/job_templates.toml`，然后环境变量 `VASP_STD/VASP_GAM/POTCAR_DIR/VASP_MPI_PROCS`，最后内置默认；示例 JSON 需显式通过 `--json` 传入。
- 提交前核查 `potcar_lib/` 与 `--potcar-type`，批量任务先用少量结构试跑；根据集群调整模板头部，避免提交真实队列账号或路径。
- 产出目录包含 `pipeline_checkpoint.json`/`pipeline_report.txt`/`finished` 标记，可用于排查；异常时先查 `VASP_error_collecting.md`，再最小化复现场景。 

## 协作与记录
- 开工前先做需求优化：澄清 do / don't、输入输出、目标队列与资源；确认后再执行，沟通全程使用中文。
- 每次改动代码或配置，同时更新对应子目录的 AGENTS 记录，描述需求、修改方向、依赖冲突与验证方式；无子目录则在根 AGENTS 补充。
- 文档与 README 按层级清晰排版，示例命令使用反引号或代码块；提交前再次 `git status`，避免遗漏或误提交二进制产物。 

## 需求与修改记录
- 2025-12-15：新增 `pyproject.toml` 与 `setup.cfg` 支持 `pip install -e .`，映射包名 `vasp` 并打包默认模板；`parser_vasp` 在缺少 structuregenerator 时给出可选依赖提示。
- 2025-12-15：添加 `.gitignore` 并清理已跟踪的 `__pycache__/*.pyc` 等构建/缓存产物，避免再次提交；本地工作目录可保留 `test/` 示例结构自测。
- 2025-12-15：新增 `vasp_cli_entry.py` 与入口脚本映射 `vasp/vaspflow -> vasp_cli_entry:main`，在运行时强制优先加载当前项目的 `vasp` 包，避免与其他同名包冲突导致 CLI 无法导入 pipelines。
- 2025-12-15：完善 `config_example.json` 使用说明，指导复制为自定义配置并通过 `--json` 显式加载，明确 CLI 覆盖优先级与示例命令。
- 2025-12-16：简化 CLI，只保留 `-i/--input`、`-p/--pressure` 与 `--config`，模块与其他参数全部由配置文件提供；更新 `config_example.json` 添加 `modules`/`submit`，重写 README 用法说明。
