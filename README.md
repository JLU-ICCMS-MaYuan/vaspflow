# VASP 计算工具使用指南

## 配置与优先级
- 运行入口已简化：命令行仅接受 `-i/--input`（文件或目录）、`-p/--pressure`（可多值，可省略使用配置）以及 `--config`（必填，指定 JSON）。所有模块与参数由配置文件决定。
- 程序与赝势来源（高→低）：① 环境变量 `VASP_STD`/`VASP_GAM`、`POTCAR_DIR`、`VASP_MPI_PROCS`，脚本头可用 `JOB_HEADER_BASH/SLURM/PBS/LSF` 覆盖；② `vasp/config/job_templates.toml` 的 `defaults` 与 `templates`（按 bash/slurm/pbs/lsf 定义头）；③ 兼容旧 `~/.my_scriptrc.py` (`vaspstd_path/vaspgam_path/potcar_dir/*title/default_mpi_procs`)；④ 内置默认。
- 配置优先：模块与参数来自 `--config` 指定的 JSON；命令行仅覆盖输入路径与压强。
- 赝势管理：优先使用当前工作目录下的 `potcar_lib`（支持 `potcar_lib/元素` 或 `potcar_lib/元素/POTCAR`），若缺元素且提供 `potcar_dir` 则仅在 `potcar_dir`（含 `potcar_type` 子目录）中寻找唯一候选并复制到 `potcar_lib`。若找不到或候选不唯一会报错，请手动将所需 POTCAR 放入 `potcar_lib` 后重试。

### config_example.json 使用指引
- 位置：`config_example.json`，可复制为 `my_config.json` 后按实际集群与 POTCAR 路径修改。
- 加载：在命令中通过 `--config my_config.json` 显式指定，未指定不会自动读取。
- 配置字段：`modules`（如 `["relax","scf","dos"]`）、`potcar_dir/potcar_type`、`kspacing/encut`、`job_system/mpi_procs`、`tasks/max_workers`、`submit` 等。
- 优先级：CLI 仅覆盖输入与压强；其余字段以配置为准。
- 常用示例：
  ```bash
  vasp -i POSCAR --config my_config.json -p 0 5
  vasp -i ./structures --config my_config.json
  ```
- 建议每个项目保存一份定制配置（如队列脚本头、POTCAR 路径），敏感信息勿提交仓库。

## 目录规则与执行模式
- 结构输入：`-i` 接受单文件或目录，自动判定批量；`--structure-ext` 过滤扩展名（如 `vasp,cif,res,xsf`，默认 vasp）。
- 压强分层：`-p/--pressure` 支持多值（GPa，写入 INCAR `PSTRESS=pressure*10`），目录为 `<结构前缀>/<pressure>_GPa/01_relax...`，结构前缀为文件去后缀或目录下文件名去后缀。
- 批量与并行：目录输入自动批量；`--tasks N` 控制同时处理的结构数（>1 并行），每个压强/结构生成 `batch_summary.txt`。
- 提交策略：默认仅生成输入与 `run_*.sh`；加 `--submit` 时提交并等待本命令覆盖的全部步骤完成（含批量/多压强、phonon 子位移），结束后在工作目录写 `finished`。
- 状态与重跑：`prepare_only` 模式步骤标记为 `PREPARED`，不会记为完成；重复运行会按产物检查（OUTCAR/CHGCAR 等）决定是否重跑，缺产物或仅 `PREPARED` 一律视为未完成。失败不重试，直接报错退出。

## 步骤依赖与并行编排
- 基础步骤：`relax`；`scf` 依赖 relax；`dos/band/elf/cohp/bader/fermisurface` 依赖 relax+scf；`phonon` 依赖 relax；`md` 依赖 relax。
- 组合执行：`combo` 支持串写模块（如 `relax phonon dos`），程序自动补齐依赖；relax 完成后并行启动 phonon 与电子链路（scf→dos/band/elf/cohp/bader/fermisurface），md 也可并行，统一由 Python 监控。
- 任务提交：bash 直接在目录内启动；队列模式通过 `run_*.sh` 调用 `sbatch/qsub/bsub`，任务号消失后再查 OUTCAR 判定完成。

## 输出、断点与检查
- 每条子链有独立 `pipeline_checkpoint.json` 与 `pipeline_report.txt`；成功结束写 `finished`。已完成且产物存在才会跳过，若缺少关键文件则降级为未完成重跑。
- 批量模式在压强目录下生成 `batch_summary.txt` 汇总成功/失败。

## 常用示例（自动建目录）
- 建议将 `config_example.json` 复制为 `my_config.json`，填写 `modules`（如 `["relax","scf","dos"]`）、`potcar_dir`、`job_system`、`mpi_procs`、`tasks/max_workers`、`submit` 等。
- 单结构：`vasp -i POSCAR --config my_config.json -p 0 5`
- 目录批量 + 多压强：`vasp -i ./structures --config my_config.json -p 0 10`
- 运行前准备赝势：在当前工作目录创建 `potcar_lib`，放入选定赝势（如 `potcar_lib/Si/POTCAR`）；缺失元素时若提供 `potcar_dir` 会复制唯一匹配项到 `potcar_lib`，否则报错提示手动补齐。

这些命令覆盖单结构/批量、多压强场景；其他参数均在配置文件中维护，更新配置即可复用同一命令。
