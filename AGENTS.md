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
    1. 新增 `wannier90flow/prepare_wannier90.py`：读取 POSCAR 与 TOML/JSON 配置，生成包含 `unit_cell_cart`、`atoms_*`、`projections`、`mp_grid`、显式 `kpoints` 及可选 `kpoint_path` 的 `.win` 文件。
    2. 添加示例配置 (`wannier90flow/sample_wannier90.toml`，另有 `wannier90flow/input.toml` 可直接使用)，复用 Ce1Sc2H24 参数，演示投影、能窗与网格设置；`kpoint_path` 通过 vaspkit 303 生成的 `KPATH.in` 自动读取。
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
