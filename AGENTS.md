# DFTFlow AGENTS.md

## 需求处理记录

### [2026-01-19] 功能扩展：Quantum ESPRESSO 支持
- **需求**：在 `qeflow` 目录下增加 `qe_scf.py`，仿照 `vaspflow/scf` 的功能，用于生成 QE 的 SCF 输入文件。
- **方案**：
    1. 在 `qeflow/scf/` 目录下创建 `qe_scf.py`。
    2. 实现 `QESetup` 类，支持解析 `POSCAR` 结构文件。
    3. 在 `qeflow/input.toml` 中配置 `qe_params`（Namelists）和 `pseudo_map`。
    4. 自动计算 K 点网格（基于 `kmesh` 密度）。
    5. 生成 `qe.in` 和 Slurm 提交脚本 `run_qe.sh`。
    6. 在 `pyproject.toml` 中注册 `qe_scf` 脚本。
- **状态**：已完成。
