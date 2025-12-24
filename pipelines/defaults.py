"""
集中管理 pipelines 与 CLI 的最低优先级默认值。
所有硬编码默认应引用此模块，便于统一调整与审计。
"""

from __future__ import annotations

from typing import Dict, List, Tuple

# 全局/控制
DEFAULT_PRESSURES: List[float] = [0.0]
DEFAULT_SUBMIT: bool = False
DEFAULT_MAX_WORKERS: int = 2
DEFAULT_STRUCTURE_EXT: str = "vasp"
DEFAULT_DOS_TYPE: str = "element"
DEFAULT_MODULES: List[str] = ["scf", "dos"]
DEFAULT_KSPACING: float = 0.2
DEFAULT_ENCUT: float = 520
DEFAULT_MPI_PROCS = 8  # 支持字符串或数字；此处用数字兜底

# 文件命名
DEFAULT_RELAX_CHECKPOINT_NAME = "relax_checkpoint.json"
DEFAULT_GENERIC_CHECKPOINT_NAME = "pipeline_checkpoint.json"
DEFAULT_REPORT_NAME = "pipeline_report.txt"

# 多阶段结构优化参数
DEFAULT_RELAX_COMMON: Dict[str, str | int | float] = {
    "isym": 0,
    "symprec": "1e-5",
    "ncore": 4,
    "ismear": 0,
    "nelm": 200,
    "nelmin": 6,
}

DEFAULT_RELAX_STAGES: List[Dict[str, str | int | float]] = [
    {
        "suffix": "1",
        "encut": 300,
        "prec": "Normal",
        "kspacing": 0.8,
        "sigma": 0.2,
        "ediff": "1e-3",
        "ediffg": "-0.2",
        "nsw": 100,
        "ibrion": 2,
        "isif": 2,
        "potim": 0.3,
    },
    {
        "suffix": "2",
        "encut": 400,
        "prec": "Normal",
        "kspacing": 0.6,
        "sigma": 0.1,
        "ediff": "1e-4",
        "ediffg": "-0.1",
        "nsw": 200,
        "ibrion": 2,
        "isif": 4,
        "potim": 0.1,
    },
    {
        "suffix": "3",
        "encut": 500,
        "prec": "Accurate",
        "kspacing": 0.4,
        "sigma": 0.05,
        "ediff": "1e-5",
        "ediffg": "-0.05",
        "nsw": 300,
        "ibrion": 2,
        "isif": 3,
        "potim": 0.05,
    },
    {
        "suffix": "4",
        "encut": DEFAULT_ENCUT,
        "prec": "Accurate",
        "kspacing": DEFAULT_KSPACING,
        "sigma": 0.05,
        "ediff": "1e-6",
        "ediffg": "-0.01",
        "nsw": 500,
        "ibrion": 2,
        "isif": 3,
        "potim": 0.1,
    },
]

# 电子性质 INCAR 默认
DEFAULT_RELAX_INCAR = {
    "prec": "Accurate",
    "ediff": "1E-6",
    "ismear": 0,
    "sigma": 0.05,
    "nsw": 200,
    "isif": 3,
    "ediffg": "-0.01",
    "lwav": False,
    "lcharg": True,
}

DEFAULT_SCF_INCAR = {
    "prec": "Accurate",
    "ediff": "1E-8",
    "ismear": 0,
    "sigma": 0.05,
    "lwav": False,
    "lcharg": True,
}

DEFAULT_DOS_INCAR = {
    "prec": "Accurate",
    "icharg": 11,
    "ismear": -5,
    "lorbit": 11,
    "nedos": 2000,
    "lwav": False,
    "lcharg": False,
}

DEFAULT_BAND_INCAR = {
    "prec": "Accurate",
    "icharg": 11,
    "ismear": 0,
    "sigma": 0.05,
    "lorbit": 11,
    "lwav": False,
    "lcharg": False,
}

DEFAULT_ELF_INCAR = {
    "prec": "Accurate",
    "icharg": 11,
    "lelf": True,
    "lwav": False,
    "lcharg": False,
}

DEFAULT_COHP_INCAR = {
    "prec": "Accurate",
    "ismear": -5,
    "lorbit": 11,
    "lwav": False,
    "lcharg": False,
}

DEFAULT_FERMI_INCAR = {
    "prec": "Accurate",
    "icharg": 11,
    "ismear": 0,
    "sigma": 0.05,
    "lorbit": 11,
    "lwav": False,
    "lcharg": False,
}

# MD 默认
DEFAULT_MD = {
    "potim": 1.0,
    "tebeg": 300.0,
    "teend": 300.0,
    "nsw": 200,
    "encut": DEFAULT_ENCUT,
    "ediff": "1E-6",
    "ismear": 0,
    "sigma": 0.05,
}

DEFAULT_MD_RELAX = {
    "prec": "Accurate",
    "encut": DEFAULT_ENCUT,
    "ediff": "1E-6",
    "ismear": 0,
    "sigma": 0.05,
    "ibrion": 2,
    "nsw": 120,
    "isif": 3,
    "ediffg": "-0.02",
}

# 声子默认
DEFAULT_PHONON_SUPERCELL: Tuple[int, int, int] = (2, 2, 2)
DEFAULT_PHONON_METHOD = "disp"
DEFAULT_PHONON_STRUCTURE = "primitive"
DEFAULT_PHONON_RELAX_INCAR = {
    "prec": "Accurate",
    "encut": DEFAULT_ENCUT,
    "ediff": "1E-8",
    "ismear": 0,
    "sigma": 0.01,
    "ibrion": 2,
    "nsw": 200,
    "isif": 3,
    "ediffg": "-1E-5",
}
DEFAULT_PHONON_INCAR = {
    "prec": "Accurate",
    "encut": DEFAULT_ENCUT,
    "ediff": "1E-8",
    "ismear": 0,
    "sigma": 0.01,
    "ibrion": -1,
    "nsw": 0,
}

# 队列与脚本兜底
DEFAULT_QUEUE: str = "bash"
DEFAULT_QUEUE_HEADERS: Dict[str, str] = {
    "bash": "#!/bin/sh\n",
    "slurm": "\n".join(
        [
            "#!/bin/sh",
            "#SBATCH --job-name=vasp",
            "#SBATCH --output=log.out",
            "#SBATCH --error=log.err",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=8",
            "#SBATCH --cpus-per-task=1",
        ]
    ),
    "pbs": "\n".join(
        [
            "#!/bin/sh",
            "#PBS -N vasp",
            "#PBS -j oe",
            "#PBS -l nodes=1:ppn=8",
        ]
    ),
    "lsf": "\n".join(
        [
            "#!/bin/bash",
            "#BSUB -n 8",
            "#BSUB -q normal",
            "#BSUB -J vasp",
            "#BSUB -o log.out",
        ]
    ),
}
