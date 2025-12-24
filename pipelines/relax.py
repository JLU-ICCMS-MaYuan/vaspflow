import logging
import shutil
from pathlib import Path
from typing import Optional

from vasp.pipelines.base import BasePipeline
from vasp.pipelines.utils import prepare_potcar, ensure_poscar, find_symmetry, check_vasp_completion
from vasp.utils.job import load_job_config, submit_job, JobConfig, select_job_header

logger = logging.getLogger(__name__)


class RelaxPipeline(BasePipeline):
    """仅执行结构优化的 Pipeline。"""

    def __init__(
        self,
        structure_file: Path,
        work_dir: Path,
        kspacing: float = 0.2,
        encut: Optional[float] = None,
        queue_system: Optional[str] = None,
        mpi_procs: Optional[str] = None,
        pressure: float = 0.0,
        potcar_map: Optional[dict[str, str]] = None,
        job_cfg: Optional[JobConfig] = None,
        config_path: Optional[Path] = None,
        **kwargs,
    ):
        super().__init__(
            structure_file,
            work_dir,
            checkpoint_file=Path(work_dir) / "relax_checkpoint.json",
            pressure=pressure,
            **kwargs,
        )

        self.job_cfg = job_cfg or load_job_config(config_path)
        self.kspacing = kspacing
        self.encut = encut
        self.queue_system = queue_system or (self.job_cfg.default_queue if self.job_cfg else None)
        if not self.queue_system:
            raise ValueError("未找到队列配置，请在 job_templates.local.toml 的 [templates] 中至少提供一个队列头")
        self.mpi_procs = mpi_procs
        self.pressure = pressure
        self.potcar_map = potcar_map or {}

        self.relax_dir = self.work_dir / "01_relax"

    def get_steps(self):
        return ["relax"]

    def execute_step(self, step_name: str) -> bool:
        if step_name == "relax":
            return self._run_relax()
        logger.error(f"未知步骤: {step_name}")
        return False

    def _run_relax(self) -> bool:
        logger.info("执行结构优化...")

        self.relax_dir.mkdir(parents=True, exist_ok=True)

        # 若已有收敛结果，直接复用并视为完成
        if check_vasp_completion(self.relax_dir):
            logger.info("检测到已有收敛的 relax 结果，直接复用")
            contcar = self.relax_dir / "CONTCAR"
            if contcar.exists():
                relaxed = self.work_dir / "POSCAR_relaxed"
                shutil.copy(contcar, relaxed)
                self.steps_data["relaxed_structure"] = str(relaxed)
                prim, std, sg = find_symmetry(relaxed, self.work_dir, symprec=1e-3)
                if prim:
                    self.steps_data["primitive_structure"] = str(prim)
                if std:
                    self.steps_data["conventional_structure"] = str(std)
                if sg:
                    self.steps_data["spacegroup"] = sg
            return True

        ensure_poscar(self.structure_file, self.relax_dir / "POSCAR")

        self._write_relax_incar_set(self.relax_dir)
        self._write_kpoints(self.relax_dir / "KPOINTS", self.relax_dir / "POSCAR", self.kspacing)

        if not self.potcar_map:
            logger.error("缺少 [potcar] 配置，无法生成 POTCAR")
            return False
        if not prepare_potcar(
            self.relax_dir / "POSCAR",
            self.potcar_map,
            self.relax_dir / "POTCAR",
        ):
            logger.error("POTCAR准备失败")
            return False

        job_script = self._write_multi_stage_job_script(self.relax_dir, "relax")

        if self.prepare_only:
            logger.info("submit=false，仅准备四阶段结构优化输入和脚本，不提交任务")
            return True

        job_id = self._submit_job(self.relax_dir, job_script)

        if not self._wait_for_job(job_id, self.relax_dir, self.queue_system):
            return False

        if not self._check_convergence(self.relax_dir):
            logger.error("结构优化未收敛")
            return False

        contcar = self.relax_dir / "CONTCAR"
        if contcar.exists():
            relaxed = self.work_dir / "POSCAR_relaxed"
            shutil.copy(contcar, relaxed)
            self.steps_data["relaxed_structure"] = str(relaxed)

            # 对称性分析，生成原胞/标准晶胞（放在压强目录）
            prim, std, sg = find_symmetry(relaxed, self.work_dir, symprec=1e-3)
            if prim:
                self.steps_data["primitive_structure"] = str(prim)
            if std:
                self.steps_data["conventional_structure"] = str(std)
            if sg:
                self.steps_data["spacegroup"] = sg
            self._save_checkpoint()

        logger.info("结构优化完成")
        return True

    def _write_relax_incar_set(self, work_dir: Path) -> None:
        """生成四阶段优化的 INCAR_1..4。"""
        params = [
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
                "encut": self.encut if self.encut else 520,
                "prec": "Accurate",
                "kspacing": self.kspacing,
                "sigma": 0.05,
                "ediff": "1e-6",
                "ediffg": "-0.01",
                "nsw": 500,
                "ibrion": 2,
                "isif": 3,
                "potim": 0.1,
            },
        ]

        for p in params:
            incar_file = work_dir / f"INCAR_{p['suffix']}"
            with open(incar_file, "w") as f:
                f.write("ISTART   = 0\n")
                f.write("ICHARG   = 2\n")
                f.write("ISYM     = 0\n")
                f.write(f"ENCUT    = {p['encut']}\n")
                f.write(f"PREC     = {p['prec']}\n")
                f.write("SYMPREC  = 1e-5\n")
                f.write("NCORE    = 4\n")
                f.write(f"KSPACING = {p['kspacing']}\n")
                f.write("ISMEAR   = 0\n")
                f.write(f"SIGMA    = {p['sigma']}\n")
                f.write("NELM     = 200\n")
                f.write("NELMIN   = 6\n")
                f.write(f"EDIFF    = {p['ediff']}\n")
                f.write(f"EDIFFG   = {p['ediffg']}\n")
                f.write(f"NSW      = {p['nsw']}\n")
                f.write(f"IBRION   = {p['ibrion']}\n")
                f.write(f"ISIF     = {p['isif']}\n")
                f.write(f"POTIM    = {p['potim']}\n")
                f.write("LWAVE    = .FALSE.\n")
                f.write("LCHARG   = .FALSE.\n")
                f.write(f"PSTRESS  = {self.pressure_kbar}\n")

    def _write_kpoints(self, kpoints_file: Path, poscar_file: Path, kspacing: float):
        """写入自动 K 点（由 KSPACING 计算网格数）。"""
        from vasp.pipelines.utils import kspacing_to_mesh
        n1, n2, n3 = kspacing_to_mesh(poscar_file, kspacing)
        with open(kpoints_file, "w") as f:
            f.write("Automatic mesh\n")
            f.write("0\n")
            f.write("Gamma\n")
            f.write(f"{n1} {n2} {n3}\n")
            f.write("0 0 0\n")

    def _write_multi_stage_job_script(self, work_dir: Path, job_name: str) -> str:
        """生成四阶段串行优化的提交脚本。"""
        header = select_job_header(self.queue_system, self.job_cfg)
        run_line = self._build_run_line(self.job_cfg.vasp_std)

        lines = [
            header,
            "set -e",
            "for i in 1 2 3 4; do",
            "  cp INCAR_${i} INCAR",
            "  if [ -f CONTCAR ]; then cp CONTCAR POSCAR; fi",
            "  if command -v killall >/dev/null 2>&1; then killall -9 vasp_std >/dev/null 2>&1 || true; fi",
            f"  {run_line} > vasp.log_${i} 2>&1",
            "  cp CONTCAR CONTCAR_${i} 2>/dev/null || true",
            "  cp OUTCAR  OUTCAR_${i} 2>/dev/null || true",
            "done",
        ]

        script_file = work_dir / f"run_{job_name}.sh"
        script_file.write_text("\n".join(lines) + "\n")
        script_file.chmod(0o755)
        return str(script_file)

    def _build_run_line(self, binary: str) -> str:
        """生成执行 VASP 的运行行，兼容数字或完整 MPI 命令。"""
        if self.mpi_procs is None:
            mpi_default = self.job_cfg.default_mpi_procs
            if isinstance(mpi_default, str):
                return f"{mpi_default} {binary}"
            return f"mpirun -np {mpi_default or 8} {binary}"

        if isinstance(self.mpi_procs, str):
            mpi_str = self.mpi_procs.strip()
            if mpi_str.isdigit():
                return f"mpirun -np {mpi_str} {binary}"
            return f"{mpi_str} {binary}"

        return f"mpirun -np {self.mpi_procs} {binary}"

    def _submit_job(self, work_dir: Path, job_script: str) -> str:
        if self.prepare_only:
            logger.info("submit=false，未提交结构优化任务")
            return "prepare_only"
        return submit_job(Path(job_script), self.queue_system)
