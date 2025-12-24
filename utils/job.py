import importlib.util
import logging
import os
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

try:  # Python 3.11+
    import tomllib  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - 兼容旧版
    try:
        import tomli as tomllib  # type: ignore
    except ModuleNotFoundError:
        tomllib = None

logger = logging.getLogger(__name__)

# 兜底默认
from vasp.pipelines import defaults as dft


@dataclass
class JobConfig:
    """封装作业提交相关的用户配置。"""

    vasp_std: str
    vasp_gam: str
    bashtitle: Optional[str]
    slurmtitle: Optional[str]
    pbstitle: Optional[str]
    lsftitle: Optional[str]
    default_mpi_procs: Union[int, str] = 8
    default_queue: str = "bash"


def _load_templates_from_toml(toml_path: Path) -> tuple[dict, dict]:
    """从 TOML 配置读取 defaults 与 queue header 模板。"""
    if not tomllib or not toml_path.exists():
        return {}, {}

    with toml_path.open("rb") as f:
        data = tomllib.load(f)
    defaults = data.get("defaults", {}) or {}
    templates_raw = data.get("templates", {}) or {}
    templates = {k.lower(): (v.get("header") or "").strip() for k, v in templates_raw.items()}
    return defaults, templates


def load_job_config(toml_path: Optional[Path] = None) -> JobConfig:
    """
    读取 VASP 运行配置，仅保留：
    1) 显式传入的 --config
    2) 当前工作目录 job_templates.local.toml
    其他历史兼容路径全部移除。
    """
    if not tomllib:
        raise RuntimeError("当前环境缺少 tomllib/tomli，请使用 Python 3.11+ 或安装 tomli")

    cfg: dict[str, Union[str, int]] = {}
    target = Path(toml_path) if toml_path else (Path.cwd() / "job_templates.local.toml")
    if not target.exists():
        raise ValueError(f"未找到配置文件 {target}，请通过 --config 指定或在工作目录提供 job_templates.local.toml")

    defaults, templates = _load_templates_from_toml(target)

    if defaults.get("vasp_std"):
        cfg["vasp_std"] = defaults["vasp_std"]
    if defaults.get("vasp_gam"):
        cfg["vasp_gam"] = defaults["vasp_gam"]
    if defaults.get("mpi_procs"):
        cfg["default_mpi_procs"] = defaults["mpi_procs"]
    if defaults.get("mpi_cmd"):
        cfg["default_mpi_procs"] = defaults["mpi_cmd"]
    if "default_mpi_procs" not in cfg:
        cfg["default_mpi_procs"] = dft.DEFAULT_MPI_PROCS

    # 兜底队列模板：若未提供 [templates] 则使用最低优先级默认
    if not templates:
        templates = {dft.DEFAULT_QUEUE: dft.DEFAULT_QUEUE_HEADERS.get(dft.DEFAULT_QUEUE, "")}

    # 仅读取第一个定义的队列模板，其余忽略
    first_queue, first_header = next(iter(templates.items()))
    allowed_queues = {"bash", "slurm", "pbs", "lsf"}
    if first_queue not in allowed_queues:
        raise ValueError(f"[templates.{first_queue}] 队列类型不支持，仅允许 {allowed_queues}")
    if not first_header:
        first_header = dft.DEFAULT_QUEUE_HEADERS.get(first_queue, "")
    if not first_header:
        raise ValueError(f"[templates.{first_queue}] 缺少 header 内容，且默认兜底为空")

    cfg["default_queue"] = first_queue or dft.DEFAULT_QUEUE
    cfg["bashtitle"] = first_header if first_queue == "bash" else (dft.DEFAULT_QUEUE_HEADERS.get("bash") if first_queue != "bash" else first_header)
    cfg["slurmtitle"] = first_header if first_queue == "slurm" else (dft.DEFAULT_QUEUE_HEADERS.get("slurm") if first_queue != "slurm" else first_header)
    cfg["pbstitle"] = first_header if first_queue == "pbs" else (dft.DEFAULT_QUEUE_HEADERS.get("pbs") if first_queue != "pbs" else first_header)
    cfg["lsftitle"] = first_header if first_queue == "lsf" else (dft.DEFAULT_QUEUE_HEADERS.get("lsf") if first_queue != "lsf" else first_header)

    required = ["vasp_std", "vasp_gam"]
    missing = [k for k in required if not cfg.get(k)]
    if missing:
        raise ValueError(f"job_templates 配置缺少必需字段: {missing}，请在 {target} 补全")

    return JobConfig(
        vasp_std=str(cfg["vasp_std"]),
        vasp_gam=str(cfg["vasp_gam"]),
        bashtitle=str(cfg.get("bashtitle") or "").strip() if cfg.get("bashtitle") else None,
        slurmtitle=str(cfg["slurmtitle"]).strip() if cfg.get("slurmtitle") else None,
        pbstitle=str(cfg["pbstitle"]).strip() if cfg.get("pbstitle") else None,
        lsftitle=str(cfg["lsftitle"]).strip() if cfg.get("lsftitle") else None,
        default_mpi_procs=cfg.get("default_mpi_procs") or dft.DEFAULT_MPI_PROCS,
        default_queue=str(cfg.get("default_queue") or dft.DEFAULT_QUEUE),
    )


def select_job_header(queue_system: str, cfg: JobConfig) -> str:
    """脚本头选择器，仅使用配置中定义的第一个队列模板。"""
    queue = (queue_system or cfg.default_queue or "").lower()
    if queue not in {"bash", "slurm", "pbs", "lsf"}:
        raise ValueError(f"不支持的队列类型: {queue}")
    header_map = {
        "bash": cfg.bashtitle,
        "slurm": cfg.slurmtitle,
        "pbs": cfg.pbstitle,
        "lsf": cfg.lsftitle,
    }
    header = header_map.get(queue)
    if not header:
        raise ValueError(f"未为队列 {queue} 提供脚本头，请在 [templates.{queue}] 中配置或调整 default_queue")
    return header.strip()


def write_job_script(
    work_dir: Path,
    job_name: str,
    queue_system: str,
    cfg: JobConfig,
    use_gamma: bool = False,
    mpi_procs: Optional[Union[int, str]] = None,
) -> Path:
    """根据 queue_system 生成作业脚本。"""
    work_dir = Path(work_dir)
    queue = (queue_system or cfg.default_queue or "").lower()
    script_file = work_dir / f"run_{job_name}.sh"

    header = select_job_header(queue, cfg)
    binary = cfg.vasp_gam if use_gamma else cfg.vasp_std
    # 支持字符串或数字：字符串直接作为启动命令前缀，数字走 mpirun -np
    run_line: str
    if mpi_procs is None:
        mpi_default = cfg.default_mpi_procs
        if isinstance(mpi_default, str):
            run_line = f"{mpi_default} {binary} > vasp.log"
        else:
            mpi = mpi_default or 8
            run_line = f"mpirun -np {mpi} {binary} > vasp.log"
    elif isinstance(mpi_procs, str):
        mpi_str = mpi_procs.strip()
        if mpi_str.isdigit():
            run_line = f"mpirun -np {mpi_str} {binary} > vasp.log"
        else:
            run_line = f"{mpi_str} {binary} > vasp.log"
    else:
        run_line = f"mpirun -np {mpi_procs} {binary} > vasp.log"

    lines = [header, run_line]

    script_file.write_text("\n".join(lines) + "\n")
    script_file.chmod(0o755)
    return script_file


def submit_job(script_path: Path, queue_system: str) -> str:
    """提交作业，不同队列系统自动选择命令；失败时回落到本地 bash 执行。"""
    queue = (queue_system or "bash").lower()
    script_path = Path(script_path).resolve()
    workdir = script_path.parent
    jobid_file = workdir / "job_id.txt"

    # 如果已有作业ID且仍在运行，直接复用，不再重复提交
    if jobid_file.exists():
        try:
            existing_id = jobid_file.read_text().strip()
            if existing_id:
                from vasp.utils.job import is_job_active  # 避免循环导入
                active = is_job_active(existing_id, queue_system)
                if active:
                    logger.info(f"检测到目录已有运行中的任务 {existing_id}，跳过重复提交")
                    return existing_id
        except Exception:
            pass

    attempts = 0
    while True:
        attempts += 1
        try:
            if queue == "bash":
                subprocess.Popen(["bash", script_path.name], cwd=workdir)
                return "bash"
            if queue == "slurm":
                output = subprocess.check_output(
                    ["sbatch", script_path.name],
                    cwd=workdir,
                    text=True,
                    stderr=subprocess.STDOUT,
                ).strip()
                # 典型输出: "Submitted batch job 123456"
                parts = output.split()
                job_id = parts[-1] if parts else output
                jobid_file.write_text(str(job_id))
                return job_id
            if queue == "pbs":
                output = subprocess.check_output(["qsub", script_path.name], cwd=workdir, text=True).strip()
                jobid_file.write_text(str(output))
                return output
            if queue == "lsf":
                cmd = f"bsub < {shlex.quote(str(script_path))}"
                output = subprocess.check_output(["bash", "-lc", cmd], cwd=workdir, text=True).strip()
                # 典型输出: "Job <123456> is submitted ..."
                import re
                m = re.search(r"<(\d+)>", output)
                job_id = m.group(1) if m else output
                jobid_file.write_text(str(job_id))
                return job_id
        except subprocess.CalledProcessError as exc:
            msg = exc.output if hasattr(exc, "output") else str(exc)
            if queue == "slurm" and "AssocMaxSubmitJobLimit" in msg:
                logger.warning("提交超出作业数限制，等待20秒后重试（第 %d 次）...", attempts)
                time.sleep(20)
                continue
            logger.error("队列提交失败，请检查集群环境或脚本", exc_info=exc)
            raise
        except Exception as exc:  # pragma: no cover - 运行时容错
            logger.error("队列提交失败，请检查集群环境或脚本", exc_info=exc)
            raise

    # 未知队列，回退（理论不会到这里）
    logger.warning("未知队列系统 %s，回退到 bash", queue_system)
    subprocess.Popen(["bash", str(script_path)], cwd=script_path.parent)
    return "bash_unknown"


def is_job_active(job_id: str, queue_system: Optional[str]) -> Optional[bool]:
    """
    查询作业是否仍在队列中。

    返回值：True=仍在运行/排队，False=不在队列，None=无法判断或使用bash。
    """
    queue = (queue_system or "bash").lower()

    if queue == "bash":
        return None

    try:
        if queue == "slurm":
            output = subprocess.check_output(["squeue", "-j", str(job_id)], text=True)
            # squeue 输出含表头 + 任务行
            return len(output.strip().splitlines()) > 1

        if queue == "pbs":
            output = subprocess.check_output(["qstat", str(job_id)], text=True)
            return job_id in output

        if queue == "lsf":
            output = subprocess.check_output(["bjobs", str(job_id)], text=True)
            return job_id in output
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        logger.warning("未找到队列命令，跳过队列状态检查")
        return None
    except Exception as exc:  # pragma: no cover
        logger.warning("队列状态检查异常，跳过", exc_info=exc)
        return None

    return None
