"""Run the WAPE cases locally (one after the other) through a ``run.sh`` script."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor


def write_run_script(case_dirs: list[Path], work_dir: str | Path) -> Path:
    """Write ``<work_dir>/run.sh`` that runs ``2D_WAPE`` in every case folder."""
    work_dir = Path(work_dir)
    lines = ["#!/bin/bash"]
    for d in case_dirs:
        rel = os.path.relpath(d, work_dir)
        lines += [f"cd {rel}/", "time ./2D_WAPE >out.out", "cd - >/dev/null"]
    script = work_dir / "run.sh"
    script.write_text("\n".join(lines) + "\n")
    script.chmod(0o755)
    return script

def _run_case(case_dir: Path) -> None:
     exe = next(case_dir.glob("2D_*WAPE"))              # 2D_WAPE or 2D_EWAPE
     with open(case_dir / "out.out", "w") as out, open(case_dir / "err.err", "w") as err:
         subprocess.run([f"./{exe.name}"], cwd=case_dir, stdout=out, stderr=err, check=True)


def run_local(work_dir: str | Path, n_jobs: int = 1) -> None:
    """Run every case of ``work_dir``, ``n_jobs`` at the same time."""
    if n_jobs == 1:
        subprocess.run(["bash", "run.sh"], cwd=work_dir, check=True)
        return
    case_dirs = sorted(p.parent for p in Path(work_dir).glob("t*/*/2D_*WAPE"))
    with ThreadPoolExecutor(n_jobs) as pool:
        list(pool.map(_run_case, case_dirs))




