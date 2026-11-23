"""Run the WAPE cases locally (one after the other) through a ``run.sh`` script."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path


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


def run_local(work_dir: str | Path) -> None:
    """Run ``<work_dir>/run.sh`` and wait for all cases to finish."""
    subprocess.run(["bash", "run.sh"], cwd=work_dir, check=True)
