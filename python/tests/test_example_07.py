"""Example 07 (water from scratch, written with the setup API) reproduces the committed run."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from epsr_cases import PYTHON_DIR, compare_dirs
from epsr_definitions import WATER_RUN

SCRIPT = PYTHON_DIR / "examples" / "07_water_from_scratch" / "water_from_scratch.py"
#: Excepted from the comparison, as documented in the example's README: the trend
#: file, which every read of the folder rewrites. (The date line of the .out is masked.)
EXCEPTED = {"waterbox.EPSR.qdr"}
#: Files the setup API keeps in the folder that the validated run does not have.
EXTRA = {"Water.mol", "epsr_project.json"}


def _run_example(work: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    environment = {**os.environ, "EPSR_EXAMPLES_WORK": str(work), "MPLBACKEND": "Agg"}
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--no-3d", *arguments],
        cwd=PYTHON_DIR,
        env=environment,
        capture_output=True,
        text=True,
        timeout=3600,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    return completed


@pytest.mark.slow
@pytest.mark.experimental
@pytest.mark.xdist_group("water")
def test_example_07_reproduces_the_committed_run(tmp_path: Path) -> None:
    if not WATER_RUN.is_dir():
        pytest.skip("the committed example output is not available")
    completed = _run_example(tmp_path)
    run_folder = tmp_path / "output" / "epsr_run"
    committed = shutil.copytree(WATER_RUN, tmp_path / "committed")
    for name in EXCEPTED:
        (committed / name).unlink()
    comparison = compare_dirs(committed, run_folder)
    assert comparison.all_identical, comparison.report()
    produced = {path.name for path in run_folder.iterdir() if path.is_file()}
    expected = {path.name for path in WATER_RUN.iterdir() if path.is_file()}
    assert produced - expected == EXTRA
    assert "OW-OW first peak at 2.79 Å" in completed.stdout
    assert (tmp_path / "output" / "figures" / "fits.png").is_file()
    assert (tmp_path / "output" / "summary.txt").is_file()


@pytest.mark.slow
def test_example_07_quick(tmp_path: Path) -> None:
    completed = _run_example(tmp_path, "--quick", "--seed", "42")
    assert "seed 42" in completed.stdout
    assert "31 EPSR iterations recorded" in completed.stdout
