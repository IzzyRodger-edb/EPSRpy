"""The scripts in ``python/examples`` run end to end.

Each example is a folder ``NN_name/`` with one script. Example 07 has its own tests
(``test_example_07.py``: the full run against the committed output, and a quick run).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from epsr_cases import PYTHON_DIR

EXAMPLE_FOLDERS = sorted(
    folder
    for folder in (PYTHON_DIR / "examples").glob("[0-9][0-9]_*")
    if folder.is_dir() and folder.name != "07_water_from_scratch"
)
#: Command-line arguments that keep the longer examples short under test:
#: 10 iterations per EPSR run, no 3D image (pyvista needs an OpenGL context).
SHORT_RUN_ARGUMENTS = {
    "02_refinement_loop": ["--runs", "3"],
    "05_water_existing_epsr_files": ["--quick"],
    "06_sio2_existing_epsr_files": ["--quick"],
    "08_cu_btc_from_scratch": ["--quick"],
    "09_sio2_from_scratch": ["--quick", "--no-3d"],
    "10_benzene_from_scratch": ["--quick", "--no-3d"],
}
#: Examples that take over a minute even in their short mode.
SLOW = {"08_cu_btc_from_scratch", "09_sio2_from_scratch", "10_benzene_from_scratch"}


def _script(folder: Path) -> Path:
    (script,) = folder.glob("*.py")  # one script per example folder
    return script


@pytest.mark.parametrize(
    "folder",
    [
        pytest.param(folder, marks=[pytest.mark.slow] if folder.name in SLOW else [])
        for folder in EXAMPLE_FOLDERS
    ],
    ids=[folder.name for folder in EXAMPLE_FOLDERS],
)
def test_example_runs(folder: Path, tmp_path: Path) -> None:
    environment = {
        **os.environ,
        "EPSR_EXAMPLES_WORK": str(tmp_path),
        "MPLBACKEND": "Agg",
    }
    command = [
        sys.executable,
        str(_script(folder)),
        *SHORT_RUN_ARGUMENTS.get(folder.name, []),
    ]
    completed = subprocess.run(
        command,
        cwd=PYTHON_DIR,
        env=environment,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert completed.stdout.strip()


def test_examples_found() -> None:
    assert [folder.name for folder in EXAMPLE_FOLDERS] == [
        "01_single_run",
        "02_refinement_loop",
        "03_perturb_potential",
        "04_tools",
        "05_water_existing_epsr_files",
        "06_sio2_existing_epsr_files",
        "08_cu_btc_from_scratch",
        "09_sio2_from_scratch",
        "10_benzene_from_scratch",
    ]
    for folder in EXAMPLE_FOLDERS:
        assert (folder / "README.md").is_file(), folder.name
        _script(folder)
