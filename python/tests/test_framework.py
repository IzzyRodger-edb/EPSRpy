"""Crystalline systems from a PDB file: Cu-BTC with the setup API (plan Phase 7a, C20-C23)."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr import Framework, Project, System, tools
from epsr._elements import ELEMENT_MASSES

from epsr_cases import (
    HAS_CORRECTLY_ROUNDED_SINHF,
    OLD_SINHF_REASON,
    PYTHON_DIR,
    compare_dirs,
    extract_archive,
)
from epsr_definitions import (
    CU_BTC_CASE,
    CU_BTC_SETTINGS,
    CU_BTC_TYPES,
    cu_btc_samples,
    cu_btc_system,
)
from epsr_workflow import (
    CU_BTC,
    load_steps,
    prepare_case_folder,
    run_steps_with_wrapper,
)

ARCHIVE = CU_BTC_CASE / "fortran_reference.tar.gz"
EXAMPLE = PYTHON_DIR / "examples" / "08_cu_btc_from_scratch" / "cu_btc_from_scratch.py"
#: Prompt files of the validated steps; the setup API passes the answers itself.
STEP_INPUTS = {"readpdb.stdin", "changeato.stdin"}
#: The steps' EPSR runs (numbers in steps.txt) after which the archive holds the folder.
EPSR_STEPS = (4, 6, 9, 10)


def _needs_case() -> None:
    if not CU_BTC.available or not ARCHIVE.exists():
        pytest.skip("validation/cu-btc, its data or its Fortran archive not available")


def _steps_after(count: int, destination: Path) -> Path:
    """The validated steps' folder after the first ``count`` steps."""
    folder = prepare_case_folder(CU_BTC, destination)
    for _ in run_steps_with_wrapper(
        load_steps(CU_BTC)[:count], folder, CU_BTC.case_folder
    ):
        pass
    for name in STEP_INPUTS:
        (folder / name).unlink()
    return folder


def _assert_identical(reference: Path, candidate: Path) -> None:
    comparison = compare_dirs(reference, candidate)
    assert comparison.all_identical, comparison.report()
    _assert_no_extra_files(reference, candidate)


def _assert_no_extra_files(reference: Path, candidate: Path) -> None:
    extra = {p.name for p in candidate.iterdir()} - {
        p.name for p in reference.iterdir()
    }
    assert extra <= {"epsr_project.json"}, extra


def _assert_matches_fortran(snapshot: Path, candidate: Path) -> None:
    """Every file agrees with the Fortran's, as in ``test_cu_btc.py``.

    Identical, or equal to the last digits, with the tolerance policy of
    ``python/docs/testing.md``. Only meaningful with a correctly rounded ``sinhf``
    (``epsr_cases.HAS_CORRECTLY_ROUNDED_SINHF``), with which the archive was made.
    Byte-identity is required against the validated steps run by the wrapper
    instead (see the test).
    """
    comparison = compare_dirs(snapshot, candidate)
    assert comparison.passed, comparison.report()
    _assert_no_extra_files(snapshot, candidate)


def _steps_snapshots(destination: Path) -> dict[int, Path]:
    """The validated steps through the wrapper: a copy of the folder after each EPSR run."""
    folder = prepare_case_folder(CU_BTC, destination / "work")
    snapshots = {}
    steps = load_steps(CU_BTC)
    for number, result in enumerate(
        run_steps_with_wrapper(steps, folder, CU_BTC.case_folder), start=1
    ):
        if result.step.is_epsr_run:
            snapshot = shutil.copytree(folder, destination / f"after_step_{number:02d}")
            for name in STEP_INPUTS:
                (snapshot / name).unlink()
            snapshots[number] = snapshot
    return snapshots


# --------------------------------------------------------------------- C20: Framework
def test_framework_from_pdb() -> None:
    _needs_case()
    framework = Framework.from_pdb(CU_BTC_CASE / "CuBTC.pdb", CU_BTC_TYPES)
    assert framework.n_atoms == 624 and framework.box_size == 26.3336
    assert sorted(set(framework.residue_sizes)) == [2, 18]
    assert (
        framework.residue_sizes.count(2) == 24
        and framework.residue_sizes.count(18) == 32
    )
    assert str(framework).startswith(
        "CuBTC.pdb: cubic cell 26.3336 Å, 624 atoms, 56 residues"
    )
    system = cu_btc_system()
    assert system.density == pytest.approx(624 / 26.3336**3)
    assert system.n_molecules == 56 and system.box_size == 26.3336
    assert [t.label for t in system.atom_types] == ["CU", "O", "C1", "C2", "C3", "H"]
    assert system.pair_parameters()["CU-O"].terms == "LJ + Coulomb"


def test_framework_checks(tmp_path: Path) -> None:
    _needs_case()
    pdb = CU_BTC_CASE / "CuBTC.pdb"
    with pytest.raises(ValueError, match="order of first appearance"):
        Framework.from_pdb(pdb, list(reversed(CU_BTC_TYPES)))
    wrong_element = [*CU_BTC_TYPES[:5], epsr.AtomType("H", "He", 0, 0, 0)]
    with pytest.raises(ValueError, match="readpdb reads the element"):
        Framework.from_pdb(pdb, wrong_element)
    triclinic = tmp_path / "cell.pdb"
    triclinic.write_text(
        pdb.read_text().replace("90.00  90.00  90.00", "90.00  90.00 120.00", 1)
    )
    with pytest.raises(epsr.EpsrUnsupportedError, match="not cubic"):
        Framework.from_pdb(triclinic, CU_BTC_TYPES)
    framework = Framework.from_pdb(pdb, CU_BTC_TYPES)
    with pytest.raises(ValueError, match="not both"):
        System("CuBTC", 300.0, 0.03, structure=framework)
    with pytest.raises(epsr.EpsrUnsupportedError, match="framework"):
        System("CuBTC", 300.0, structure=framework).add(
            epsr.Molecule.atom(CU_BTC_TYPES[0]), 1
        )
    with pytest.raises(KeyError, match="Monte Carlo"):
        System("CuBTC", 300.0, structure=framework).monte_carlo(stepx=1.0)


# --------------------------------------------------------------------- C20, C21: tools
def test_readpdb_and_changeato_give_the_validated_box(tmp_path: Path) -> None:
    _needs_case()
    folder = tmp_path / "tools"
    folder.mkdir()
    shutil.copy2(CU_BTC_CASE / "CuBTC.pdb", folder)
    tools.readpdb(folder, "CuBTC.pdb", 0.68)
    after_readpdb = _steps_after(1, tmp_path / "readpdb")
    assert (folder / "CuBTC.ato").read_bytes() == (
        after_readpdb / "CuBTC.ato"
    ).read_bytes()
    assert (folder / "vanderWaalsRadii.txt").read_bytes() == (
        CU_BTC_CASE / "vanderWaalsRadii.txt"
    ).read_bytes()
    system = cu_btc_system()
    commands = tools.changeato_commands(
        temperature=system.temperature,
        parameters=system.monte_carlo_parameters,
        atom_types=[t.to_dict() for t in system.atom_types],
        default_masses=ELEMENT_MASSES,
    )
    assert commands == (CU_BTC_CASE / "changeato.stdin").read_text().splitlines()
    tools.changeato(folder, "CuBTC", commands)
    after_changeato = _steps_after(2, tmp_path / "changeato")
    assert (folder / "CuBTC.ato").read_bytes() == (
        after_changeato / "CuBTC.ato"
    ).read_bytes()


def test_create_cu_btc_equals_the_validated_setup(tmp_path: Path) -> None:
    _needs_case()
    project = Project.create(
        tmp_path / "CuBTC", cu_btc_system(), cu_btc_samples(), CU_BTC_SETTINGS
    )
    reference = _steps_after(3, tmp_path / "steps")
    (reference / "CuBTC.EPSR.inp").unlink()  # the steps copy a hand-written one
    comparison = compare_dirs(reference, project.folder)
    assert comparison.all_identical, comparison.report()
    rebuilt = project.system
    assert rebuilt.structure is not None
    assert rebuilt.structure.path == project.folder / "CuBTC.pdb"
    assert rebuilt.monte_carlo_parameters == cu_btc_system().monte_carlo_parameters


# --------------------------------------------------------------------- C22: readers
def test_bragg_readers(tmp_path: Path) -> None:
    _needs_case()
    snapshot = extract_archive(ARCHIVE, tmp_path / "fortran") / "after_step_04"
    bragg = epsr.io.read_bragg(snapshot / "CuBTC.EPSR.bragg")
    assert (bragg.n_reflections, bragg.n_unique) == (80726, 940)
    assert bragg.multiplicity.sum() == 80726 and 0.0 in bragg.q
    intensities = epsr.io.read_inten(snapshot / "CuBTC.EPSR.inten")
    assert intensities.n_pairs == 21 and intensities.bragg.shape == (21, 101)
    assert intensities.q[0] == 0.0 and intensities.q[-1] == pytest.approx(5.0)
    assert np.abs(intensities.convolved).max() > 0


# --------------------------------------------------------------------- C23: the schedule
@pytest.mark.slow
@pytest.mark.experimental
@pytest.mark.xdist_group("cu-btc")
def test_cu_btc_equals_the_fortran_after_every_run(tmp_path: Path) -> None:
    """create + equilibrate(1 + 99) + refine(ereq 10, 2 x 100). After each of the four
    EPSR runs every file is byte-identical to the validated steps run by the wrapper,
    and agrees with the original Fortran's (archive) as in ``test_cu_btc.py``, when
    the platform's ``sinhf`` is correctly rounded like the archive's."""
    _needs_case()
    archive = extract_archive(ARCHIVE, tmp_path / "fortran")
    snapshots = {}
    for number in EPSR_STEPS:
        snapshot = archive / f"after_step_{number:02d}"
        for name in STEP_INPUTS:
            (snapshot / name).unlink(missing_ok=True)
        snapshots[number] = snapshot
    project = Project.create(
        tmp_path / "CuBTC", cu_btc_system(), cu_btc_samples(), CU_BTC_SETTINGS
    )
    wrapper = _steps_snapshots(tmp_path / "steps")
    compared = []

    def check(step: int) -> None:
        _assert_identical(wrapper[step], project.folder)
        if HAS_CORRECTLY_ROUNDED_SINHF:
            _assert_matches_fortran(snapshots[step], project.folder)
        compared.append(step)

    equilibration = project.equilibrate(
        iterations=100, callback=lambda run, simulation: check((4, 6)[run - 1])
    )
    history = project.refine(
        ereq=10.0,
        runs=2,
        iterations=100,
        callback=lambda run, simulation: check((9, 10)[run - 1]),
    )
    assert compared == list(EPSR_STEPS)
    if not HAS_CORRECTLY_ROUNDED_SINHF:
        warnings.warn(f"Cu-BTC archive not compared: {OLD_SINHF_REASON}", stacklevel=1)
    r_factors = [*equilibration.r_factor, *history.r_factor]
    assert r_factors == pytest.approx([3.08e-2, 2.29e-2, 1.69e-2, 1.64e-2], rel=0.02)
    result = project.results()
    bragg = epsr.io.read_bragg(project.output("bragg"))
    assert (bragg.n_reflections, bragg.n_unique) == (80726, 940)
    copper_copper = epsr.analysis.first_peak(
        result.r, result.gr_intra("CU", "CU"), (2.0, 3.2)
    )
    assert copper_copper == pytest.approx(2.64, abs=0.05)
    assert result.potential.pair_params()["CU-O"].terms == "LJ + Coulomb"


@pytest.mark.slow
def test_example_08_quick(tmp_path: Path) -> None:
    _needs_case()
    environment = {
        **os.environ,
        "EPSR_EXAMPLES_WORK": str(tmp_path),
        "MPLBACKEND": "Agg",
    }
    completed = subprocess.run(
        [sys.executable, str(EXAMPLE), "--quick"],
        cwd=PYTHON_DIR,
        env=environment,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "80726 Bragg reflections at 940 unique Q values" in completed.stdout
