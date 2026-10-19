"""``Project.create`` and the refinement vocabulary (plan Phase 5, C8-C17)."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr import NeutronSample, Project, ProjectStateError

from epsr_cases import CASE_FILES_TO_SKIP, TEST_PROBLEMS, compare_dirs, copy_files
from epsr_definitions import (
    SILICA_SETTINGS,
    WATER_RUN,
    WATER_SETTINGS,
    silica_samples,
    silica_system,
    water_samples,
    water_system,
)
from epsr_workflow import (
    SILICA,
    WATER,
    load_steps,
    prepare_case_folder,
    run_steps_with_wrapper,
    shorten_schedule,
)

#: Files of the API's folder that the validated steps do not have.
API_EXTRA_FILES = {"epsr_project.json", "Water.mol"}


def _needs_water() -> None:
    if not WATER.available:
        pytest.skip("validation/water or its data not available")


def _hashes(folder: Path) -> dict[str, str]:
    return {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in folder.iterdir()
        if path.is_file()
    }


def _steps_folder(case, destination: Path, steps) -> Path:  # type: ignore[no-untyped-def]
    folder = prepare_case_folder(case, destination)
    for _ in run_steps_with_wrapper(steps, folder, case.case_folder):
        pass
    for stdin_file in folder.glob(
        "*.stdin"
    ):  # prompt files of the steps, not EPSR files
        stdin_file.unlink()
    return folder


def _setup_steps(case) -> list:  # type: ignore[no-untyped-def]
    steps = load_steps(case)
    return steps[: next(i for i, step in enumerate(steps) if step.is_epsr_run)]


def _assert_identical(reference: Path, candidate: Path, extra: set[str]) -> None:
    comparison = compare_dirs(reference, candidate)
    assert comparison.all_identical, comparison.report()
    names = {path.name for path in candidate.iterdir() if path.is_file()}
    reference_names = {path.name for path in reference.iterdir() if path.is_file()}
    assert names - reference_names <= extra, names - reference_names


@pytest.fixture
def water_project(tmp_path: Path) -> Project:
    _needs_water()
    return Project.create(
        tmp_path / "water", water_system(), water_samples(), WATER_SETTINGS
    )


# --------------------------------------------------------------------- C8: create
def test_create_builds_the_validated_folder(
    water_project: Project, tmp_path: Path
) -> None:
    reference = _steps_folder(WATER, tmp_path / "steps", _setup_steps(WATER))
    (reference / "waterbox.EPSR.inp").unlink()  # the steps copy a hand-written one
    _assert_identical(
        reference, water_project.folder, API_EXTRA_FILES | {"waterbox.EPSR.inp"}
    )
    assert not water_project.has_run
    assert not (water_project.folder / "waterbox.pcof").exists()
    assert str(water_project.status()) == (
        "not run yet; 1000 molecules; 3 datasets (D2O, HDO, H2O); seed default"
    )


def test_create_silica_skips_fmole(tmp_path: Path) -> None:
    if not SILICA.available:
        pytest.skip("validation/silica_bulk or its data not available")
    project = Project.create(
        tmp_path / "silica", silica_system(), silica_samples(), SILICA_SETTINGS
    )
    reference = _steps_folder(SILICA, tmp_path / "steps", _setup_steps(SILICA))
    (reference / "silica.EPSR.inp").unlink()
    _assert_identical(
        reference, project.folder, {"epsr_project.json", "silica.EPSR.inp"}
    )


def test_project_file(water_project: Project) -> None:
    record = json.loads(water_project.project_file.read_text())
    assert record["seed"] is None and record["relax_iterations"] == 100
    assert [sample["name"] for sample in record["samples"]] == ["D2O", "HDO", "H2O"]
    assert water_project.system.to_dict() == water_system().to_dict()
    samples = water_project.samples
    assert [sample.name for sample in samples] == ["D2O", "HDO", "H2O"]
    assert samples[0].data_file == water_project.folder / "SLS18502-D2O.mint01"
    assert samples[1].isotopes == {"HW": {1: 0.5, 2: 0.5}}


@pytest.mark.parametrize(
    "change, error, message",
    [
        ({"settings": {"ereq": 3.0}}, ValueError, "refine"),
        ({"settings": {"niter": 5}}, KeyError, "cannot be set"),
        ({"settings": {"data": {"Null": {"tweak": 1.0}}}}, KeyError, "Null"),
        ({"samples": []}, ValueError, "at least one sample"),
        ({"seed": 0}, ValueError, "seed"),
        ({"seed": 1.5}, TypeError, "seed"),
        ({"relax_iterations": -1}, ValueError, "relax_iterations"),
    ],
)
def test_create_refuses(
    change: dict[str, object], error: type[Exception], message: str, tmp_path: Path
) -> None:
    _needs_water()
    arguments: dict[str, object] = {
        "system": water_system(),
        "samples": water_samples(),
        "settings": WATER_SETTINGS,
    }
    arguments.update(change)
    with pytest.raises(error, match=message):
        Project.create(tmp_path / "water", **arguments)  # type: ignore[arg-type]
    assert not (tmp_path / "water").exists()


def test_create_refuses_unknown_isotopes(tmp_path: Path) -> None:
    _needs_water()
    samples = water_samples()
    samples[0] = NeutronSample("D2O", samples[0].data_file, isotopes={"HW": {5: 1.0}})
    with pytest.raises(ValueError, match="mass number 5"):
        Project.create(tmp_path / "water", water_system(), samples)


def test_create_overwrite(water_project: Project) -> None:
    with pytest.raises(FileExistsError):
        Project.create(water_project.folder, water_system(), water_samples())
    Project.create(
        water_project.folder, water_system(), water_samples(), overwrite=True
    )


def test_create_with_a_seed(tmp_path: Path) -> None:
    _needs_water()
    boxes = {}
    for seed in (None, 1, 42):
        project = Project.create(
            tmp_path / f"seed{seed}", water_system(), water_samples(), seed=seed
        )
        boxes[seed] = (project.folder / "waterbox.ato").read_bytes()
        assert project.seed == seed
    assert boxes[None] == boxes[1]
    assert boxes[42] != boxes[1]
    unrandomised = Project.create(
        tmp_path / "lattice",
        water_system(),
        water_samples(),
        randomise=False,
        seed=42,
        relax_iterations=0,
    )
    rng_records = [
        line.split()
        for line in (unrandomised.folder / "waterbox.ato").read_text().splitlines()
        if len(line.split()) == 35
    ]
    assert rng_records[0][:3] == ["0", "-42", "0"]


# --------------------------------------------------------------------- C9-C12: running
def test_run_and_refine_need_the_scratch_start(water_project: Project) -> None:
    with pytest.raises(ProjectStateError, match="equilibrate"):
        water_project.run(iterations=2)
    with pytest.raises(ProjectStateError, match="equilibrate"):
        water_project.refine(ereq=3.0)
    with pytest.raises(ValueError, match="more than the scratch_iterations"):
        water_project.equilibrate(iterations=1)


def test_short_schedule_equals_the_validated_steps(
    water_project: Project, tmp_path: Path
) -> None:
    """create + equilibrate(1 + 2) + refine(ereq 3, 2 x 2 iterations), recorded."""
    equilibration = water_project.equilibrate(iterations=3, record=epsr.vis.available())
    assert len(equilibration) == 2  # the scratch start, then 2 iterations
    assert list(equilibration.iteration) == [1, 2]
    history = water_project.refine(
        ereq=3.0, runs=2, iterations=2, record=epsr.vis.available()
    )
    assert len(history) == 2 and history.ep_amplitude[-1] > 0
    steps = shorten_schedule(load_steps(WATER), 2, 2)
    reference = _steps_folder(WATER, tmp_path / "steps", steps)
    _assert_identical(reference, water_project.folder, API_EXTRA_FILES)
    status = water_project.status()
    assert status.has_run and status.iterations == 7 and status.ereq == 3.0
    if epsr.vis.available():
        assert water_project.trajectory().nsteps == 7


def test_pending_pcof_settings_reach_the_scratch_start(tmp_path: Path) -> None:
    _needs_water()
    project = Project.create(
        tmp_path / "water",
        water_system(),
        water_samples(),
        {"rminpt": 8.0, "rshmin": 0.7},
    )
    assert project.status().pending == {"rminpt": 8.0}
    project.equilibrate(iterations=2)
    assert float(project.settings["rminpt"]) == 8.0
    assert project.status().pending == {}


def test_set_edits_the_files(water_project: Project) -> None:
    water_project.set(ntimes=3, data={"HDO": {"tweak": 1.05}}, psigma2=0.02)
    assert water_project.input_file().get("ntimes") == "3"
    assert water_project.input_file().blocks[1]["tweak"] == "1.05"
    assert water_project.input_file().blocks[0]["tweak"] == "1.0"
    assert water_project.status().pending == {"psigma2": 0.02}
    with pytest.raises(ValueError, match="refine"):
        water_project.set(ereq=1.0)
    water_project.equilibrate(iterations=2)
    water_project.set(ereq=2.0)
    assert float(water_project.settings["ereq"]) == 2.0
    settings = water_project.settings
    assert settings["ntimes"] == "3" and settings["data"]["HDO"]["tweak"] == "1.05"
    assert "niter" not in settings and "fnameato" not in settings


def test_results_change_no_file(water_project: Project) -> None:
    water_project.equilibrate(iterations=2)
    water_project.run(iterations=1)  # so that a .qdr trend file exists
    before = _hashes(water_project.folder)
    result = water_project.results()
    assert _hashes(water_project.folder) == before
    assert result.summary == water_project.last_summary
    assert result.summary.r_factor_for("HDO") == result.summary.r_factor_per_dataset[1]


# --------------------------------------------------------------------- C14: datasets by name
def test_datasets_by_name(water_project: Project) -> None:
    water_project.equilibrate(iterations=2)
    result = water_project.results()
    assert result.dataset_names == ["D2O", "HDO", "H2O"]
    assert result.fit("HDO").name == "HDO"
    np.testing.assert_array_equal(result.fit("HDO").data, result.fit(1).data)
    np.testing.assert_array_equal(result.data_fq("H2O"), result.data_fq(-1))
    assert result.get_parameter("datafile", "HDO") == "SLS18500-HDO5050.mint01"
    with pytest.raises(KeyError, match="Null"):
        result.fit("Null")
    by_name = water_project.simulation(overrides={"data": {"HDO": {"tweak": 1.05}}})
    by_number = water_project.simulation(overrides={"data": {2: {"tweak": 1.05}}})
    assert by_name.overrides == by_number.overrides == ["data 2", "tweak 1.05"]
    checked = water_project.simulation(settings={"data": {1: {"tweak": 1.05}}})
    assert checked.overrides == ["data 2", "tweak 1.05"]
    as_overrides = water_project.simulation(
        epsr.Settings({"data": {1: {"tweak": 1.05}}})
    )
    assert as_overrides.overrides == ["data 2", "tweak 1.05"]
    with pytest.raises(TypeError, match="not both"):
        water_project.simulation(overrides={"ntimes": 1}, settings={"ntimes": 1})


def test_dataset_names_without_a_project_file(project: Project) -> None:
    simulation = project.simulation()
    assert simulation.dataset_names == [dataset.name for dataset in simulation.datasets]


# --------------------------------------------------------------------- C15, C16, C17
def test_pair_params_water(water_project: Project) -> None:
    water_project.equilibrate(iterations=2)
    table = water_project.results().potential.pair_params()
    expected = water_system().pair_parameters()
    for row, reference in zip(table, expected, strict=True):
        assert row.types == reference.types and row.terms == reference.terms
        assert (row.epsilon, row.sigma, row.charge_product) == pytest.approx(
            (reference.epsilon, reference.sigma, reference.charge_product), rel=1e-6
        )
        assert row.coulfac == 1.0 and row.extras == ""
    assert "coulfac" in str(table).splitlines()[0]


def test_pair_params_show_pcof_extras(tmp_path: Path) -> None:
    case = TEST_PROBLEMS / "200DMSO400Water" / "EPSR-test-ex1"
    if not case.exists():
        pytest.skip("test_problems not available")
    folder = copy_files(case, tmp_path / "ex1", CASE_FILES_TO_SKIP)
    table = epsr.Simulation(folder).potential.pair_params()
    assert table["S-S"].rminex == pytest.approx(4.0)
    assert table["S-S"].terms.endswith("+ repulsive core")
    assert "rminex 4" in table["S-S"].extras and "coordex 0.1" in table["S-S"].extras
    assert sum(row.rminex > 0 for row in table) == 3  # type: ignore[operator]


def test_analysis() -> None:
    r = np.linspace(0.0, 10.0, 1001)
    g = np.exp(-((r - 2.8) ** 2) / 0.02) + 0.5 * np.exp(-((r - 1.8) ** 2) / 0.01)
    assert epsr.analysis.first_peak(r, g, (2.0, 4.0)) == pytest.approx(2.8)
    assert epsr.analysis.first_peak(r, g, (1.5, 2.3)) == pytest.approx(1.8)
    flat = np.ones_like(r)
    expected = 4 / 3 * np.pi * 3.0**3 * 0.1
    assert epsr.analysis.coordination_number(r, flat, 0.1, 3.0) == pytest.approx(
        expected, 1e-4
    )
    with pytest.raises(ValueError, match="window"):
        epsr.analysis.first_peak(r, g, (20.0, 30.0))


def test_vis_available_without_the_plot_extra() -> None:
    code = (
        "import sys; sys.modules['matplotlib'] = None\n"
        "import epsr.vis as vis\n"
        "assert vis.available() is False\n"
        "try:\n    vis.TrajectoryRecorder\nexcept ImportError as e:\n    print('ok', e)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.startswith("ok epsr.vis needs the 'plot' extra")
    assert epsr.vis.available() == (
        __import__("importlib").util.find_spec("matplotlib") is not None
    )


# --------------------------------------------------------------------- the full schedule
@pytest.mark.slow
@pytest.mark.experimental
@pytest.mark.xdist_group("water")
def test_full_water_schedule_equals_the_committed_run(tmp_path: Path) -> None:
    """C9, C10, C12, C16: create + equilibrate(100) + refine(3.0, 2 x 100), recorded."""
    _needs_water()
    record = epsr.vis.available()
    project = Project.create(
        tmp_path / "water", water_system(), water_samples(), WATER_SETTINGS
    )
    equilibration = project.equilibrate(iterations=100, record=record)
    history = project.refine(ereq=3.0, runs=2, iterations=100, record=record)
    assert len(equilibration) == 2 and len(history) == 2
    committed = tmp_path / "committed"
    shutil.copytree(WATER_RUN, committed)
    (committed / "waterbox.EPSR.qdr").unlink()  # excepted, as for the example (plan)
    _assert_identical(
        committed, project.folder, API_EXTRA_FILES | {"waterbox.EPSR.qdr"}
    )
    if record:
        assert project.trajectory().nsteps == 1 + 99 + 200
    result = project.results()
    oxygen_oxygen = epsr.analysis.first_peak(
        result.r, result.gr_inter("OW", "OW"), (2.0, 4.0)
    )
    hydrogen_bond = epsr.analysis.first_peak(
        result.r, result.gr_inter("OW", "HW"), (1.5, 2.3)
    )
    assert oxygen_oxygen == pytest.approx(2.79, abs=0.03)
    assert hydrogen_bond == pytest.approx(1.8, abs=0.1)
