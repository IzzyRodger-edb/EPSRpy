"""The typed setup helpers of ``epsr.tools`` (plan Phase 1: makemole, mixato,
randomise with a seed, fmole, nwts from components)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

import epsr
from epsr import tools

from epsr_cases import REPOSITORY
from epsr_workflow import WATER, load_steps, prepare_case_folder, run_steps_with_wrapper

WATER_CASE = REPOSITORY / "validation" / "water"
#: The committed example 07 run, identical to validation/water (plan constraint 1).
WATER_RUN = (
    REPOSITORY / "python" / "examples" / "07_water_from_scratch" / "output" / "epsr_run"
)
#: Isotopic composition of HW in each water sample, in validation/water's dataset order.
HYDROGEN = {"D2O": {2: 1.0}, "HDO": {1: 0.5, 2: 0.5}, "H2O": {0: 1.0}}
#: The files the setup programs write (before the first EPSR run).
SETUP_FILES = [
    "Water.ato",
    "waterbox.ato",
    *(f"{sample}.NWTS.dat" for sample in HYDROGEN),
    *(f"{sample}.NWTStot.wts" for sample in HYDROGEN),
]
#: Setup files that no EPSR run rewrites, so they can be compared with the final run.
FILES_UNCHANGED_BY_EPSR = [name for name in SETUP_FILES if name != "waterbox.ato"]


def _water_components(sample: str) -> list[dict[str, object]]:
    return [
        {"atom": "OW", "atsymbol": "O", "iexchange": 0, "abundances": {0: 1.0}},
        {
            "atom": "HW",
            "atsymbol": "H",
            "iexchange": True,
            "abundances": HYDROGEN[sample],
        },
    ]


def _build_with_typed_tools(folder: Path, seed: int | None = None) -> Path:
    folder.mkdir(parents=True)
    shutil.copy2(WATER_CASE / "Water.mol", folder)
    tools.makemole(folder, "Water")
    tools.mixato(folder, [("Water", 1000)], density=0.1, box="waterbox")
    tools.randomise(folder, "waterbox", seed=seed)
    tools.fmole(folder, "waterbox", 100)
    for sample in HYDROGEN:
        tools.nwts(
            folder,
            f"{sample}.NWTS.dat",
            box="waterbox",
            components=_water_components(sample),
        )
    return folder


@pytest.fixture(scope="module")
def typed_water(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if not WATER.available:
        pytest.skip("validation/water or its data not available")
    return _build_with_typed_tools(tmp_path_factory.mktemp("typed") / "water")


def test_water_setup_matches_the_validated_steps(
    typed_water: Path, tmp_path: Path
) -> None:
    folder = prepare_case_folder(WATER, tmp_path / "steps")
    setup_steps = []
    for step in load_steps(WATER):
        if step.is_epsr_run:
            break
        setup_steps.append(step)
    for _ in run_steps_with_wrapper(setup_steps, folder, WATER.case_folder):
        pass
    for name in SETUP_FILES:
        assert (typed_water / name).read_bytes() == (folder / name).read_bytes(), name


def test_water_setup_matches_the_committed_run(typed_water: Path) -> None:
    for name in FILES_UNCHANGED_BY_EPSR:
        assert (typed_water / name).read_bytes() == (
            WATER_RUN / name
        ).read_bytes(), name


def test_fmole_puts_the_mol_file_back(typed_water: Path) -> None:
    assert (typed_water / "Water.mol").is_file()
    assert not list(typed_water.glob(".fmole-*"))


def test_randomise_with_a_seed(typed_water: Path, tmp_path: Path) -> None:
    folder = tmp_path / "box"
    folder.mkdir()
    shutil.copy2(typed_water / "Water.ato", folder)
    tools.mixato(folder, [("Water", 1000)], density=0.1, box="waterbox")
    mixed = (folder / "waterbox.ato").read_bytes()
    boxes = {}
    for seed in (None, 1, 42):
        (folder / "waterbox.ato").write_bytes(mixed)
        boxes[seed] = tools.randomise(folder, "waterbox", seed=seed).read_bytes()
    assert boxes[None] == boxes[1]
    assert boxes[42] != boxes[1]
    with pytest.raises(TypeError):
        tools.randomise(folder, "waterbox", seed=True)
    with pytest.raises(ValueError, match="seed"):
        tools.randomise(folder, "waterbox", seed=0)


# --------------------------------------------------------------------- errors
def test_makemole_needs_the_mol_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="Water.mol"):
        tools.makemole(tmp_path, "Water")


def test_mixato_checks_its_input(typed_water: Path, tmp_path: Path) -> None:
    shutil.copy2(typed_water / "Water.ato", tmp_path)
    with pytest.raises(FileNotFoundError, match="Ice.ato"):
        tools.mixato(tmp_path, [("Ice", 10)], density=0.1, box="box")
    with pytest.raises(ValueError, match="overwrite the template"):
        tools.mixato(tmp_path, [("Water", 10)], density=0.1, box="Water")
    with pytest.raises(ValueError, match="count"):
        tools.mixato(tmp_path, [("Water", 0)], density=0.1, box="box")
    with pytest.raises(ValueError, match="density"):
        tools.mixato(tmp_path, [("Water", 10)], density=0.0, box="box")
    with pytest.raises(ValueError, match="one word"):
        tools.mixato(tmp_path, [("Water", 10)], density=0.1, box="water box")


def test_fmole_checks_its_input(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        tools.fmole(tmp_path, "waterbox", 10)
    (tmp_path / "waterbox.ato").write_text("")
    with pytest.raises(ValueError, match="iterations"):
        tools.fmole(tmp_path, "waterbox", 0)


def test_abundances_text() -> None:
    assert tools.abundances_text({0: 1.0}) == "0 1"
    assert tools.abundances_text({2: 1}) == "2 1"
    assert tools.abundances_text({1: 0.5, 2: 0.5}, "H") == "1 0.5 2 0.5"
    assert tools.abundances_text({2: 0.36, 1: 0.64}) == "2 0.36 1 0.64"
    with pytest.raises(ValueError, match="sum to 0.9"):
        tools.abundances_text({1: 0.4, 2: 0.5})
    with pytest.raises(ValueError, match="mass number 4"):
        tools.abundances_text({4: 1.0}, "H")
    with pytest.raises(ValueError, match="not in the isotope table"):
        tools.abundances_text({0: 1.0}, "Xx")
    with pytest.raises(ValueError, match="at most 5"):
        tools.abundances_text({mass: 1 / 6 for mass in range(6)})
    with pytest.raises(TypeError):
        tools.abundances_text({"D": 1.0})  # type: ignore[dict-item]


def test_isotopes() -> None:
    hydrogen = tools.isotopes("H")
    assert sorted(hydrogen) == [0, 1, 2, 3]
    assert hydrogen[2] == pytest.approx(0.6671)
    assert tools.isotopes("Xx") == {}


def test_nwts_with_components_checks_its_input(tmp_path: Path) -> None:
    components = _water_components("D2O")
    with pytest.raises(ValueError, match="box="):
        tools.nwts(tmp_path, "D2O.NWTS.dat", components=components)
    with pytest.raises(ValueError, match=r"<sample>\.NWTS\.dat"):
        tools.nwts(tmp_path, "D2O.dat", box="waterbox", components=components)
    with pytest.raises(ValueError, match="mass number 7"):
        tools.nwts(
            tmp_path,
            "D2O.NWTS.dat",
            box="waterbox",
            components=[{"atom": "HW", "atsymbol": "H", "abundances": {7: 1.0}}],
        )
    with pytest.raises(epsr.EpsrIOError, match="waterbox.ato"):
        tools.nwts(tmp_path, "D2O.NWTS.dat", box="waterbox", components=components)
