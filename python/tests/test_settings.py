"""``epsr.settings``: EPSR's keyword table, ``Settings`` and the ``.inp`` (Phase 4, C5)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

import epsr
from epsr import tools
from epsr.samples import NeutronSample
from epsr.settings import (
    MANAGED_KEYWORDS,
    Settings,
    keywords,
    suggest_rshmin,
    value_text,
)
from epsr.system import AtomType, Molecule, System

from epsr_cases import REPOSITORY, compare_dirs
from epsr_workflow import WATER, load_steps, prepare_case_folder, run_steps_with_wrapper

WATER_DATA = REPOSITORY / "experimental_data" / "water" / "data"


# --------------------------------------------------------------------- the table
def test_keyword_table() -> None:
    table = keywords()
    assert len(table) == 73
    assert list(table)[:3] == ["feedback", "potfac", "ereq"]
    assert table["ereq"].file == "pcof" and table["ntimes"].file == "inp"
    assert table["rshmin"].file == "data" and table["rshmin"].default == "0.7"
    assert table["qmin"].position == 46
    assert Settings.describe("feedback") == (
        "feedback (.EPSR.inp, real, default 0.9): Confidence factor - should be < 1. [0.8]"
    )
    assert Settings.file_of("rminpt") == "pcof"


# --------------------------------------------------------------------- values
def test_value_text() -> None:
    assert value_text("ntimes", 10) == "10"
    assert value_text("ereq", 3) == "3.0"
    assert value_text("rshmin", 0.7) == "0.7"
    assert value_text("qwidthqmax", (0.2, 5.0)) == "0.2 5.0"
    assert value_text("hklqmin", [0.0, 10, 0.0]) == "0.0 10.0 0.0"
    assert value_text("sfreq", "50 1") == "50 1"
    assert value_text("control-p", "none 0 0") == "none 0 0"
    assert value_text("reppottype", "harmonic") == "harmonic"


@pytest.mark.parametrize(
    "name, value, error, message",
    [
        ("ntimes", 1.5, TypeError, "integer"),
        ("ntimes", "many", ValueError, "integer"),
        ("ntimes", True, TypeError, "not a valid value"),
        ("feedback", "high", ValueError, "real"),
        ("feedback", (), ValueError, "empty"),
        ("feedbak", 0.9, KeyError, "did you mean 'feedback'"),
    ],
)
def test_bad_values(
    name: str, value: object, error: type[Exception], message: str
) -> None:
    with pytest.raises(error, match=message):
        Settings({name: value})


@pytest.mark.parametrize("name", sorted(MANAGED_KEYWORDS))
def test_managed_keywords_are_refused(name: str) -> None:
    with pytest.raises(KeyError, match="cannot be set"):
        Settings({name: 1})


def test_settings_mapping() -> None:
    settings = Settings({"ntimes": 10, "rshmin": 0.7}, ereq=0.0, qwidthqmax=[0.2, 5.0])
    assert settings["qwidthqmax"] == (0.2, 5.0)
    assert settings.part("inp") == {"ntimes": 10, "qwidthqmax": (0.2, 5.0)}
    assert settings.part("pcof") == {"ereq": 0.0}
    assert settings.part("data") == {"rshmin": 0.7}
    settings["data"] = {"D2O": {"tweak": 1.05}, 2: {"rshmin": 0.9}}
    assert set(settings) == {"ntimes", "rshmin", "ereq", "qwidthqmax", "data"}
    assert settings.dataset_values(0, "D2O") == {"rshmin": 0.7, "tweak": 1.05}
    assert settings.dataset_values(2, "H2O") == {"rshmin": 0.9}
    with pytest.raises(KeyError, match="not a per-dataset keyword"):
        settings["data"] = {"D2O": {"ntimes": 3}}
    del settings["ntimes"]
    assert "ntimes" not in settings
    assert Settings.from_dict(settings.to_dict()).to_dict() == settings.to_dict()


def test_override_lines() -> None:
    settings = Settings({"ntimes": 10, "ereq": 3.0})
    assert settings.override_lines() == ["ntimes 10", "ereq 3.0"]
    settings["rshmin"] = 0.8
    settings["data"] = {"HDO": {"tweak": 1.05}}
    with pytest.raises(ValueError, match="dataset names"):
        settings.override_lines()
    assert settings.override_lines(["D2O", "HDO"]) == [
        "ntimes 10",
        "ereq 3.0",
        "data 1",
        "rshmin 0.8",
        "data 2",
        "rshmin 0.8",
        "tweak 1.05",
    ]
    settings["data"] = {"Null": {"tweak": 1.0}}
    with pytest.raises(KeyError, match="no dataset named 'Null'"):
        settings.override_lines(["D2O", "HDO"])


def test_render_inp_layout() -> None:
    datasets = [{"datafile": "a.mint01", "wtsfile": "A.NWTStot.wts", "rshmin": 0.9}]
    text = Settings({"ntimes": 7}).render_inp("box.EPSR", "box", datasets)
    lines = text.splitlines()
    assert lines[0] == "box.EPSR               Title of this file"
    names = [line.split()[0] for line in lines[1:] if line.strip()]
    inp_names = [k.name for k in keywords().values() if k.file == "inp"]
    data_names = [k.name for k in keywords().values() if k.file == "data"]
    assert names == [*inp_names, "data", *data_names, "q"]
    values = {
        line.split()[0]: line.split()[1] for line in lines[1:] if len(line.split()) > 1
    }
    assert values["ntimes"] == "7" and values["fnameato"] == "box.ato"
    assert values["fnamepcof"] == "box.pcof" and values["ndata"] == "1"
    assert (values["ireset"], values["iinit"], values["nsumt"], values["niter"]) == (
        "1",
        "1",
        "-1",
        "1",
    )
    assert values["rshmin"] == "0.9" and values["tweak"] == "1.0"
    assert "ereq" not in values  # .pcof keywords are not in the .inp


# --------------------------------------------------------------------- acceptance
def test_scratch_run_from_the_rendered_inp_equals_the_validated_one(
    tmp_path: Path,
) -> None:
    if not WATER.available:
        pytest.skip("validation/water or its data not available")
    oxygen = AtomType("OW", "O", 0.6503, 3.165492, -0.82, coarse_grain={"rmax": 1.55})
    hydrogen = AtomType("HW", "H", 0.0, 0.0, 0.41, 2.0, exchangeable=True)
    water = Molecule.from_geometry(
        "Water",
        centre=("O", oxygen),
        arms=[("H1", hydrogen, 0.976), ("H2", hydrogen, 0.976)],
        angle=107.134,
    )
    system = System("waterbox", 300.0, 0.1).add(water, 1000)
    samples = [
        NeutronSample(
            "D2O", WATER_DATA / "SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}
        ),
        NeutronSample(
            "HDO",
            WATER_DATA / "SLS18500-HDO5050.mint01",
            isotopes={"HW": {1: 0.5, 2: 0.5}},
        ),
        NeutronSample("H2O", WATER_DATA / "SLS18498-H2O.mint01"),
    ]
    settings = Settings({"rshmin": 0.7})

    folder = tmp_path / "new"
    folder.mkdir()
    for sample in samples:
        shutil.copy2(sample.data_file, folder)
    (folder / "Water.mol").write_text(water.render_mol(system.temperature))
    tools.makemole(folder, "Water")
    tools.mixato(folder, [("Water", 1000)], system.density, system.name)
    tools.randomise(folder, system.name)
    (folder / "Water.mol").unlink()  # as validation/water/steps.txt
    tools.fmole(folder, system.name, 100)
    for sample in samples:
        tools.nwts(
            folder,
            sample.nwts_file,
            box=system.name,
            components=sample.nwts_components(system),
        )
    datasets = [sample.dataset_keywords(settings.part("data")) for sample in samples]
    (folder / "waterbox.EPSR.inp").write_text(
        settings.render_inp("waterbox.EPSR", system.name, datasets)
    )
    epsr.Simulation(folder, "waterbox.EPSR.inp").run()

    reference = prepare_case_folder(WATER, tmp_path / "steps")
    steps = load_steps(WATER)
    first_run = next(index for index, step in enumerate(steps) if step.is_epsr_run)
    for _ in run_steps_with_wrapper(
        steps[: first_run + 1], reference, WATER.case_folder
    ):
        pass
    (
        reference / "mixato.stdin"
    ).unlink()  # a prompt file of the steps, not an EPSR file
    comparison = compare_dirs(reference, folder)
    assert comparison.all_identical, comparison.report()
    assert sorted(path.name for path in folder.iterdir()) == sorted(
        path.name for path in reference.iterdir()
    )


def test_suggest_rshmin() -> None:
    hydrogen = AtomType("H", "H", 0, 0, 0)
    argon = AtomType("Ar", "Ar", 0.99, 3.4, 0)
    assert (
        suggest_rshmin(System("argon", 85.0, 0.02).add(Molecule.atom(argon), 10)) == 1.2
    )
    dimer = Molecule(
        "H2", [("A", hydrogen, (0, 0, 0)), ("B", hydrogen, (0.6, 0, 0))], [(0, 1)]
    )
    assert suggest_rshmin(System("hydrogen", 20.0, 0.04).add(dimer, 10)) == 0.4
    assert epsr.keywords() == keywords() and epsr.Settings is Settings
