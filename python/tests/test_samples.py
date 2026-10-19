"""``epsr.samples``: neutron samples, their weights and the ``.mint01`` checks (Phase 3, C4)."""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr import tools
from epsr.samples import NeutronSample, check_sample_names
from epsr.system import AtomType, Molecule, System

from epsr_cases import REPOSITORY

WATER_CASE = REPOSITORY / "validation" / "water"
SILICA_CASE = REPOSITORY / "validation" / "silica_bulk"
WATER_DATA = REPOSITORY / "experimental_data" / "water" / "data"
SILICA_DATA = REPOSITORY / "experimental_data" / "silica_bulk" / "data"
WATER_RUN = (
    REPOSITORY / "python" / "examples" / "07_water_from_scratch" / "output" / "epsr_run"
)

OW = AtomType("OW", "O", 0.6503, 3.165492, -0.82, coarse_grain={"rmax": 1.55})
HW = AtomType("HW", "H", 0.0, 0.0, 0.41, 2.0, exchangeable=True)


def water_system() -> System:
    water = Molecule.from_geometry(
        "Water",
        centre=("O", OW),
        arms=[("H1", HW, 0.976), ("H2", HW, 0.976)],
        angle=107.134,
    )
    return System("waterbox", temperature=300.0, density=0.1).add(water, 1000)


def water_samples() -> list[NeutronSample]:
    return [
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


def _box(system: System, folder: Path) -> Path:
    """The box of ``system`` (makemole + mixato; nwts only needs the atom types)."""
    folder.mkdir(parents=True)
    for molecule in system.molecules:
        (folder / f"{molecule.name}.mol").write_text(
            molecule.render_mol(system.temperature)
        )
        tools.makemole(folder, molecule.name)
    tools.mixato(
        folder,
        [(molecule.name, count) for molecule, count in system.components],
        system.density,
        system.name,
    )
    return folder


@pytest.fixture(scope="module")
def water_box(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if not (WATER_CASE / "Water.mol").exists() or not WATER_DATA.exists():
        pytest.skip("validation/water or its data not available")
    return _box(water_system(), tmp_path_factory.mktemp("water") / "box")


# --------------------------------------------------------------------- acceptance
def test_water_weights_equal_the_validated_ones(
    water_box: Path, tmp_path: Path
) -> None:
    system = water_system()
    generated = shutil.copytree(water_box, tmp_path / "generated")
    validated = shutil.copytree(water_box, tmp_path / "validated")
    for sample in water_samples():
        weights = tools.nwts(
            generated,
            sample.nwts_file,
            box=system.name,
            components=sample.nwts_components(system),
        )
        shutil.copy2(WATER_CASE / sample.nwts_file, validated)
        tools.nwts(validated, sample.nwts_file)
        reference = validated / sample.weights_file
        assert weights.read_bytes() == reference.read_bytes(), sample.name
        assert weights.read_bytes() == (WATER_RUN / sample.weights_file).read_bytes()
        # nwts rewrites its input; the result is the validated file itself
        rewritten = (generated / sample.nwts_file).read_bytes()
        assert rewritten == (WATER_CASE / sample.nwts_file).read_bytes(), sample.name


def test_render_nwts_is_what_nwts_reads(water_box: Path, tmp_path: Path) -> None:
    system = water_system()
    folder = shutil.copytree(water_box, tmp_path / "rendered")
    sample = water_samples()[1]
    (folder / sample.nwts_file).write_text(sample.render_nwts(system))
    tools.nwts(folder, sample.nwts_file)
    assert (folder / sample.nwts_file).read_bytes() == (
        WATER_CASE / sample.nwts_file
    ).read_bytes()
    text = sample.render_nwts(system)
    assert "abundances  1 0.5 2 0.5" in text and "iexchange   1" in text


def test_silica_weights(tmp_path: Path) -> None:
    if not (SILICA_CASE / "SiO2.NWTS.dat").exists() or not SILICA_DATA.exists():
        pytest.skip("validation/silica_bulk or its data not available")
    si = AtomType("Si", "Si", 0.175, 1.03, 2.0, coarse_grain={"rmax": 1.11})
    oxygen = AtomType("O", "O", 0.165, 3.5, -1.0, coarse_grain={"rmax": 1.55})
    system = System("silica", 298.0, 0.0664)
    system.add(Molecule.atom(si), 1500).add(Molecule.atom(oxygen), 3000)
    sample = NeutronSample("SiO2", SILICA_DATA / "SiO2_NIMROD.mint01", rshmin=1.2)
    box = _box(system, tmp_path / "box")
    validated = shutil.copytree(box, tmp_path / "validated")
    weights = tools.nwts(
        box, sample.nwts_file, box="silica", components=sample.nwts_components(system)
    )
    shutil.copy2(SILICA_CASE / "SiO2.NWTS.dat", validated)
    tools.nwts(validated, "SiO2.NWTS.dat")
    assert weights.read_bytes() == (validated / "SiO2.NWTStot.wts").read_bytes()
    assert (box / "SiO2.NWTS.dat").read_bytes() == (
        validated / "SiO2.NWTS.dat"
    ).read_bytes()


# --------------------------------------------------------------------- the sample
def test_sample_defaults_and_keywords() -> None:
    d2o, _, h2o = water_samples()
    assert d2o.nwts_file == "D2O.NWTS.dat" and d2o.weights_file == "D2O.NWTStot.wts"
    assert h2o.composition("HW") == {0: 1.0}
    assert d2o.dataset_keywords() == {
        "datafile": "SLS18502-D2O.mint01",
        "wtsfile": "D2O.NWTStot.wts",
        "nrtype": 5,
        "rshmin": 0.7,
        "szeros": 0.0,
        "tweak": 1.0,
        "efilereq": 1.0,
        "porodadd": 0.0,
    }
    assert d2o.dataset_keywords(rshmin=0.9)["rshmin"] == 0.9
    own = NeutronSample("D2O", WATER_DATA / "SLS18502-D2O.mint01", rshmin=1.1)
    assert own.dataset_keywords(rshmin=0.9)["rshmin"] == 1.1
    assert d2o.data is not None and d2o.data.n_points == 997
    assert "D2O" in repr(d2o) and "HW 2 1" in repr(d2o)


def test_exchange_comes_from_the_atom_type_unless_overridden() -> None:
    system = water_system()
    components = water_samples()[0].nwts_components(system)
    assert [component["iexchange"] for component in components] == [0, 1]
    not_exchanging = NeutronSample(
        "D2O", WATER_DATA / "SLS18502-D2O.mint01", exchangeable=[]
    )
    assert [c["iexchange"] for c in not_exchanging.nwts_components(system)] == [0, 0]
    with pytest.raises(TypeError, match="list"):
        NeutronSample("D2O", WATER_DATA / "SLS18502-D2O.mint01", exchangeable="HW")


@pytest.mark.parametrize(
    "isotopes, error, message",
    [
        ({"HW": {1: 0.4, 2: 0.5}}, ValueError, "sum to 0.9"),
        ({"HW": {"D": 1.0}}, TypeError, "mass numbers are ints"),
        ({"HW": {m: 1 / 6 for m in range(6)}}, ValueError, "at most 5"),
        ({"HW": {2: 1.5}}, ValueError, "between 0 and 1"),
    ],
)
def test_bad_isotopes_raise_on_creation(
    isotopes: dict[str, dict[object, float]], error: type[Exception], message: str
) -> None:
    with pytest.raises(error, match=message):
        NeutronSample("D2O", WATER_DATA / "SLS18502-D2O.mint01", isotopes=isotopes)  # type: ignore[arg-type]


def test_isotopes_are_checked_against_the_table() -> None:
    system = water_system()
    tritium = NeutronSample(
        "T2O", WATER_DATA / "SLS18502-D2O.mint01", isotopes={"HW": {3: 1}}
    )
    tritium.check(system)  # H-3 is in the table
    unknown_mass = NeutronSample(
        "X", WATER_DATA / "SLS18502-D2O.mint01", isotopes={"HW": {4: 1}}
    )
    with pytest.raises(ValueError, match="mass number 4"):
        unknown_mass.check(system)
    oxygen_18 = NeutronSample(
        "X", WATER_DATA / "SLS18502-D2O.mint01", isotopes={"OW": {18: 1}}
    )
    oxygen_18.check(system)
    unknown_label = NeutronSample(
        "X", WATER_DATA / "SLS18502-D2O.mint01", isotopes={"DW": {2: 1}}
    )
    with pytest.raises(KeyError, match="not in the system"):
        unknown_label.nwts_components(system)


def test_sample_names_and_files(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="sample name"):
        NeutronSample("D2O sample", WATER_DATA / "SLS18502-D2O.mint01")
    blank = tmp_path / "my data.mint01"
    shutil.copy2(WATER_DATA / "SLS18502-D2O.mint01", blank)
    with pytest.raises(ValueError, match="blank"):
        NeutronSample("D2O", blank)
    with pytest.raises(FileNotFoundError):
        NeutronSample("D2O", tmp_path / "missing.mint01")
    with pytest.raises(ValueError, match="nrtype"):
        NeutronSample("D2O", WATER_DATA / "SLS18502-D2O.mint01", nrtype=7)
    d2o, hdo, _ = water_samples()
    check_sample_names([d2o, hdo])
    with pytest.raises(ValueError, match="two samples"):
        check_sample_names([d2o, d2o])
    same_file = NeutronSample("Other", WATER_DATA / "SLS18502-D2O.mint01")
    with pytest.raises(ValueError, match="same name"):
        check_sample_names([d2o, same_file])


def test_notes_on_the_q_range() -> None:
    d2o = water_samples()[0]
    assert d2o.notes() == []
    assert "nq * qstep = 30" in d2o.notes(q_max=30.0)[0]


def test_sample_dict_round_trip() -> None:
    hdo = water_samples()[1]
    rebuilt = NeutronSample.from_dict(hdo.to_dict())
    assert rebuilt.to_dict() == hdo.to_dict()
    assert rebuilt.isotopes == {"HW": {1: 0.5, 2: 0.5}}


# --------------------------------------------------------------------- .mint01
def _write_rows(path: Path, rows: list[str]) -> Path:
    path.write_text("# source\n# title\n" + "\n".join(rows) + "\n")
    return path


def test_read_mint01() -> None:
    data = epsr.io.read_mint01(WATER_DATA / "SLS18498-H2O.mint01")
    assert data.title.startswith("H2O in 1mm can")
    assert data.n_points == 997
    assert data.q_range == pytest.approx((0.175, 49.925))
    assert data.warnings() == []
    assert data.q.dtype == np.float64


@pytest.mark.parametrize(
    "rows, message",
    [
        ([f"{0.1 * i} 1.0" for i in range(20)], "three numbers"),
        ([f"{0.1 * i} 1.0 0.1" for i in range(5)], "only 5 data rows"),
        ([f"{0.01 * i} 1.0 0.1" for i in range(5001)], "at most 5000"),
        (["0 0 0"] * 20, "every row is zero"),
        ([f"{0.1 * i} x 0.1" for i in range(20)], "line 3"),
    ],
)
def test_read_mint01_refuses(rows: list[str], message: str, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=message):
        epsr.io.read_mint01(_write_rows(tmp_path / "bad.mint01", rows))


def test_mint01_warnings(tmp_path: Path) -> None:
    rows = [f"{0.1 * (20 - i)} 1.0 0.1" for i in range(20)]
    warnings = epsr.io.read_mint01(
        _write_rows(tmp_path / "odd.mint01", rows)
    ).warnings()
    assert any("increasing" in warning for warning in warnings)
    assert any("zero row" in warning for warning in warnings)
