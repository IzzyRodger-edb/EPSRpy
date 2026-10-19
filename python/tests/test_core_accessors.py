"""Engine accessors used by the setup API (plan Phase 0): keyword table, isotope
table, RNG seeding and pair parameters, at the level of ``epsr._core``."""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

import epsr
from epsr import _core, tools

from epsr_cases import CASE_FILES_TO_SKIP, REPOSITORY, TEST_PROBLEMS, copy_files

WATER_CASE = REPOSITORY / "validation" / "water"
WATER_RUN = (
    REPOSITORY / "python" / "examples" / "07_water_from_scratch" / "output" / "epsr_run"
)
#: EPSR-test with repulsive-core terms on three pairs of the ``.pcof``.
EX1_CASE = TEST_PROBLEMS / "200DMSO400Water" / "EPSR-test-ex1"

#: RNG state ``(mcall, idum, iy)`` after ``randomise`` of the water box (plan §1.7).
DEFAULT_STATE_AFTER_RANDOMISE = (6000, 1715014484, 745054426)
SEED_42_STATE_AFTER_RANDOMISE = (6000, 1163647977, 1007254037)
RNG_RECORD_LENGTH = 35  # mcall, idum, iy and the 32-entry shuffle table


# --------------------------------------------------------------------- keyword table
def _keyword_names(path: Path) -> list[str]:
    """The keywords of an EPSR keyword file, in file order (no title, ``data n`` or ``q``)."""
    names = []
    for line in path.read_text(encoding="latin-1").splitlines()[1:]:
        words = line.split()
        if not words or words[0] in ("data", "q"):
            continue
        names.append(words[0])
    return names


def _pcof_header_names(path: Path) -> list[str]:
    names = []
    for line in path.read_text(encoding="latin-1").splitlines():
        words = line.split()
        if words == ["q"]:
            break
        names.append(words[0])
    return names


def _routed(file: str) -> list[str]:
    return [
        symbol[0] for symbol in _core.epsr_keywords()["symbols"] if symbol[5] == file
    ]


def test_epsr_keyword_table_defaults_parse() -> None:
    table = _core.epsr_keywords()
    assert len(table["symbols"]) == table["marks"]["nsymbol"] == 73
    for name, value_type, _, default, description, file in table["symbols"]:
        assert file in ("inp", "pcof", "data"), name
        assert description, name
        words = default.split()
        assert words, name
        if value_type == "i":
            [int(word) for word in words]
        elif value_type == "r":
            [float(word) for word in words]
        else:
            assert value_type == "c", name


def test_epsr_keyword_routing_matches_the_files_epsr_writes() -> None:
    if not (WATER_CASE / "waterbox.EPSR.inp").exists():
        pytest.skip("validation/water not available")
    inp_names = _keyword_names(WATER_CASE / "waterbox.EPSR.inp")
    n_datasets = 3
    assert inp_names == _routed("inp") + _routed("data") * n_datasets
    assert _pcof_header_names(WATER_RUN / "waterbox.pcof") == _routed("pcof")
    assert "ereq" in _routed("pcof") and "rminpt" in _routed("pcof")


# --------------------------------------------------------------------- isotope table
def test_isotope_table() -> None:
    table = _core.isotope_table()
    lengths = {(symbol, mass_number): length for symbol, mass_number, length in table}
    assert lengths[("H", 0)] == pytest.approx(-0.3739)
    assert lengths[("H", 1)] == pytest.approx(-0.37406)
    assert lengths[("H", 2)] == pytest.approx(0.6671)
    assert ("O", 0) in lengths and ("O", 99) not in lengths
    assert len(table) == 364


# --------------------------------------------------------------------- RNG seed
@pytest.fixture(scope="module")
def mixed_water_box(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The water box after ``makemole`` and ``mixato`` (no random draw yet)."""
    if not (WATER_CASE / "Water.mol").exists():
        pytest.skip("validation/water not available")
    folder = tmp_path_factory.mktemp("water") / "mixed"
    folder.mkdir()
    shutil.copy2(WATER_CASE / "Water.mol", folder)
    tools.run(folder, "makemole", ["Water"])
    tools.run(folder, "mixato", answers=(WATER_CASE / "mixato.stdin").read_text())
    return folder


def _randomised_box(mixed: Path, destination: Path, seed: int | None) -> bytes:
    folder = shutil.copytree(mixed, destination)
    if seed is not None:
        _core.set_rng_seed(folder / "waterbox.ato", seed)
    tools.run(folder, "randomise", ["waterbox"])
    return (folder / "waterbox.ato").read_bytes()


def _rng_state(ato: bytes) -> tuple[int, int, int]:
    for line in ato.decode("latin-1").splitlines():
        words = line.split()
        if len(words) == RNG_RECORD_LENGTH and all(
            word.lstrip("-").isdigit() for word in words
        ):
            return int(words[0]), int(words[1]), int(words[2])
    raise AssertionError("no RNG record in the .ato file")


def test_seed_rewrites_only_the_rng_record(
    mixed_water_box: Path, tmp_path: Path
) -> None:
    folder = shutil.copytree(mixed_water_box, tmp_path / "box")
    before = (folder / "waterbox.ato").read_text().splitlines()
    _core.set_rng_seed(folder / "waterbox.ato", 7)
    after = (folder / "waterbox.ato").read_text().splitlines()
    assert len(before) == len(after)
    changed = [
        index
        for index, (old, new) in enumerate(zip(before, after, strict=True))
        if old != new
    ]
    assert len(changed) == 1
    assert after[changed[0]].split()[:3] == ["0", "-7", "0"]


def test_seed_default_is_seed_one(mixed_water_box: Path, tmp_path: Path) -> None:
    default = _randomised_box(mixed_water_box, tmp_path / "default", None)
    seed_one = _randomised_box(mixed_water_box, tmp_path / "one", 1)
    assert default == seed_one
    assert _rng_state(default) == DEFAULT_STATE_AFTER_RANDOMISE


def test_seed_42_is_reproducible_and_different(
    mixed_water_box: Path, tmp_path: Path
) -> None:
    first = _randomised_box(mixed_water_box, tmp_path / "first", 42)
    second = _randomised_box(mixed_water_box, tmp_path / "second", 42)
    default = _randomised_box(mixed_water_box, tmp_path / "default", None)
    assert first == second
    assert first != default
    assert _rng_state(first) == SEED_42_STATE_AFTER_RANDOMISE


@pytest.mark.parametrize("seed", [0, -1, -42, 2**31 - 1])
def test_invalid_seeds_are_rejected(
    mixed_water_box: Path, tmp_path: Path, seed: int
) -> None:
    folder = shutil.copytree(mixed_water_box, tmp_path / "box")
    contents = (folder / "waterbox.ato").read_bytes()
    with pytest.raises(ValueError, match="seed must be between 1 and 2147483646"):
        _core.set_rng_seed(folder / "waterbox.ato", seed)
    assert (folder / "waterbox.ato").read_bytes() == contents


def test_seed_of_a_missing_file(tmp_path: Path) -> None:
    with pytest.raises(epsr.EpsrError):
        _core.set_rng_seed(tmp_path / "missing.ato", 1)


# --------------------------------------------------------------------- pair parameters
def test_pair_params_are_lorentz_berthelot(project: epsr.Project) -> None:
    engine = project.simulation()._engine
    types = {values["label"]: values for values in engine.type_params()}
    pairs = engine.pair_params()
    assert [(pair["type1"], pair["type2"]) for pair in pairs] == engine.pair_labels()
    for pair in pairs:
        first, second = types[pair["type1"]], types[pair["type2"]]
        epsilon = np.sqrt(np.float32(first["epsilon"]) * np.float32(second["epsilon"]))
        sigma = np.float32(0.5) * (
            np.float32(first["sigma"]) + np.float32(second["sigma"])
        )
        charge_product = (
            np.float32(first["charge"]) * np.float32(second["charge"])
        ) * np.float32(1390.0)
        assert np.float32(pair["ecor"]) == epsilon
        assert np.float32(pair["dcor"]) == sigma
        assert np.float32(pair["qcor"]) == charge_product
        assert pair["coulfac"] == 1.0
        assert pair["gaussians"] == []


def test_pair_params_show_the_pcof_extras(tmp_path: Path) -> None:
    if not EX1_CASE.exists():
        pytest.skip("test_problems not available")
    folder = copy_files(EX1_CASE, tmp_path / "ex1", CASE_FILES_TO_SKIP)
    pairs = {
        (pair["type1"], pair["type2"]): pair
        for pair in epsr.Simulation(folder)._engine.pair_params()
    }
    core = {
        labels: (pair["rminex"], pair["pcoefex"], pair["coordex"])
        for labels, pair in pairs.items()
        if pair["rminex"] > 0
    }
    assert core[("S", "S")] == pytest.approx((4.0, 0.0, 0.1))
    assert len(core) == 3
    assert sorted(
        values for labels, values in core.items() if labels != ("S", "S")
    ) == [
        pytest.approx((1.5, 0.2, 0.05)),
        pytest.approx((2.6, 0.5, 0.0)),
    ]
    assert [pair["pcof_entry"] for pair in pairs.values()] == list(range(len(pairs)))
