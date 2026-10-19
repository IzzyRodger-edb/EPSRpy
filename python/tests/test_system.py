"""``epsr.system``: atom types, molecules and the box (plan Phase 2, C1-C3)."""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import pytest

import epsr
from epsr import tools
from epsr.system import (
    MAX_ATOM_TYPES,
    MAX_MOLECULES,
    AtomType,
    Molecule,
    PairParameters,
    System,
)

from epsr_cases import REPOSITORY, TEMPLATE

WATER_CASE = REPOSITORY / "validation" / "water"
SILICA_CASE = REPOSITORY / "validation" / "silica_bulk"

OW = AtomType(
    "OW",
    element="O",
    epsilon=0.6503,
    sigma=3.165492,
    charge=-0.82,
    coarse_grain={"rmax": 1.55},
)
HW = AtomType(
    "HW", element="H", epsilon=0.0, sigma=0.0, charge=0.41, mass=2.0, exchangeable=True
)
SI = AtomType("Si", "Si", 0.175, 1.03, 2.0, coarse_grain={"rmax": 1.11})
OX = AtomType("O", "O", 0.165, 3.5, -1.0, coarse_grain={"rmax": 1.55})


def water() -> Molecule:
    return Molecule.from_geometry(
        "Water",
        centre=("O", OW),
        arms=[("H1", HW, 0.976), ("H2", HW, 0.976)],
        angle=107.134,
    )


def _needs(path: Path) -> Path:
    if not path.exists():
        pytest.skip(f"{path} not available")
    return path


def _makemole_template(mol_text: str, name: str, folder: Path) -> bytes:
    folder.mkdir(parents=True)
    (folder / f"{name}.mol").write_text(mol_text)
    return tools.makemole(folder, name).read_bytes()


# --------------------------------------------------------------------- AtomType
def test_atom_type_defaults() -> None:
    assert SI.mass == pytest.approx(28.0855)
    assert AtomType("H", "H", 0, 0, 0).mass == pytest.approx(1.00797)
    assert HW.mass == 2.0 and HW.exchangeable and not OW.exchangeable
    assert OW.sphere == (1.0, 0.0, 1.55, 1.0)
    assert HW.sphere == (1.0, 0.0, 0.0, 1.0)
    # values equal to makemole's defaults are dropped, so equal types compare equal
    assert AtomType("OW", "O", 0.6503, 3.165492, -0.82, coarse_grain={"nsph": 1}) == (
        AtomType("OW", "O", 0.6503, 3.165492, -0.82)
    )
    assert len({OW, HW, OW}) == 2


def test_atom_type_docstring_states_the_mixing_rules() -> None:
    assert "Lorentz" in (AtomType.__doc__ or "")
    assert "Coulomb-only" in (AtomType.__doc__ or "")


@pytest.mark.parametrize(
    "arguments, message",
    [
        (("OWAT", "O", 0.1, 3.0, 0.0), "1-3 letters"),
        (("1W", "O", 0.1, 3.0, 0.0), "1-3 letters"),
        (("OW", "Xx", 0.1, 3.0, 0.0), "not a chemical symbol"),
        (("OW", "O", -0.1, 3.0, 0.0), "epsilon"),
        (("OW", "O", 0.1, -3.0, 0.0), "sigma"),
        (("OW", "O", 0.1, 3.0, 0.0, 0.0), "mass"),
    ],
)
def test_atom_type_validation(arguments: tuple[object, ...], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        AtomType(*arguments)  # type: ignore[arg-type]


def test_atom_type_coarse_grain_keys() -> None:
    with pytest.raises(ValueError, match="coarse_grain"):
        AtomType("OW", "O", 0.1, 3.0, 0.0, coarse_grain={"radius": 1.0})


# --------------------------------------------------------------------- Molecule: .mol files
def test_water_mol_is_byte_identical_to_the_validated_file() -> None:
    reference = _needs(WATER_CASE / "Water.mol").read_text()
    assert water().render_mol(300.0) == reference


def test_water_template_from_makemole(tmp_path: Path) -> None:
    reference = _needs(WATER_CASE / "Water.mol").read_text()
    generated = _makemole_template(water().render_mol(300.0), "Water", tmp_path / "new")
    validated = _makemole_template(reference, "Water", tmp_path / "validated")
    assert generated == validated


@pytest.mark.parametrize("atom_type, file_name", [(SI, "Si"), (OX, "O")])
def test_silica_atoms(atom_type: AtomType, file_name: str, tmp_path: Path) -> None:
    reference = _needs(SILICA_CASE / f"{file_name}.mol").read_text()
    molecule = Molecule.atom(atom_type)
    assert molecule.name == file_name and molecule.atom_names == [file_name]
    assert molecule.render_mol(298.0) == reference
    generated = _makemole_template(
        molecule.render_mol(298.0), file_name, tmp_path / "new"
    )
    validated = _makemole_template(reference, file_name, tmp_path / "validated")
    assert generated == validated


def test_from_geometry_coordinates() -> None:
    molecule = water()
    assert molecule.atom("O").position == (0.0, 0.0, 0.0)
    assert molecule.atom("H1").position == (0.976, 0.0, 0.0)
    angle = math.radians(107.134)
    expected = (0.976 * math.cos(angle), 0.976 * math.sin(angle), 0.0)
    assert molecule.atom("H2").position == expected
    assert molecule.angle_value("H1", "O", "H2") == 107.134
    assert molecule.bond_length("O", "H2") == 0.976


@pytest.mark.parametrize("n_arms, angle", [(3, 106.7), (4, None)])
def test_from_geometry_pyramid_and_tetrahedron(
    n_arms: int, angle: float | None
) -> None:
    hydrogen = AtomType("H", "H", 0, 0, 0.1)
    centre = AtomType("C", "C" if n_arms == 4 else "N", 0.3, 3.4, -0.4)
    molecule = Molecule.from_geometry(
        "Centre",
        centre=("X", centre),
        arms=[(f"H{i}", hydrogen, 1.09) for i in range(n_arms)],
        angle=angle,
    )
    expected = angle if angle is not None else math.degrees(math.acos(-1 / 3))
    measured = []
    for i in range(n_arms):
        assert math.dist(molecule.atom(f"H{i}").position, (0, 0, 0)) == pytest.approx(
            1.09
        )
        for j in range(i + 1, n_arms):
            raw = Molecule(  # measured from the coordinates, not the stored targets
                "Probe",
                [(atom.name, atom.type, atom.position) for atom in molecule.atoms],
                molecule.bonds,
            ).angle_value(f"H{i}", "X", f"H{j}")
            measured.append(raw)
    assert measured == pytest.approx([expected] * len(measured))
    assert len(molecule.angles) == len(measured)


def test_explicit_molecule_by_name_and_index() -> None:
    by_name = Molecule(
        "Water",
        atoms=[
            ("O", OW, (0, 0, 0)),
            ("H1", HW, water().atom("H1").position),
            ("H2", HW, water().atom("H2").position),
        ],
        bonds=[("O", "H1"), ("O", "H2")],
        angles=[("H1", "O", "H2")],
    )
    by_index = Molecule(
        "Water", [(a.name, a.type, a.position) for a in by_name.atoms], [(0, 1), (0, 2)]
    )
    assert by_name.bonds == by_index.bonds == ((0, 1), (0, 2))
    assert by_name.angles == by_index.angles == ((1, 0, 2),)
    lines = by_name.render_mol(300.0).splitlines()
    assert "bond 001 002    0.97600" in lines
    assert "angle 002 001 003  107.13400" in lines


def test_rotations_render() -> None:
    hydrogen = AtomType("H", "H", 0, 0, 0)
    carbon = AtomType("C", "C", 0.3, 3.4, 0)
    ethane_like = Molecule(
        "Chain",
        [
            ("C1", carbon, (0, 0, 0)),
            ("C2", carbon, (1.5, 0, 0)),
            ("H", hydrogen, (2.0, 1.0, 0)),
        ],
        [("C1", "C2"), ("C2", "H")],
        rotations=[("C1", "C2")],
    )
    assert "rot   1   2" in ethane_like.render_mol(300.0).splitlines()


@pytest.mark.parametrize(
    "build, error, message",
    [
        (lambda: Molecule("Bad name", [("O", OW, (0, 0, 0))]), ValueError, "file name"),
        (lambda: Molecule("W", []), ValueError, "no atoms"),
        (
            lambda: Molecule(
                "W", [("O", OW, (0, 0, 0)), ("O", HW, (1, 0, 0))], [(0, 1)]
            ),
            ValueError,
            "used twice",
        ),
        (
            lambda: Molecule(
                "W", [("O", OW, (0, 0, 0)), ("H", HW, (0.1, 0, 0))], [(0, 1)]
            ),
            ValueError,
            "closer than",
        ),
        (
            lambda: Molecule("W", [("O", OW, (0, 0, 0)), ("H", HW, (1, 0, 0))]),
            ValueError,
            "no bonds",
        ),
        (
            lambda: Molecule(
                "W",
                [("O", OW, (0, 0, 0)), ("H1", HW, (1, 0, 0)), ("H2", HW, (0, 3, 0))],
                [("O", "H1")],
            ),
            ValueError,
            "do not connect",
        ),
        (
            lambda: Molecule(
                "W",
                [("O", OW, (0, 0, 0)), ("H1", HW, (1, 0, 0)), ("H2", HW, (0, 1, 0))],
                [("O", "H1"), ("O", "H2")],
                angles=[],
            ),
            ValueError,
            "missing",
        ),
        (
            lambda: Molecule(
                "W",
                [("O", OW, (0, 0, 0)), ("H1", HW, (1, 0, 0)), ("H2", HW, (0, 1, 0))],
                [("O", "H1"), ("O", "H2")],
                angles=[("O", "H1", "H2")],
            ),
            ValueError,
            "bonded to the middle",
        ),
        (
            lambda: Molecule(
                "W",
                [("O", OW, (0, 0, 0)), ("H1", HW, (1, 0, 0)), ("H2", HW, (0, 1, 0))],
                [("O", "H1"), ("O", "H2")],
                rotations=[("H1", "H2")],
            ),
            ValueError,
            "not about a bond",
        ),
        (
            lambda: Molecule("W", [("O", OW, (0, 0, 0))], [("O", "X")]),
            KeyError,
            "no atom named",
        ),
        (
            lambda: Molecule("W", [("O", OW, (0, 0, 0))], dihedrals=[(0, 1, 2, 3)]),
            epsr.EpsrUnsupportedError,
            "dihedral",
        ),
        (
            lambda: Molecule(
                "W",
                [(f"A{i}", HW, (0.4 * i, 0, 0)) for i in range(1000)],
                [(i, i + 1) for i in range(999)],
            ),
            ValueError,
            "three digits",
        ),
    ],
)
def test_molecule_validation(
    build: object, error: type[Exception], message: str
) -> None:
    with pytest.raises(error, match=message):
        build()  # type: ignore[operator]


def test_molecule_atom_on_class_and_instance() -> None:
    assert Molecule.atom(SI, name="Silicon").name == "Silicon"
    assert water().atom("H1").name == "H1"
    assert water().atom(2).name == "H2"


def test_restraint_model_is_documented() -> None:
    assert "vibtemp" in (Molecule.__doc__ or "")
    assert "not force constants" in (Molecule.__doc__ or "")


def test_from_xyz(tmp_path: Path) -> None:
    molecule = water()
    lines = ["3", "water"]
    for atom in molecule.atoms:
        lines.append(f"{atom.type.element} " + " ".join(repr(v) for v in atom.position))
    path = tmp_path / "water.xyz"
    path.write_text("\n".join(lines) + "\n")
    from_file = Molecule.from_xyz(path, types={"O": OW, "H": HW}, name="Water")
    assert from_file.atom_names == ["O1", "H1", "H2"]
    assert from_file.bonds == ((0, 1), (0, 2))
    assert from_file.angle_value("H1", "O1", "H2") == pytest.approx(107.134)
    with pytest.raises(KeyError, match="no atom type"):
        Molecule.from_xyz(path, types={"O": OW})


@pytest.mark.parametrize(
    "path, temperature",
    [
        (WATER_CASE / "Water.mol", 300.0),
        (SILICA_CASE / "Si.mol", 298.0),
        (TEMPLATE / "WaterMolecule.mol", 300.0),
    ],
)
def test_from_mol_round_trip(path: Path, temperature: float) -> None:
    text = _needs(path).read_text()
    assert Molecule.from_mol(path).render_mol(temperature) == text


def test_from_mol_dmso_gives_the_same_template(tmp_path: Path) -> None:
    """DMSO.mol lists its potential lines in another order; makemole does not care."""
    path = _needs(TEMPLATE / "DMSO.mol")
    molecule = Molecule.from_mol(path)
    assert [t.label for t in molecule.atom_types] == ["S", "C", "O", "H"]
    assert len(molecule.rotations) == 2
    generated = _makemole_template(molecule.render_mol(300.0), "DMSO", tmp_path / "new")
    original = _makemole_template(path.read_text(), "DMSO", tmp_path / "original")
    assert generated == original


@pytest.mark.parametrize("name", ["WaterMolecule", "DMSO"])
def test_from_ato_matches_from_mol(name: str) -> None:
    from_mol = Molecule.from_mol(_needs(TEMPLATE / f"{name}.mol"))
    from_ato = Molecule.from_ato(TEMPLATE / f"{name}.ato")
    assert from_ato.atom_types == from_mol.atom_types
    assert from_ato.bonds == from_mol.bonds
    assert {frozenset(r) for r in from_ato.rotations} == {
        frozenset(r) for r in from_mol.rotations
    }
    for bond in from_mol.bonds:  # three significant digits in the .ato file
        assert from_ato.bond_length(*bond) == pytest.approx(
            from_mol.bond_length(*bond), rel=5e-3
        )
    for angle in from_mol.angles:
        assert from_ato.angle_value(*angle) == pytest.approx(
            from_mol.angle_value(*angle), abs=1.0
        )


def test_molecule_dict_round_trip() -> None:
    molecule = water()
    rebuilt = Molecule.from_dict(molecule.to_dict(), {"OW": OW, "HW": HW})
    assert rebuilt == molecule
    assert rebuilt.render_mol(300.0) == molecule.render_mol(300.0)


# --------------------------------------------------------------------- System
def test_water_system() -> None:
    system = System("waterbox", temperature=300.0, density=0.1).add(water(), 1000)
    assert system.n_atoms == 3000 and system.n_molecules == 1000
    assert system.box_size == pytest.approx(31.0723, abs=1e-4)
    assert [t.label for t in system.atom_types] == ["OW", "HW"]
    assert str(system).startswith("waterbox: 1000 × Water, 3000 atoms, box 31.07 Å")


def test_water_pair_parameters() -> None:
    table = System("waterbox", 300.0, 0.1).add(water(), 1000).pair_parameters()
    assert [row.label for row in table] == ["OW-OW", "OW-HW", "HW-HW"]
    oo, oh, hh = table
    assert (oo.epsilon, oo.sigma) == pytest.approx((0.6503, 3.165492))
    assert oo.charge_product == pytest.approx(934.636)
    assert (oh.epsilon, oh.sigma, oh.charge_product) == pytest.approx(
        (0.0, 1.582746, -467.318)
    )
    assert (hh.epsilon, hh.sigma, hh.charge_product) == pytest.approx(
        (0.0, 0.0, 233.659)
    )
    assert [row.terms for row in table] == [
        "LJ + Coulomb",
        "Coulomb only",
        "Coulomb only",
    ]
    assert table["HW-OW"] is oh and table[("OW", "HW")] is oh
    text = str(table).splitlines()
    assert text[0].split() == ["pair", "eps_ij", "sigma_ij", "1390*qi*qj", "terms"]
    assert text[1].split() == [
        "OW-OW",
        "0.6503",
        "3.165",
        "+934.6",
        "LJ",
        "+",
        "Coulomb",
    ]
    assert text[2].split() == ["OW-HW", "0", "1.583", "-467.3", "Coulomb", "only"]


def test_pair_terms() -> None:
    assert PairParameters.terms_of(0.1, 0.0) == "LJ only"
    assert PairParameters.terms_of(0.0, 0.0) == "none"


def test_pair_table_dataframe() -> None:
    pytest.importorskip("pandas")
    frame = (
        System("waterbox", 300.0, 0.1).add(water(), 10).pair_parameters().to_dataframe()
    )
    assert list(frame["terms"]) == ["LJ + Coulomb", "Coulomb only", "Coulomb only"]


def test_atom_type_order_follows_add_order() -> None:
    system = System("silica", 298.0, 0.0664)
    system.add(Molecule.atom(SI), 1500).add(Molecule.atom(OX), 3000)
    assert [t.label for t in system.atom_types] == ["Si", "O"]
    assert [row.label for row in system.pair_parameters()] == ["Si-Si", "Si-O", "O-O"]
    assert system.box_size == pytest.approx(40.7707, abs=1e-4)


def test_system_validation() -> None:
    with pytest.raises(ValueError, match="density"):
        System("box", 300.0, 0.0)
    with pytest.raises(ValueError, match="temperature"):
        System("box", -1.0, 0.1)
    with pytest.raises(epsr.EpsrUnsupportedError, match="non-cubic"):
        System("box", 300.0, 0.1, box=(10, 10, 12))
    system = System("box", 300.0, 0.1).add(water(), 10)
    with pytest.raises(ValueError, match="already in the system"):
        system.add(water(), 10)
    with pytest.raises(ValueError, match="count"):
        system.add(Molecule.atom(SI), 0)
    with pytest.raises(ValueError, match="system's name"):
        system.add(Molecule.atom(SI, name="box"), 1)
    other_oxygen = AtomType("OW", "O", 0.5, 3.0, -0.8)
    with pytest.raises(ValueError, match="differs"):
        system.add(Molecule.atom(other_oxygen, name="Ion"), 1)
    with pytest.raises(ValueError, match=str(MAX_MOLECULES)):
        system.add(Molecule.atom(SI), MAX_MOLECULES)


def test_atom_type_limit() -> None:
    system = System("box", 300.0, 0.1)
    for index in range(MAX_ATOM_TYPES):
        label = f"A{index}"
        system.add(Molecule.atom(AtomType(label, "Ar", 0.1, 3.4, 0.0)), 1)
    with pytest.raises(ValueError, match=f"at most {MAX_ATOM_TYPES}"):
        system.add(Molecule.atom(AtomType("B0", "Ar", 0.1, 3.4, 0.0)), 1)


def test_system_dict_round_trip() -> None:
    system = System("waterbox", 300.0, 0.1).add(water(), 1000)
    rebuilt = System.from_dict(system.to_dict())
    assert rebuilt.to_dict() == system.to_dict()
    assert rebuilt.molecules[0] == system.molecules[0]


def test_top_level_exports() -> None:
    assert (
        epsr.AtomType is AtomType
        and epsr.Molecule is Molecule
        and epsr.System is System
    )


def test_makemole_with_the_rendered_file_in_a_fresh_folder(tmp_path: Path) -> None:
    """The template written from Python feeds mixato like the validated one."""
    _needs(WATER_CASE / "Water.mol")
    folder = tmp_path / "box"
    folder.mkdir()
    (folder / "Water.mol").write_text(water().render_mol(300.0))
    tools.makemole(folder, "Water")
    tools.mixato(folder, [("Water", 1000)], 0.1, "waterbox")
    reference = tmp_path / "reference"
    reference.mkdir()
    shutil.copy2(WATER_CASE / "Water.mol", reference)
    tools.makemole(reference, "Water")
    tools.mixato(reference, [("Water", 1000)], 0.1, "waterbox")
    assert (folder / "waterbox.ato").read_bytes() == (
        reference / "waterbox.ato"
    ).read_bytes()
