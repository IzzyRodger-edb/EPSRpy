"""Defining a system in Python: atom types, molecules and the simulation box.

These are plain value objects. They describe the system; :meth:`epsr.Project.create`
turns them into the files EPSR's own setup programs read (``.mol`` for
``makemole``, the answers of ``mixato``) and runs those programs, so the result is
byte for byte what EPSR would build from the same numbers.

>>> OW = AtomType("OW", element="O", epsilon=0.6503, sigma=3.165492, charge=-0.82,
...               coarse_grain={"rmax": 1.55})
>>> HW = AtomType("HW", element="H", epsilon=0.0, sigma=0.0, charge=0.41, mass=2.0,
...               exchangeable=True)
>>> water = Molecule.from_geometry("Water", centre=("O", OW),
...                                arms=[("H1", HW, 0.976), ("H2", HW, 0.976)],
...                                angle=107.134)
>>> system = System("waterbox", temperature=300.0, density=0.1).add(water, 1000)

**How EPSR models the interactions** (see ``docs/plans/python_setup_api.md`` §1.8):

* *Between molecules*, every pair of atom types interacts through the reference
  potential mixed from the two types (Lorentz–Berthelot, Lennard-Jones + Coulomb;
  :meth:`System.pair_parameters`), plus the empirical potential once ``ereq > 0``.
* *Within a molecule*, bonds and angles are harmonic restraints to the target
  geometry, with a stiffness set per molecule by ``vibtemp``/``angtemp``; see
  :class:`Molecule`.
"""

from __future__ import annotations

import functools
import math
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import Any

from . import _format
from ._elements import COVALENT_RADII, ELEMENT_MASSES
from .errors import EpsrUnsupportedError

__all__ = [
    "AtomType",
    "Atom",
    "Molecule",
    "System",
    "PairParameters",
    "PairTable",
    "LABEL_MAX_LENGTH",
    "MAX_ATOM_TYPES",
    "MAX_MOLECULES",
    "MAX_ATOMS_PER_MOLECULE",
]

#: ``makemole`` and the ``.ato`` file keep three characters of an atom-type label.
LABEL_MAX_LENGTH = 3
#: Engine limits (``state.rs``): atom types (``MCOMP``) and molecules (``MMOL``).
MAX_ATOM_TYPES = 40
MAX_MOLECULES = 50_000
#: ``makemole`` labels the atoms of a molecule ``001`` to ``999``.
MAX_ATOMS_PER_MOLECULE = 999
#: Coulomb prefactor of the engine: 1390 q_i q_j / r in kJ/mol with r in Å.
COULOMB_FACTOR = 1390.0
#: Atoms of a molecule closer than this (Å) are refused as a likely input error.
MIN_ATOM_SEPARATION = 0.3
#: ``Molecule.from_xyz(bonds="auto")``: bonded if closer than r_i + r_j + this (Å).
BOND_TOLERANCE = 0.4
#: Density written to a molecule template (``.mol`` ``density`` line); the box has its own.
TEMPLATE_DENSITY = 0.001

_LABEL_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9]{0,2}$")
_NAME_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_-]{0,23}$")
#: The four trailing values of a ``.mol`` ``potential`` line and their names.
COARSE_GRAIN_KEYS = ("nsph", "rmin", "rmax", "cg")

#: An atom of a molecule is referred to by its name or its 0-based index.
AtomRef = str | int


def _check_file_name(name: str, what: str) -> str:
    if not isinstance(name, str) or not _NAME_PATTERN.match(name):
        raise ValueError(
            f"{what} {name!r}: use 1-24 letters, digits, '-' or '_', starting with "
            "a letter (it becomes a file name)"
        )
    return name


def _positive(value: float, what: str, *, allow_zero: bool = False) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        relation = ">= 0" if allow_zero else "> 0"
        raise ValueError(f"{what} must be {relation}, got {value!r}")
    return value


# --------------------------------------------------------------------------- AtomType
@dataclass(frozen=True)
class AtomType:
    """An EPSR atom type and its reference potential.

    Parameters
    ----------
    label:
        At most three letters or digits, starting with a letter (``makemole``
        and the ``.ato`` file keep three characters, so longer labels would
        silently collide). Unique within a :class:`System`.
    element:
        Chemical symbol (``"O"``, ``"Cu"``). It selects the neutron scattering
        lengths and the default mass.
    epsilon, sigma:
        Lennard-Jones well depth (kJ/mol) and diameter (Å).
    charge:
        Point charge (e).
    mass:
        Mass in amu; by default EPSR's atomic weight of ``element``.
    exchangeable:
        True for labile hydrogens that exchange between molecules (``iexchange 1``
        in every sample's weights; a :class:`~epsr.samples.NeutronSample` may
        override it).
    coarse_grain:
        The coarse-graining sphere of the ``.mol`` ``potential`` line: any of
        ``nsph``, ``rmin``, ``rmax`` and ``cg``. Missing values take ``makemole``'s
        defaults ``nsph=1, rmin=0, rmax=0.5*sigma, cg=1``.

    Notes
    -----
    **Pairs are always mixed from the two types.** EPSR derives the reference
    potential of every pair of types i, j by Lorentz–Berthelot:
    ε_ij = √(ε_i ε_j), σ_ij = (σ_i + σ_j)/2, plus the Coulomb term
    1390 q_i q_j / r (kJ/mol, r in Å). There is no per-pair Lennard-Jones input.
    The Lennard-Jones term is only computed when ε_ij > 0, so **a type with
    ε = 0 has no Lennard-Jones site**: all its pairs are Coulomb-only (the SPC
    convention for water hydrogens). :meth:`System.pair_parameters` shows the
    resulting terms of each pair. The only per-pair adjustments EPSR has
    (``coulfac``, an extra repulsive core, additional Gaussian terms) live in the
    ``.pcof`` file and are not set from here.
    """

    label: str
    element: str
    epsilon: float
    sigma: float
    charge: float
    mass: float = None  # type: ignore[assignment]  # filled in from the element
    exchangeable: bool = field(default=False, kw_only=True)
    coarse_grain: Mapping[str, float] = field(default_factory=dict, kw_only=True)

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not _LABEL_PATTERN.match(self.label):
            raise ValueError(
                f"atom type label {self.label!r}: use 1-{LABEL_MAX_LENGTH} letters or "
                "digits starting with a letter (EPSR keeps three characters)"
            )
        if self.element not in ELEMENT_MASSES:
            raise ValueError(
                f"atom type {self.label}: {self.element!r} is not a chemical symbol"
            )
        set_field = functools.partial(object.__setattr__, self)
        set_field(
            "epsilon",
            _positive(self.epsilon, f"{self.label}: epsilon", allow_zero=True),
        )
        set_field(
            "sigma", _positive(self.sigma, f"{self.label}: sigma", allow_zero=True)
        )
        set_field("charge", float(self.charge))
        if not math.isfinite(self.charge):
            raise ValueError(f"{self.label}: charge must be finite")
        mass = ELEMENT_MASSES[self.element] if self.mass is None else self.mass
        set_field("mass", _positive(mass, f"{self.label}: mass"))
        set_field("exchangeable", bool(self.exchangeable))
        set_field("coarse_grain", self._normalised_coarse_grain(self.coarse_grain))

    def _normalised_coarse_grain(
        self, given: Mapping[str, float] | None
    ) -> dict[str, float]:
        """Only the values that differ from makemole's defaults (so equal types compare equal)."""
        given = dict(given or {})
        unknown = set(given) - set(COARSE_GRAIN_KEYS)
        if unknown:
            raise ValueError(
                f"atom type {self.label}: unknown coarse_grain keys {sorted(unknown)}; "
                f"valid: {list(COARSE_GRAIN_KEYS)}"
            )
        defaults = self._coarse_grain_defaults()
        return {
            key: float(given[key])
            for key in COARSE_GRAIN_KEYS
            if key in given and _format.f32(given[key]) != _format.f32(defaults[key])
        }

    def _coarse_grain_defaults(self) -> dict[str, float]:
        return {"nsph": 1.0, "rmin": 0.0, "rmax": 0.5 * self.sigma, "cg": 1.0}

    @property
    def sphere(self) -> tuple[float, float, float, float]:
        """``(nsph, rmin, rmax, cg)`` as written to the ``.mol`` file (defaults filled in)."""
        values = {**self._coarse_grain_defaults(), **self.coarse_grain}
        return tuple(values[key] for key in COARSE_GRAIN_KEYS)  # type: ignore[return-value]

    def __hash__(self) -> int:
        return hash(
            (
                self.label,
                self.element,
                self.epsilon,
                self.sigma,
                self.charge,
                self.mass,
                self.exchangeable,
                tuple(sorted(self.coarse_grain.items())),
            )
        )

    def potential_line(self) -> str:
        """The ``.mol`` ``potential`` line of this type (makemole's own layout)."""
        nsph, rmin, rmax, cg = self.sphere
        return (
            f"potential {self.label:<3} {_format.fixed(self.epsilon, 12, 6)} "
            f"{_format.fixed(self.sigma, 12, 6)} {_format.fixed(self.mass, 12, 6)} "
            f"{_format.fixed(self.charge, 12, 6)} {self.element:<3} "
            f"{_format.exponent(nsph)} {_format.fixed(rmin, 12, 6)} "
            f"{_format.fixed(rmax, 12, 6)} {_format.exponent(cg)}"
        )

    def to_dict(self) -> dict[str, Any]:
        """A JSON-serialisable dict (``AtomType(**d)`` rebuilds it)."""
        values: dict[str, Any] = {
            "label": self.label,
            "element": self.element,
            "epsilon": self.epsilon,
            "sigma": self.sigma,
            "charge": self.charge,
            "mass": self.mass,
            "exchangeable": self.exchangeable,
        }
        if self.coarse_grain:
            values["coarse_grain"] = dict(self.coarse_grain)
        return values


# --------------------------------------------------------------------------- Molecule
@dataclass(frozen=True)
class Atom:
    """One atom of a :class:`Molecule`: its name, atom type and position (Å)."""

    name: str
    type: AtomType
    position: tuple[float, float, float]


class _ClassOrInstanceMethod:
    """``Molecule.atom(atom_type)`` builds a one-atom molecule; ``molecule.atom("H1")``
    returns an atom of an existing molecule."""

    def __init__(
        self, on_class: Callable[..., Any], on_instance: Callable[..., Any]
    ) -> None:
        self._on_class = on_class
        self._on_instance = on_instance
        self.__doc__ = f"{on_class.__doc__}\n\nOn a molecule: {on_instance.__doc__}"

    def __get__(self, instance: Any, owner: type) -> Callable[..., Any]:
        if instance is None:
            return functools.partial(self._on_class, owner)
        return functools.partial(self._on_instance, instance)


def _pair_key(first: int, second: int) -> tuple[int, int]:
    return (first, second) if first < second else (second, first)


def _angle_key(first: int, vertex: int, second: int) -> tuple[int, int, int]:
    return (first, vertex, second) if first < second else (second, vertex, first)


class Molecule:
    """A molecule: named atoms with their types and positions, bonds, angles, rotations.

    Parameters
    ----------
    name:
        The molecule type's name (it becomes the template file ``<name>.mol``).
    atoms:
        ``(atom_name, atom_type, (x, y, z))`` per atom, positions in Å. Atom names
        are unique within the molecule and are only used in Python; the ``.mol``
        file numbers the atoms ``001``, ``002``... in this order.
    bonds:
        Bonded pairs, each atom given by its name or 0-based index:
        ``[("O", "H1"), ("O", "H2")]``.
    angles:
        Angle restraints ``(end, vertex, end)``. ``makemole`` needs one for every
        two atoms bonded to a common atom (and not bonded to each other); ``None``
        (the default) generates exactly those. An explicit list must contain them all.
    rotations:
        Bonds about which EPSR makes internal rotation moves (``rot`` lines).
    bond_lengths, angle_values:
        Target values (Å, degrees) keyed by the bond or angle, where they differ
        from the geometry of ``atoms``. By default they are measured from it.
    vibtemp, angtemp, dihtemp:
        Restraint "temperatures" (stiffness), see the notes.
    ecore, dcore:
        Soft repulsive core (WCA, ε and σ) between atoms of the molecule further
        apart than 1-3; 0 0 switches it off.
    dihedrals:
        Not supported: the Rust ``makemole`` has no ``dihedral`` lines. Giving any
        raises :class:`epsr.EpsrUnsupportedError`.

    Notes
    -----
    **The intramolecular model is EPSR's restraint model, not a force field.**
    Bonds and angles take *target geometry only*, not force constants. Each
    bond is restrained harmonically about its target length, with an energy
    proportional to ``vibtemp`` · RT · √μ / r₀ · (r − r₀)² (μ the reduced mass of
    the two atom types). Angles are converted by ``makemole`` into 1-3 distance
    restraints of the same form, additionally scaled by ``angtemp``. The stiffness
    is therefore set once per molecule by ``vibtemp``, ``angtemp`` and ``dihtemp``
    and scaled by the masses, the bond length and the temperature; per-bond force
    constants (e.g. from Dissolve or OPLS) have no EPSR equivalent.
    """

    def __init__(
        self,
        name: str,
        atoms: Sequence[tuple[str, AtomType, Sequence[float]] | Atom],
        bonds: Iterable[Sequence[AtomRef]] = (),
        angles: Iterable[Sequence[AtomRef]] | None = None,
        rotations: Iterable[Sequence[AtomRef]] = (),
        *,
        bond_lengths: Mapping[Any, float] | None = None,
        angle_values: Mapping[Any, float] | None = None,
        vibtemp: float = 65.0,
        angtemp: float = 1.0,
        dihtemp: float = 10.0,
        ecore: float = 0.0,
        dcore: float = 0.0,
        dihedrals: Iterable[Any] | None = None,
    ) -> None:
        if dihedrals:
            raise EpsrUnsupportedError(
                f"molecule {name}: dihedral restraints are not supported (the Rust "
                "makemole cannot generate DIH groups)"
            )
        self.name = _check_file_name(name, "molecule name")
        self.atoms: tuple[Atom, ...] = tuple(self._atom(entry) for entry in atoms)
        self._check_atoms()
        self._index_by_name = {
            atom.name: index for index, atom in enumerate(self.atoms)
        }
        self.bonds: tuple[tuple[int, int], ...] = self._bonds(bonds)
        self._check_connected()
        required = self.required_angles()
        self.angles: tuple[tuple[int, int, int], ...] = (
            required if angles is None else self._angles(angles, required)
        )
        self.rotations: tuple[tuple[int, int], ...] = self._rotations(rotations)
        self._bond_lengths = self._targets(bond_lengths or {}, 2)
        self._angle_values = self._targets(angle_values or {}, 3)
        self.vibtemp = _positive(vibtemp, "vibtemp")
        self.angtemp = _positive(angtemp, "angtemp")
        self.dihtemp = _positive(dihtemp, "dihtemp")
        self.ecore = _positive(ecore, "ecore", allow_zero=True)
        self.dcore = _positive(dcore, "dcore", allow_zero=True)

    # ---------------------------------------------------------------- construction
    def _atom(self, entry: tuple[str, AtomType, Sequence[float]] | Atom) -> Atom:
        if isinstance(entry, Atom):
            return entry
        try:
            atom_name, atom_type, position = entry
        except (TypeError, ValueError):
            raise TypeError(
                f"molecule {self.name}: atoms are (name, AtomType, (x, y, z)), got {entry!r}"
            ) from None
        if not isinstance(atom_type, AtomType):
            raise TypeError(
                f"molecule {self.name}, atom {atom_name!r}: the type must be an AtomType"
            )
        x, y, z = (float(value) for value in position)
        return Atom(str(atom_name), atom_type, (x, y, z))

    def _check_atoms(self) -> None:
        if not self.atoms:
            raise ValueError(f"molecule {self.name}: no atoms")
        if len(self.atoms) > MAX_ATOMS_PER_MOLECULE:
            raise ValueError(
                f"molecule {self.name}: {len(self.atoms)} atoms; makemole numbers atoms "
                f"with three digits (at most {MAX_ATOMS_PER_MOLECULE})"
            )
        names: set[str] = set()
        types: dict[str, AtomType] = {}
        for index, atom in enumerate(self.atoms):
            if not atom.name or any(character.isspace() for character in atom.name):
                raise ValueError(f"molecule {self.name}: bad atom name {atom.name!r}")
            if atom.name in names:
                raise ValueError(
                    f"molecule {self.name}: atom name {atom.name!r} is used twice"
                )
            names.add(atom.name)
            other = types.setdefault(atom.type.label, atom.type)
            if other != atom.type:
                raise ValueError(
                    f"molecule {self.name}: two different atom types are labelled "
                    f"{atom.type.label!r}"
                )
            for other_atom in self.atoms[:index]:
                if math.dist(atom.position, other_atom.position) < MIN_ATOM_SEPARATION:
                    raise ValueError(
                        f"molecule {self.name}: atoms {other_atom.name} and {atom.name} "
                        f"are closer than {MIN_ATOM_SEPARATION} Å"
                    )

    def _index(self, reference: AtomRef) -> int:
        """0-based index of an atom given by name or index."""
        if isinstance(reference, bool):
            raise TypeError(f"molecule {self.name}: bad atom reference {reference!r}")
        if isinstance(reference, int):
            if not 0 <= reference < len(self.atoms):
                raise IndexError(
                    f"molecule {self.name}: atom index {reference} out of range "
                    f"(0..{len(self.atoms) - 1})"
                )
            return reference
        try:
            return self._index_by_name[reference]
        except KeyError:
            raise KeyError(
                f"molecule {self.name}: no atom named {reference!r}; atoms are "
                f"{list(self._index_by_name)}"
            ) from None

    def _bonds(self, bonds: Iterable[Sequence[AtomRef]]) -> tuple[tuple[int, int], ...]:
        result: list[tuple[int, int]] = []
        seen: set[tuple[int, int]] = set()
        for bond in bonds:
            first, second = (self._index(reference) for reference in bond)
            if first == second:
                raise ValueError(f"molecule {self.name}: atom {first} bonded to itself")
            if _pair_key(first, second) in seen:
                raise ValueError(
                    f"molecule {self.name}: bond {self._label(first)}-{self._label(second)} "
                    "given twice"
                )
            seen.add(_pair_key(first, second))
            result.append((first, second))
        return tuple(result)

    def _check_connected(self) -> None:
        if len(self.atoms) == 1:
            return
        if not self.bonds:
            raise ValueError(
                f"molecule {self.name}: {len(self.atoms)} atoms but no bonds"
            )
        neighbours = self.neighbours()
        seen = {0}
        stack = [0]
        while stack:
            for other in neighbours[stack.pop()]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        if len(seen) != len(self.atoms):
            loose = [
                self.atoms[index].name
                for index in range(len(self.atoms))
                if index not in seen
            ]
            raise ValueError(
                f"molecule {self.name}: the bonds do not connect atoms {loose} to the rest"
            )

    def _angles(
        self,
        angles: Iterable[Sequence[AtomRef]],
        required: tuple[tuple[int, int, int], ...],
    ) -> tuple[tuple[int, int, int], ...]:
        bonded = {_pair_key(*bond) for bond in self.bonds}
        result: list[tuple[int, int, int]] = []
        seen: set[tuple[int, int, int]] = set()
        for angle in angles:
            first, vertex, second = (self._index(reference) for reference in angle)
            if (
                first == second
                or _pair_key(first, vertex) not in bonded
                or _pair_key(vertex, second) not in bonded
            ):
                raise ValueError(
                    f"molecule {self.name}: angle {self._triple(first, vertex, second)} "
                    "needs both ends bonded to the middle atom"
                )
            key = _angle_key(first, vertex, second)
            if key not in seen:
                seen.add(key)
                result.append((first, vertex, second))
        missing = [triple for triple in required if _angle_key(*triple) not in seen]
        if missing:
            raise ValueError(
                f"molecule {self.name}: makemole needs a restraint for every angle "
                f"between two bonds; missing {[self._triple(*t) for t in missing]} "
                "(pass angles=None to generate them all from the geometry)"
            )
        return tuple(result)

    def _rotations(
        self, rotations: Iterable[Sequence[AtomRef]]
    ) -> tuple[tuple[int, int], ...]:
        bonded = {_pair_key(*bond) for bond in self.bonds}
        result = []
        for rotation in rotations:
            first, second = (self._index(reference) for reference in rotation)
            if _pair_key(first, second) not in bonded:
                raise ValueError(
                    f"molecule {self.name}: rotation {self._label(first)}-"
                    f"{self._label(second)} is not about a bond"
                )
            result.append((first, second))
        return tuple(result)

    def _targets(
        self, given: Mapping[Any, float], size: int
    ) -> dict[tuple[int, ...], float]:
        """Target bond lengths (size 2) or angles (size 3) keyed by sorted indices."""
        known = (
            {_pair_key(*bond) for bond in self.bonds}
            if size == 2
            else {_angle_key(*angle) for angle in self.angles}
        )
        what = "bond" if size == 2 else "angle"
        targets: dict[tuple[int, ...], float] = {}
        for references, value in given.items():
            indices = [self._index(reference) for reference in references]
            if len(indices) != size:
                raise ValueError(f"molecule {self.name}: bad {what} key {references!r}")
            key = _pair_key(*indices) if size == 2 else _angle_key(*indices)
            if key not in known:
                raise ValueError(
                    f"molecule {self.name}: {references!r} is not one of the {what}s"
                )
            value = _positive(value, f"molecule {self.name}: {what} {references!r}")
            if size == 3 and value > 180:
                raise ValueError(f"molecule {self.name}: angle {value} > 180 degrees")
            targets[key] = value
        return targets

    # ---------------------------------------------------------------- alternative constructors
    @classmethod
    def from_geometry(
        cls,
        name: str,
        centre: tuple[str, AtomType],
        arms: Sequence[tuple[str, AtomType, float]],
        angle: float | None = None,
        **options: Any,
    ) -> Molecule:
        """A molecule of one central atom bonded to 1-4 others (water, NH₃, CH₄).

        ``arms`` gives ``(name, type, bond length)`` of each outer atom and
        ``angle`` the angle between any two arms, in degrees (omit it for one arm;
        four arms are tetrahedral). The centre sits at the origin and the first arm
        on +x; with two arms the second lies in the xy plane. The restraint targets
        are exactly the given lengths and angle. ``options`` are passed to
        :class:`Molecule` (``vibtemp``, ``rotations``...).
        """
        centre_name, centre_type = centre
        directions = _arm_directions(len(arms), angle, name)
        atoms: list[tuple[str, AtomType, Sequence[float]]] = [
            (centre_name, centre_type, (0.0, 0.0, 0.0))
        ]
        for (arm_name, arm_type, length), direction in zip(
            arms, directions, strict=True
        ):
            length = _positive(length, f"molecule {name}: bond length of {arm_name}")
            atoms.append((arm_name, arm_type, tuple(length * c for c in direction)))
        bonds = [(centre_name, arm[0]) for arm in arms]
        tetrahedral = len(arms) == 4 and angle is None
        arm_angle = _TETRAHEDRAL_ANGLE if tetrahedral else angle
        angle_values = {
            (arms[i][0], centre_name, arms[j][0]): arm_angle
            for i in range(len(arms))
            for j in range(i + 1, len(arms))
        }
        options.setdefault(
            "bond_lengths", {(centre_name, arm[0]): arm[2] for arm in arms}
        )
        options.setdefault("angle_values", angle_values if len(arms) > 1 else {})
        return cls(name, atoms, bonds, **options)

    @classmethod
    def _single_atom(
        cls, atom_type: AtomType, name: str | None = None, **options: Any
    ) -> Molecule:
        """A one-atom "molecule" of ``atom_type``, e.g. for ions or for Si and O in silica.

        The molecule is named after the type label unless ``name`` is given; its
        atom is named after the type label too.
        """
        return cls(
            name or atom_type.label,
            [(atom_type.label, atom_type, (0.0, 0.0, 0.0))],
            **options,
        )

    def _atom_by_reference(self, reference: AtomRef) -> Atom:
        """the atom with this name (or 0-based index)."""
        return self.atoms[self._index(reference)]

    atom = _ClassOrInstanceMethod(_single_atom.__func__, _atom_by_reference)  # type: ignore[attr-defined]

    @classmethod
    def from_xyz(
        cls,
        path: str | PathLike[str],
        types: Mapping[str, AtomType],
        bonds: str | Iterable[Sequence[AtomRef]] = "auto",
        *,
        name: str | None = None,
        tolerance: float = BOND_TOLERANCE,
        **options: Any,
    ) -> Molecule:
        """A molecule from an XYZ file (atom count, comment line, ``symbol x y z`` lines).

        ``types`` maps the symbol of the first column to an :class:`AtomType`.
        Atoms are named by that symbol and a running number (``O1``, ``H1``,
        ``H2``). ``bonds="auto"`` bonds atoms closer than the sum of their
        covalent radii plus ``tolerance`` (Å); otherwise give the bonds. The
        molecule is named after the file unless ``name`` is given.
        """
        path = Path(path)
        lines = path.read_text(encoding="latin-1").splitlines()
        try:
            n_atoms = int(lines[0].split()[0])
            rows = [line.split() for line in lines[2 : 2 + n_atoms]]
            entries = [
                (row[0], tuple(float(value) for value in row[1:4])) for row in rows
            ]
        except (IndexError, ValueError):
            raise ValueError(f"{path}: not an XYZ file") from None
        if len(entries) != n_atoms:
            raise ValueError(f"{path}: {n_atoms} atoms announced, {len(entries)} found")
        counts: dict[str, int] = {}
        atoms = []
        for symbol, position in entries:
            if symbol not in types:
                raise KeyError(
                    f"{path}: no atom type given for {symbol!r} (types={list(types)})"
                )
            counts[symbol] = counts.get(symbol, 0) + 1
            atoms.append((f"{symbol}{counts[symbol]}", types[symbol], position))
        if isinstance(bonds, str):
            if bonds != "auto":
                raise ValueError(
                    f"bonds must be 'auto' or a list of pairs, got {bonds!r}"
                )
            bonds = _detect_bonds(atoms, tolerance, path)
        return cls(name or path.stem, atoms, bonds, **options)

    @classmethod
    def from_mol(
        cls,
        path: str | PathLike[str],
        types: Mapping[str, AtomType] | None = None,
        *,
        name: str | None = None,
    ) -> Molecule:
        """A molecule from an EPSR ``.mol`` (gmol) file, as read by ``makemole``.

        Atom types come from the ``potential`` lines unless ``types`` (keyed by
        label) replaces them. Bond lengths, angles, rotations and ``vibtemp``,
        ``angtemp``, ``dihtemp``, ``ecoredcore`` are taken from the file; the
        ``temperature`` line is not (the :class:`System` sets it). Atoms are named
        by element and a running number. ``.mol`` lines this API cannot represent
        (``dihedral``, ``extra``, ``wholerot``, ``introt``, ``qradius``, a
        ``pseudofac`` other than 1) raise :class:`epsr.EpsrUnsupportedError`.
        """
        path = Path(path)
        return _MolFile(path).molecule(name or path.stem, types or {})

    @classmethod
    def from_ato(
        cls, path: str | PathLike[str], *, name: str | None = None
    ) -> Molecule:
        """A molecule from a one-molecule template ``.ato`` file (e.g. ``EPSR-test/DMSO.ato``).

        Atom types, positions, bonds, rotations and the restraint temperatures are
        read from the file. **Approximate:** the ``.ato`` file stores restraint
        distances with three significant digits, so bond lengths and the angles
        rebuilt from the 1-3 distances are rounded; prefer :meth:`from_mol` when
        the ``.mol`` file exists.
        """
        path = Path(path)
        return _AtoFile(path).molecule(name or path.stem)

    # ---------------------------------------------------------------- properties
    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def atom_names(self) -> list[str]:
        return [atom.name for atom in self.atoms]

    @property
    def atom_types(self) -> list[AtomType]:
        """The atom types of the molecule, in order of first appearance."""
        return list({atom.type.label: atom.type for atom in self.atoms}.values())

    def _label(self, index: int) -> str:
        return self.atoms[index].name

    def _triple(self, first: int, vertex: int, second: int) -> str:
        return f"{self._label(first)}-{self._label(vertex)}-{self._label(second)}"

    def neighbours(self) -> list[list[int]]:
        """Bonded neighbours of every atom (0-based), in bond order."""
        result: list[list[int]] = [[] for _ in self.atoms]
        for first, second in self.bonds:
            result[first].append(second)
            result[second].append(first)
        return result

    def required_angles(self) -> tuple[tuple[int, int, int], ...]:
        """The angles ``makemole`` needs: two atoms bonded to a common atom, not to each other.

        Ordered by the middle atom, then by its bond order (0-based indices).
        """
        bonded = {_pair_key(*bond) for bond in self.bonds}
        triples = []
        for vertex, around in enumerate(self.neighbours()):
            for position, first in enumerate(around):
                for second in around[position + 1 :]:
                    if _pair_key(first, second) not in bonded:
                        triples.append((first, vertex, second))
        return tuple(triples)

    def bond_length(self, first: AtomRef, second: AtomRef) -> float:
        """Target length (Å) of a bond: as given, or measured from the positions."""
        i, j = self._index(first), self._index(second)
        key = _pair_key(i, j)
        if key in self._bond_lengths:
            return self._bond_lengths[key]
        return math.dist(self.atoms[i].position, self.atoms[j].position)

    def angle_value(self, first: AtomRef, vertex: AtomRef, second: AtomRef) -> float:
        """Target angle (degrees): as given, or measured from the positions."""
        i, v, k = (self._index(reference) for reference in (first, vertex, second))
        key = _angle_key(i, v, k)
        if key in self._angle_values:
            return self._angle_values[key]
        a, b, c = (self.atoms[index].position for index in (i, v, k))
        u = [a[n] - b[n] for n in range(3)]
        w = [c[n] - b[n] for n in range(3)]
        cosine = sum(u[n] * w[n] for n in range(3)) / (math.hypot(*u) * math.hypot(*w))
        return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

    def __repr__(self) -> str:
        return (
            f"<epsr.Molecule {self.name}: {self.n_atoms} atoms "
            f"({', '.join(self.atom_names)}), {len(self.bonds)} bonds, "
            f"{len(self.angles)} angles, {len(self.rotations)} rotations>"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Molecule):
            return NotImplemented
        return self.to_dict() == other.to_dict() and self.atom_types == other.atom_types

    __hash__ = None  # type: ignore[assignment]

    # ---------------------------------------------------------------- the .mol file
    def render_mol(self, temperature: float) -> str:
        """The ``.mol`` (gmol) text ``makemole`` reads, in the layout it writes itself.

        ``temperature`` (K) is that of the :class:`System`. Coordinates and values
        are rounded to single precision as EPSR stores them.
        """
        temperature = _positive(temperature, "temperature")
        neighbours = self.neighbours()
        lines = [f" .gmol{0:12d}"]
        for index, atom in enumerate(self.atoms):
            around = neighbours[index]
            # makemole's own layout: a list-directed record after the 3-digit label
            items = "".join(
                f" {_format.list_real_field(value)}" for value in atom.position
            )
            integers = "".join(
                f" {value:11d}" for value in [len(around), *(n + 1 for n in around)]
            )
            lines.append(
                f" atom{index + 1:13d}  {index + 1:03d} {items}{integers}".rstrip()
            )
        for first, second in self.bonds:
            length = _format.fixed(self.bond_length(first, second), 10, 5)
            lines.append(f"bond {first + 1:03d} {second + 1:03d} {length}")
        for first, vertex, second in self.angles:
            value = _format.fixed(self.angle_value(first, vertex, second), 10, 5)
            lines.append(
                f"angle {first + 1:03d} {vertex + 1:03d} {second + 1:03d} {value}"
            )
        for first, second in self.rotations:
            lines.append(f"rot {first + 1:3d} {second + 1:3d}")
        for index, atom in enumerate(self.atoms):
            lines.append(f"changelabel {index + 1:03d} {atom.type.label:<3}")
        lines.append(f"pseudofac {_format.fixed(1.0, 8, 2)}")
        lines += [atom_type.potential_line() for atom_type in self.atom_types]
        lines.append(f"temperature {_format.exponent(temperature)}")
        lines.append(f"vibtemp {_format.exponent(self.vibtemp)}")
        lines.append(f"angtemp {_format.exponent(self.angtemp)}")
        lines.append(f"dihtemp {_format.exponent(self.dihtemp)}")
        lines.append(
            f"ecoredcore {_format.fixed(self.ecore, 10, 5)} {_format.fixed(self.dcore, 10, 5)}"
        )
        lines.append(f"density {_format.exponent(TEMPLATE_DENSITY)}")
        return "".join(f"{line}\n" for line in lines)

    # ---------------------------------------------------------------- serialisation
    def to_dict(self) -> dict[str, Any]:
        """A JSON-serialisable dict; atom types by label (see :meth:`System.to_dict`)."""
        return {
            "name": self.name,
            "atoms": [
                [atom.name, atom.type.label, list(atom.position)] for atom in self.atoms
            ],
            "bonds": [list(bond) for bond in self.bonds],
            "angles": [list(angle) for angle in self.angles],
            "rotations": [list(rotation) for rotation in self.rotations],
            "bond_lengths": [
                [*key, value] for key, value in sorted(self._bond_lengths.items())
            ],
            "angle_values": [
                [*key, value] for key, value in sorted(self._angle_values.items())
            ],
            "vibtemp": self.vibtemp,
            "angtemp": self.angtemp,
            "dihtemp": self.dihtemp,
            "ecore": self.ecore,
            "dcore": self.dcore,
        }

    @classmethod
    def from_dict(
        cls, data: Mapping[str, Any], types: Mapping[str, AtomType]
    ) -> Molecule:
        """Rebuild from :meth:`to_dict`, with the atom types keyed by label."""
        return cls(
            data["name"],
            [
                (atom_name, types[label], position)
                for atom_name, label, position in data["atoms"]
            ],
            [tuple(bond) for bond in data["bonds"]],
            [tuple(angle) for angle in data["angles"]],
            [tuple(rotation) for rotation in data["rotations"]],
            bond_lengths={tuple(entry[:2]): entry[2] for entry in data["bond_lengths"]},
            angle_values={tuple(entry[:3]): entry[3] for entry in data["angle_values"]},
            vibtemp=data["vibtemp"],
            angtemp=data["angtemp"],
            dihtemp=data["dihtemp"],
            ecore=data["ecore"],
            dcore=data["dcore"],
        )


_TETRAHEDRAL_ANGLE = math.degrees(math.acos(-1.0 / 3.0))


def _arm_directions(
    n_arms: int, angle: float | None, name: str
) -> list[tuple[float, float, float]]:
    """Unit vectors of ``n_arms`` bonds from a central atom, ``angle`` apart."""
    if n_arms == 1:
        if angle is not None:
            raise ValueError(f"molecule {name}: one arm has no angle")
        return [(1.0, 0.0, 0.0)]
    if n_arms == 4:
        if angle is not None and abs(angle - _TETRAHEDRAL_ANGLE) > 1e-3:
            raise ValueError(
                f"molecule {name}: four equal arms are tetrahedral "
                f"({_TETRAHEDRAL_ANGLE:.4f} degrees); give the atoms explicitly otherwise"
            )
        s = 1.0 / math.sqrt(3.0)
        return [(s, s, s), (s, -s, -s), (-s, s, -s), (-s, -s, s)]
    if n_arms not in (2, 3):
        raise ValueError(
            f"molecule {name}: from_geometry takes 1 to 4 arms, got {n_arms}"
        )
    if angle is None or not 0 < angle < (180 if n_arms == 2 else 120):
        raise ValueError(f"molecule {name}: give the angle between the arms (degrees)")
    theta = math.radians(angle)
    if n_arms == 2:
        return [(1.0, 0.0, 0.0), (math.cos(theta), math.sin(theta), 0.0)]
    # three arms on a cone about z, 120 degrees apart in azimuth
    cos_polar = math.sqrt((math.cos(theta) + 0.5) / 1.5)
    sin_polar = math.sqrt(1.0 - cos_polar**2)
    return [
        (
            sin_polar * math.cos(math.radians(120 * arm)),
            sin_polar * math.sin(math.radians(120 * arm)),
            cos_polar,
        )
        for arm in range(3)
    ]


def _detect_bonds(
    atoms: Sequence[tuple[str, AtomType, Sequence[float]]], tolerance: float, path: Path
) -> list[tuple[int, int]]:
    radii = []
    for _, atom_type, _ in atoms:
        if atom_type.element not in COVALENT_RADII:
            raise ValueError(
                f"{path}: no covalent radius for {atom_type.element}; give the bonds explicitly"
            )
        radii.append(COVALENT_RADII[atom_type.element])
    return [
        (i, j)
        for i in range(len(atoms))
        for j in range(i + 1, len(atoms))
        if math.dist(atoms[i][2], atoms[j][2]) < radii[i] + radii[j] + tolerance
    ]


def _running_names(elements: Sequence[str]) -> list[str]:
    counts: dict[str, int] = {}
    names = []
    for element in elements:
        counts[element] = counts.get(element, 0) + 1
        names.append(f"{element}{counts[element]}")
    return names


class _MolFile:
    """The parts of a gmol ``.mol`` file that :meth:`Molecule.from_mol` uses."""

    _UNSUPPORTED = ("dihedral", "extra", "wholerot", "introt", "qradius", "unitcell")

    def __init__(self, path: Path) -> None:
        self.path = path
        self.labels: list[str] = []  # 3-digit atom labels, "001"...
        self.positions: list[tuple[float, float, float]] = []
        self.neighbours: list[list[int]] = []  # 1-based, as written
        self.bonds: dict[frozenset[str], float] = {}
        self.angles: list[tuple[tuple[str, str, str], float]] = []
        self.rotations: list[tuple[int, int]] = []
        self.type_of: dict[str, str] = {}
        self.potentials: dict[str, AtomType] = {}
        self.settings: dict[str, float] = {}
        lines = path.read_text(encoding="latin-1").splitlines()
        if not lines or not lines[0].split() or lines[0].split()[0] != ".gmol":
            raise ValueError(f"{path}: not a gmol .mol file (first line must be .gmol)")
        for line in lines[1:]:
            words = line.split()
            if words:
                self._read(words)

    def _read(self, words: list[str]) -> None:
        keyword = words[0]
        if keyword in self._UNSUPPORTED:
            raise EpsrUnsupportedError(
                f"{self.path}: '{keyword}' lines are not supported"
            )
        if keyword == "atom":
            self.labels.append(words[2])
            self.positions.append((float(words[3]), float(words[4]), float(words[5])))
            self.neighbours.append([int(word) for word in words[7 : 7 + int(words[6])]])
        elif keyword == "bond":
            self.bonds[frozenset(words[1:3])] = float(words[3])
        elif keyword == "angle":
            self.angles.append(((words[1], words[2], words[3]), float(words[4])))
        elif keyword == "rot":
            self.rotations.append((int(words[1]), int(words[2])))
        elif keyword == "changelabel":
            self.type_of[words[1]] = words[2]
        elif keyword == "potential":
            label, element = words[1], words[6]
            sphere = (
                dict(
                    zip(
                        COARSE_GRAIN_KEYS,
                        (float(word) for word in words[7:11]),
                        strict=True,
                    )
                )
                if len(words) > 10
                else {}
            )
            self.potentials[label] = AtomType(
                label,
                element,
                epsilon=float(words[2]),
                sigma=float(words[3]),
                charge=float(words[5]),
                mass=float(words[4]),
                coarse_grain=sphere,
            )
        elif keyword == "pseudofac":
            if float(words[1]) != 1.0:
                raise EpsrUnsupportedError(f"{self.path}: pseudofac other than 1")
        elif keyword in ("vibtemp", "angtemp", "dihtemp"):
            self.settings[keyword] = float(words[1])
        elif keyword == "ecoredcore":
            self.settings["ecore"] = float(words[1])
            self.settings["dcore"] = float(words[2]) if len(words) > 2 else 0.0

    def molecule(self, name: str, types: Mapping[str, AtomType]) -> Molecule:
        atom_types = []
        for label in self.labels:
            type_label = self.type_of.get(label, label)
            atom_type = types.get(type_label) or self.potentials.get(type_label)
            if atom_type is None:
                raise ValueError(
                    f"{self.path}: no potential line for atom type {type_label!r}"
                )
            atom_types.append(atom_type)
        names = _running_names([atom_type.element for atom_type in atom_types])
        index_of_label = {label: index for index, label in enumerate(self.labels)}
        bonds: list[tuple[int, int]] = []
        seen: set[tuple[int, int]] = set()
        for index, around in enumerate(self.neighbours):
            for other in around:
                key = _pair_key(index, other - 1)
                if key not in seen:
                    seen.add(key)
                    bonds.append((index, other - 1))
        bond_lengths = {}
        for first, second in bonds:
            labels = frozenset((self.labels[first], self.labels[second]))
            if labels in self.bonds:
                bond_lengths[(first, second)] = self.bonds[labels]
        angle_values = {
            tuple(index_of_label[label] for label in labels): value
            for labels, value in self.angles
        }
        return Molecule(
            name,
            [
                (n, t, p)
                for n, t, p in zip(names, atom_types, self.positions, strict=True)
            ],
            bonds,
            angles=list(angle_values) or None,
            rotations=[(first - 1, second - 1) for first, second in self.rotations],
            bond_lengths=bond_lengths,
            angle_values=angle_values,
            vibtemp=self.settings.get("vibtemp", 65.0),
            angtemp=self.settings.get("angtemp", 1.0),
            dihtemp=self.settings.get("dihtemp", 10.0),
            ecore=self.settings.get("ecore", 0.0),
            dcore=self.settings.get("dcore", 0.0),
        )


class _AtoFile:
    """The first molecule of a template ``.ato`` file, read token by token (``moloutc`` layout)."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.tokens = path.read_text(encoding="latin-1").split()
        self.position = 0

    def _next(self) -> str:
        if self.position >= len(self.tokens):
            raise ValueError(f"{self.path}: unexpected end of file")
        token = self.tokens[self.position]
        self.position += 1
        return token

    def _int(self) -> int:
        return int(self._next())

    def _real(self) -> float:
        return float(self._next())

    def molecule(self, name: str) -> Molecule:
        path = self.path
        header = self.tokens[:3]
        if len(header) < 3 or int(header[0]) != 1:
            raise ValueError(f"{path}: a one-molecule template .ato file is needed")
        self.position = 3  # nmol, box side, temperature (cubic box)
        mc = [
            self._real() for _ in range(10)
        ]  # tol, step sizes, vibtemp, angtemp, dihtemp...
        n_atoms = self._int()
        self.position += 6  # centre and orientation of the molecule
        self._next()  # tethered flag
        self.position += 3  # rejected steps, scale factor, molecule number
        labels, positions, restraints = [], [], []
        for _ in range(n_atoms):
            labels.append(self._next())
            self.position += 2  # atom number, image atom
            positions.append((self._real(), self._real(), self._real()))
            restraints.append(
                [(self._int() - 1, self._real()) for _ in range(self._int())]
            )
        rotations = []
        for _ in range(self._int()):
            kind = self._next()
            if kind != "ROT":
                raise EpsrUnsupportedError(
                    f"{path}: {kind} rotational groups are not supported"
                )
            rotations.append((self._int() - 1, self._int() - 1))
            n_moving = self._int()
            self.position += n_moving  # the atoms that move with the rotation
        types: dict[str, AtomType] = {}
        while self.position < len(self.tokens) and self.tokens[self.position] in labels:
            label, element, exchange = self._next(), self._next(), self._int()
            epsilon, sigma, mass, charge, qradius = (self._real() for _ in range(5))
            nsph, rmin, rmax, cg = (self._real() for _ in range(4))
            if qradius != 0.0:
                raise EpsrUnsupportedError(
                    f"{path}: charge radius (qradius) of {label}"
                )
            types[label] = AtomType(
                label,
                element,
                epsilon,
                sigma,
                charge,
                mass,
                exchangeable=bool(exchange),
                coarse_grain=dict(
                    zip(COARSE_GRAIN_KEYS, (nsph, rmin, rmax, cg), strict=True)
                ),
            )
        ecore, dcore = self._real(), self._real()
        atom_types = [types[label] for label in labels]
        names = _running_names([atom_type.element for atom_type in atom_types])
        bonds, bond_lengths, distances = [], {}, {}
        for index, entries in enumerate(restraints):
            for other, value in entries:
                key = _pair_key(index, other)
                if value > 0 and key not in bond_lengths:
                    bonds.append(key)
                    bond_lengths[key] = value
                elif value < 0:
                    distances[key] = -value
        molecule = Molecule(
            name,
            list(zip(names, atom_types, positions, strict=True)),
            bonds,
            rotations=rotations,
            bond_lengths=bond_lengths,
            vibtemp=mc[5],
            angtemp=mc[6],
            dihtemp=mc[7],
            ecore=ecore,
            dcore=dcore,
        )
        angle_values = {}
        for first, vertex, second in molecule.angles:
            distance = distances.get(_pair_key(first, second))
            if distance is None:
                continue
            a = molecule.bond_length(first, vertex)
            b = molecule.bond_length(vertex, second)
            cosine = (a * a + b * b - distance * distance) / (2 * a * b)
            angle_values[(first, vertex, second)] = math.degrees(
                math.acos(max(-1.0, min(1.0, cosine)))
            )
        molecule._angle_values = molecule._targets(angle_values, 3)
        return molecule


# --------------------------------------------------------------------------- pairs
@dataclass(frozen=True)
class PairParameters:
    """The interactions of one pair of atom types.

    ``epsilon`` and ``sigma`` are the Lorentz–Berthelot mixed values,
    ``charge_product`` the Coulomb prefactor 1390 q_i q_j (kJ mol⁻¹ Å) and
    ``terms`` the interactions that are actually present: ``"LJ + Coulomb"``,
    ``"LJ only"``, ``"Coulomb only"`` or ``"none"``, followed by any extra
    per-pair terms of the ``.pcof`` file.

    Rows read from a running engine (``result.potential.pair_params()``) also
    have the per-pair terms of the ``.pcof`` file: ``coulfac`` (scale factor on
    the Coulomb term), the extra repulsive core (``rminex``, ``pcoefex``,
    ``coordex``; ``rmaxex``, ``pcoefexmax``) and additional Gaussian or
    modified-Morse terms (``gaussians``, ``(position, width, height)``). Rows of
    :meth:`System.pair_parameters` leave them ``None``.
    """

    types: tuple[str, str]
    epsilon: float
    sigma: float
    charge_product: float
    terms: str
    coulfac: float | None = None
    rminex: float | None = None
    pcoefex: float | None = None
    coordex: float | None = None
    rmaxex: float | None = None
    pcoefexmax: float | None = None
    gaussians: tuple[tuple[float, float, float], ...] = ()

    @property
    def label(self) -> str:
        return f"{self.types[0]}-{self.types[1]}"

    @property
    def extras(self) -> str:
        """The non-zero ``.pcof`` per-pair terms as text (``""`` if there are none)."""
        values = [
            f"{name} {value:g}"
            for name in ("rminex", "pcoefex", "coordex", "rmaxex", "pcoefexmax")
            if (value := getattr(self, name))
        ]
        if self.gaussians:
            values.append(f"{len(self.gaussians)} additional term(s)")
        return ", ".join(values)

    @staticmethod
    def terms_of(
        epsilon: float,
        charge_product: float,
        *,
        repulsive_core: bool = False,
        n_additional: int = 0,
    ) -> str:
        """``terms`` for a mixed epsilon and (scaled) Coulomb prefactor."""
        lennard_jones, coulomb = epsilon > 0, charge_product != 0
        if lennard_jones and coulomb:
            terms = "LJ + Coulomb"
        elif lennard_jones:
            terms = "LJ only"
        elif coulomb:
            terms = "Coulomb only"
        else:
            terms = "none"
        if repulsive_core:
            terms += " + repulsive core"
        if n_additional:
            terms += f" + {n_additional} additional"
        return terms


class PairTable(Sequence[PairParameters]):
    """Rows of :class:`PairParameters` in EPSR's pair order, printable as a table."""

    COLUMNS: tuple[tuple[str, str], ...] = (
        ("pair", "label"),
        ("eps_ij", "epsilon"),
        ("sigma_ij", "sigma"),
        ("1390*qi*qj", "charge_product"),
        ("terms", "terms"),
    )
    #: Extra columns shown when the rows come from the engine.
    ENGINE_COLUMNS: tuple[tuple[str, str], ...] = (
        ("coulfac", "coulfac"),
        (".pcof extras", "extras"),
    )

    def __init__(self, rows: Iterable[PairParameters]) -> None:
        self._rows = list(rows)

    @property
    def columns(self) -> tuple[tuple[str, str], ...]:
        """``(heading, attribute)`` of the printed columns."""
        from_engine = any(row.coulfac is not None for row in self._rows)
        return self.COLUMNS + (self.ENGINE_COLUMNS if from_engine else ())

    def __getitem__(self, index: Any) -> Any:
        if isinstance(index, (str, tuple)):
            return self.pair(index)
        return self._rows[index]

    def __len__(self) -> int:
        return len(self._rows)

    def pair(self, pair: str | tuple[str, str]) -> PairParameters:
        """The row of a pair, as ``"OW-HW"`` or ``("HW", "OW")`` (either order)."""
        first, second = pair.split("-") if isinstance(pair, str) else pair
        for row in self._rows:
            if row.types in ((first, second), (second, first)):
                return row
        raise KeyError(
            f"no pair {first}-{second}; pairs are {[row.label for row in self._rows]}"
        )

    @staticmethod
    def _cell(value: Any, attribute: str) -> str:
        if attribute == "charge_product":
            return f"{value:+.1f}"
        if isinstance(value, float):
            return f"{value:.4g}"
        return str(value)

    def __str__(self) -> str:
        columns = self.columns
        table = [[heading for heading, _ in columns]]
        for row in self._rows:
            table.append(
                [
                    self._cell(getattr(row, attribute), attribute)
                    for _, attribute in columns
                ]
            )
        widths = [
            max(len(line[column]) for line in table) for column in range(len(columns))
        ]
        return "\n".join(
            "  ".join(
                cell.ljust(width) for cell, width in zip(line, widths, strict=True)
            ).rstrip()
            for line in table
        )

    def __repr__(self) -> str:
        return f"<epsr.PairTable {len(self)} pairs>\n{self}"

    def to_dataframe(self) -> Any:
        """The table as a ``pandas.DataFrame`` (requires the ``pandas`` extra)."""
        from .io import import_optional

        pandas = import_optional("pandas", extra="pandas")
        return pandas.DataFrame(
            [
                {attribute: getattr(row, attribute) for _, attribute in self.columns}
                for row in self._rows
            ]
        )


# --------------------------------------------------------------------------- System
class System:
    """The simulation box: molecules, how many of each, density and temperature.

    Parameters
    ----------
    name:
        The box name: the box is ``<name>.ato`` and every EPSR output file starts
        with ``<name>.EPSR``.
    temperature:
        Temperature in K.
    density:
        Atomic number density in atoms/Å³, EPSR's own unit (water: 0.1). The box
        is cubic; its side follows from the number of atoms.
    box:
        Reserved for non-cubic boxes, which are not supported yet.

    Use :meth:`add` to put molecules in, in order: **the order of the ``add``
    calls (and of the atoms within each molecule) fixes the atom-type order**,
    and with it the order of every pair-resolved result (g(r), potentials...).
    """

    def __init__(
        self, name: str, temperature: float, density: float, *, box: Any = None
    ) -> None:
        if box is not None:
            raise EpsrUnsupportedError(
                "non-cubic boxes are not supported yet; the box is cubic, set by the density"
            )
        self.name = _check_file_name(name, "system name")
        self.temperature = _positive(temperature, "temperature")
        self.density = _positive(density, "density (atoms/Å³)")
        self.components: list[tuple[Molecule, int]] = []

    def add(self, molecule: Molecule, count: int) -> System:
        """Add ``count`` copies of ``molecule``; returns the system, so calls can chain.

        Molecules added first come first in the box, and so do their atom types.
        """
        if not isinstance(molecule, Molecule):
            raise TypeError("System.add takes a Molecule")
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ValueError(
                f"count of {molecule.name} must be an int >= 1, got {count!r}"
            )
        if molecule.name == self.name:
            raise ValueError(
                f"molecule {molecule.name!r} has the system's name; the box file would "
                "overwrite the molecule template"
            )
        if any(existing.name == molecule.name for existing, _ in self.components):
            raise ValueError(
                f"a molecule named {molecule.name!r} is already in the system"
            )
        n_molecules = self.n_molecules + count
        if n_molecules > MAX_MOLECULES:
            raise ValueError(
                f"{n_molecules} molecules; the engine allows at most {MAX_MOLECULES}"
            )
        types = {atom_type.label: atom_type for atom_type in self.atom_types}
        for atom_type in molecule.atom_types:
            other = types.setdefault(atom_type.label, atom_type)
            if other != atom_type:
                raise ValueError(
                    f"atom type {atom_type.label!r} of {molecule.name} differs from the "
                    "type with that label already in the system"
                )
        if len(types) > MAX_ATOM_TYPES:
            raise ValueError(
                f"{len(types)} atom types; the engine allows at most {MAX_ATOM_TYPES}"
            )
        self.components.append((molecule, count))
        return self

    @property
    def molecules(self) -> list[Molecule]:
        return [molecule for molecule, _ in self.components]

    @property
    def atom_types(self) -> list[AtomType]:
        """The atom types in EPSR's order: by molecule (``add`` order), then by first appearance."""
        types: dict[str, AtomType] = {}
        for molecule, _ in self.components:
            for atom_type in molecule.atom_types:
                types.setdefault(atom_type.label, atom_type)
        return list(types.values())

    def atom_type(self, label: str) -> AtomType:
        for atom_type in self.atom_types:
            if atom_type.label == label:
                return atom_type
        raise KeyError(
            f"no atom type {label!r}; types are {[t.label for t in self.atom_types]}"
        )

    @property
    def n_molecules(self) -> int:
        return sum(count for _, count in self.components)

    @property
    def n_atoms(self) -> int:
        return sum(molecule.n_atoms * count for molecule, count in self.components)

    @property
    def box_size(self) -> float:
        """Side of the cubic box (Å) at the given density."""
        return float((self.n_atoms / self.density) ** (1.0 / 3.0))

    def pair_parameters(self) -> PairTable:
        """The reference potential of every pair of atom types, in EPSR's pair order.

        Lorentz–Berthelot, as the engine mixes them (``ato.rs``): ε_ij = √(ε_i ε_j),
        σ_ij = (σ_i + σ_j)/2 and the Coulomb prefactor 1390 q_i q_j, with a
        ``terms`` column that shows which interactions are present (a type with
        ε = 0 has no Lennard-Jones site). After a run,
        ``result.potential.pair_params()`` gives the same table from the engine.
        """
        types = self.atom_types
        rows = []
        for i, first in enumerate(types):
            for second in types[i:]:
                epsilon = math.sqrt(first.epsilon * second.epsilon)
                charge_product = COULOMB_FACTOR * first.charge * second.charge
                rows.append(
                    PairParameters(
                        types=(first.label, second.label),
                        epsilon=epsilon,
                        sigma=0.5 * (first.sigma + second.sigma),
                        charge_product=charge_product,
                        terms=PairParameters.terms_of(epsilon, charge_product),
                    )
                )
        return PairTable(rows)

    def __str__(self) -> str:
        molecules = ", ".join(
            f"{count} × {molecule.name}" for molecule, count in self.components
        )
        if not self.components:
            return f"{self.name}: empty, {self.density} atoms/Å³, {self.temperature} K"
        return (
            f"{self.name}: {molecules}, {self.n_atoms} atoms, box {self.box_size:.2f} Å, "
            f"{self.density} atoms/Å³, {self.temperature} K, types "
            f"{', '.join(atom_type.label for atom_type in self.atom_types)}"
        )

    def __repr__(self) -> str:
        return f"<epsr.System {self}>"

    # ---------------------------------------------------------------- serialisation
    def to_dict(self) -> dict[str, Any]:
        """A JSON-serialisable dict (:meth:`from_dict` rebuilds the system)."""
        return {
            "name": self.name,
            "temperature": self.temperature,
            "density": self.density,
            "atom_types": [atom_type.to_dict() for atom_type in self.atom_types],
            "molecules": [
                {"molecule": molecule.to_dict(), "count": count}
                for molecule, count in self.components
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> System:
        types = {values["label"]: AtomType(**values) for values in data["atom_types"]}
        system = cls(data["name"], data["temperature"], data["density"])
        for entry in data["molecules"]:
            system.add(Molecule.from_dict(entry["molecule"], types), entry["count"])
        return system
