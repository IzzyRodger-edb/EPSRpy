"""Measured datasets: neutron samples and their isotopic composition.

A :class:`NeutronSample` is one measured data file and the isotopes of the sample
it was measured on. :meth:`epsr.Project.create` turns it into the sample's
``.NWTS.dat`` file, runs ``nwts`` to get its neutron weights and adds a ``data``
block for it to the ``.EPSR.inp`` file.

>>> samples = [
...     NeutronSample("D2O", "SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}),
...     NeutronSample("HDO", "SLS18500-HDO5050.mint01", isotopes={"HW": {1: 0.5, 2: 0.5}}),
...     NeutronSample("H2O", "SLS18498-H2O.mint01"),   # natural isotopes
... ]
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping, Sequence
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any

from . import tools
from .io import Mint01, read_mint01
from .settings import keywords

if TYPE_CHECKING:
    from .system import System

__all__ = [
    "NeutronSample",
    "check_sample_names",
    "NATURAL",
    "DATASET_KEYWORDS",
]

#: The composition of an atom type that is not given: the natural element.
NATURAL: Mapping[int, float] = {0: 1.0}
#: Keywords of a ``data`` block of the ``.EPSR.inp`` file set per sample (besides
#: ``datafile`` and ``wtsfile``), in file order.
DATASET_KEYWORDS = ("nrtype", "rshmin", "szeros", "tweak", "efilereq", "porodadd")
#: ``nrtype`` values EPSR knows (see the EPSR manual); 5 is a Gudrun ``.mint01`` file.
DATA_TYPES = (1, 2, 3, 4, 5)

_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,23}$")


class NeutronSample:
    """One neutron diffraction dataset and the isotopic composition of its sample.

    Parameters
    ----------
    name:
        The sample name (``"D2O"``): it names the files ``<name>.NWTS.dat`` and
        ``<name>.NWTStot.wts`` and selects the dataset by name in the results
        (``result.fit("D2O")``). Letters, digits, ``-`` and ``_``.
    data_file:
        The measured data. With ``nrtype=5`` (the default) a Gudrun ``.mint01``
        file, checked when the sample is created. The file is copied into the
        simulation folder, so its name must not contain blanks.
    isotopes:
        ``{atom type label: {mass number: fraction}}``, exactly EPSR's
        ``abundances`` pairs: ``{"HW": {2: 1.0}}`` is pure deuterium,
        ``{"HW": {1: 0.5, 2: 0.5}}`` a 50:50 H/D mixture. Mass number 0 is the
        natural element; types that are not given are natural. Fractions must sum
        to 1 and there are at most five isotopes per type. The mass numbers are
        checked against the engine's scattering-length table (Sears 1991) as soon
        as the sample meets a :class:`~epsr.system.System` (:meth:`check`), i.e.
        before any EPSR program runs.
    exchangeable:
        Atom type labels whose atoms exchange between molecules in this sample.
        ``None`` (the default) takes each :attr:`AtomType.exchangeable
        <epsr.system.AtomType.exchangeable>`; a list replaces them for this
        sample only.
    nrtype, rshmin, szeros, tweak, efilereq, porodadd:
        The ``.EPSR.inp`` keywords of this dataset (EPSR names). ``None`` (the
        default, except for ``nrtype``) takes the value given for all datasets in
        the project settings, else EPSR's default (``rshmin`` 0.7 Å, ``tweak`` 1,
        ``efilereq`` 1, ``szeros`` 0, ``porodadd`` 0).
    """

    def __init__(
        self,
        name: str,
        data_file: str | PathLike[str],
        isotopes: Mapping[str, Mapping[int, float]] | None = None,
        *,
        exchangeable: Iterable[str] | None = None,
        nrtype: int = 5,
        rshmin: float | None = None,
        tweak: float | None = None,
        efilereq: float | None = None,
        szeros: float | None = None,
        porodadd: float | None = None,
    ) -> None:
        if not isinstance(name, str) or not _NAME_PATTERN.match(name):
            raise ValueError(
                f"sample name {name!r}: use 1-24 letters, digits, '-' or '_' "
                "(it becomes a file name)"
            )
        self.name = name
        self.data_file = Path(data_file)
        if any(character.isspace() for character in self.data_file.name):
            raise ValueError(
                f"sample {name}: data file name {self.data_file.name!r} contains a "
                "blank; EPSR would read only the part before it. Rename the file."
            )
        if nrtype not in DATA_TYPES:
            raise ValueError(f"sample {name}: nrtype must be one of {DATA_TYPES}")
        self.nrtype = int(nrtype)
        self.data: Mint01 | None = None
        if self.nrtype == 5:
            self.data = read_mint01(self.data_file)
        elif not self.data_file.is_file():
            raise FileNotFoundError(
                f"sample {name}: data file {self.data_file} does not exist"
            )
        self.isotopes: dict[str, dict[int, float]] = {}
        for label, composition in (isotopes or {}).items():
            try:
                tools.abundances_text(composition)  # structure: ints, fractions, sum
            except (TypeError, ValueError) as error:
                raise type(error)(
                    f"sample {name}, atom type {label}: {error}"
                ) from None
            self.isotopes[label] = {int(m): float(f) for m, f in composition.items()}
        if isinstance(exchangeable, str):
            raise TypeError(
                f"sample {name}: exchangeable is a list of atom type labels"
            )
        self.exchangeable: tuple[str, ...] | None = (
            None if exchangeable is None else tuple(exchangeable)
        )
        self.rshmin = _optional(rshmin, f"{name}: rshmin", minimum=0.0)
        self.tweak = _optional(tweak, f"{name}: tweak")
        self.efilereq = _optional(efilereq, f"{name}: efilereq")
        self.szeros = _optional(szeros, f"{name}: szeros")
        self.porodadd = _optional(porodadd, f"{name}: porodadd")

    def __repr__(self) -> str:
        isotopes = ", ".join(
            f"{label} {tools.abundances_text(composition)}"
            for label, composition in self.isotopes.items()
        )
        return (
            f"<epsr.NeutronSample {self.name}: {self.data_file.name}, "
            f"{isotopes or 'natural isotopes'}>"
        )

    # ---------------------------------------------------------------- files
    @property
    def nwts_file(self) -> str:
        """``<name>.NWTS.dat``, the input of ``nwts``."""
        return f"{self.name}.NWTS.dat"

    @property
    def weights_file(self) -> str:
        """``<name>.NWTStot.wts``, the weights ``nwts`` writes (the dataset's ``wtsfile``)."""
        return f"{self.name}.NWTStot.wts"

    def composition(self, label: str) -> Mapping[int, float]:
        """``{mass number: fraction}`` of an atom type in this sample (natural if not given)."""
        return self.isotopes.get(label, NATURAL)

    def is_exchangeable(self, atom_type: Any) -> bool:
        """Whether ``atom_type`` (an :class:`~epsr.system.AtomType`) exchanges in this sample."""
        if self.exchangeable is None:
            return bool(atom_type.exchangeable)
        return atom_type.label in self.exchangeable

    def check(self, system: System) -> None:
        """Check the sample against a system: labels exist, isotopes are in the table.

        Raises ``ValueError`` (``KeyError`` for an unknown label) before any EPSR
        program runs; ``nwts`` would otherwise stop and wait for typed input.
        """
        labels = [atom_type.label for atom_type in system.atom_types]
        for label in [*self.isotopes, *(self.exchangeable or ())]:
            if label not in labels:
                raise KeyError(
                    f"sample {self.name}: atom type {label!r} is not in the system "
                    f"(types: {labels})"
                )
        for atom_type in system.atom_types:
            try:
                tools.abundances_text(
                    self.composition(atom_type.label), atom_type.element
                )
            except ValueError as error:
                raise ValueError(
                    f"sample {self.name}, atom type {atom_type.label}: {error}"
                ) from None

    def nwts_components(self, system: System) -> list[dict[str, Any]]:
        """The ``component`` blocks of the ``.NWTS.dat`` file, one per atom type in system order."""
        self.check(system)
        return [
            {
                "atom": atom_type.label,
                "atsymbol": atom_type.element,
                "iexchange": int(self.is_exchangeable(atom_type)),
                "abundances": dict(self.composition(atom_type.label)),
            }
            for atom_type in system.atom_types
        ]

    def render_nwts(self, system: System) -> str:
        """The ``<name>.NWTS.dat`` text ``nwts`` reads for this sample in ``system``."""
        return tools.render_nwts(self.name, system.name, self.nwts_components(system))

    def dataset_keywords(
        self,
        defaults: Mapping[str, Any] | None = None,
        overrides: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """The keywords of this sample's ``data`` block of the ``.EPSR.inp`` file.

        Each value comes from ``overrides`` (settings given for this dataset),
        else from the sample, else from ``defaults`` (settings given for all
        datasets), else from EPSR's keyword table.
        """
        table = keywords()
        values: dict[str, Any] = {
            "datafile": self.data_file.name,
            "wtsfile": self.weights_file,
        }
        for keyword in DATASET_KEYWORDS:
            value = getattr(self, keyword)
            if overrides and keyword in overrides:
                value = overrides[keyword]
            elif value is None:
                value = (defaults or {}).get(keyword)
            if value is None:
                default = table[keyword].default
                is_integer = table[keyword].value_type == "i"
                value = int(default) if is_integer else float(default)
            values[keyword] = value
        return values

    def notes(self, q_max: float | None = None) -> list[str]:
        """Remarks on the data (not errors): see :meth:`epsr.io.Mint01.warnings`.

        With ``q_max`` (``nq * qstep`` of the settings), also says when the data
        extend beyond the Q range EPSR uses (the rest is read and ignored).
        """
        if self.data is None:
            return []
        notes = self.data.warnings()
        highest = self.data.q_range[1]
        if q_max is not None and highest > q_max:
            notes.append(
                f"data extend to Q = {highest:.1f} 1/Å; EPSR uses them up to "
                f"nq * qstep = {q_max:g} 1/Å"
            )
        return notes

    # ---------------------------------------------------------------- serialisation
    def to_dict(self) -> dict[str, Any]:
        """A JSON-serialisable dict (:meth:`from_dict` rebuilds the sample)."""
        return {
            "name": self.name,
            "data_file": str(self.data_file),
            "isotopes": {
                label: [[mass, fraction] for mass, fraction in composition.items()]
                for label, composition in self.isotopes.items()
            },
            "exchangeable": (
                None if self.exchangeable is None else list(self.exchangeable)
            ),
            **{keyword: getattr(self, keyword) for keyword in DATASET_KEYWORDS},
        }

    @classmethod
    def from_dict(
        cls, data: Mapping[str, Any], *, data_file: str | PathLike[str] | None = None
    ) -> NeutronSample:
        """Rebuild from :meth:`to_dict`; ``data_file`` replaces the stored path."""
        return cls(
            data["name"],
            data_file if data_file is not None else data["data_file"],
            {
                label: {int(mass): fraction for mass, fraction in pairs}
                for label, pairs in data["isotopes"].items()
            },
            exchangeable=data["exchangeable"],
            **{keyword: data[keyword] for keyword in DATASET_KEYWORDS},
        )


def check_sample_names(samples: Sequence[NeutronSample]) -> None:
    """Sample names and data file names must be unique within a project."""
    names: set[str] = set()
    files: dict[str, str] = {}
    for sample in samples:
        if sample.name in names:
            raise ValueError(f"two samples are named {sample.name!r}")
        names.add(sample.name)
        other = files.setdefault(sample.data_file.name, sample.name)
        if other != sample.name:
            raise ValueError(
                f"samples {other} and {sample.name} use data files with the same name "
                f"{sample.data_file.name}; they would overwrite each other"
            )


def _finite(value: float, what: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{what} must be a finite number")
    return value


def _optional(
    value: float | None, what: str, minimum: float | None = None
) -> float | None:
    if value is None:
        return None
    value = _finite(value, what)
    if minimum is not None and value < minimum:
        raise ValueError(f"{what} must be >= {minimum:g}")
    return value
