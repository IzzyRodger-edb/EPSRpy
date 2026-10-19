"""EPSR settings: the keyword table of the engine and a validated mapping of values.

Settings use **EPSR's own keyword names** (``ntimes``, ``ereq``, ``qmin``, ``rshmin``...)
so that the EPSR manual, the engine's messages and the files all use the same words.
The table behind them is the engine's (``rust/epsr-core/src/symtab.rs``); it gives
each keyword's type, default, description and the file it is saved to:

* ``"inp"``: the header of the ``.EPSR.inp`` run-control file;
* ``"pcof"``: the header of the ``.pcof`` potential file (``ereq``, ``rminpt``...);
* ``"data"``: a ``data n`` block of the ``.EPSR.inp`` file, one per dataset
  (``rshmin``, ``tweak``...). Given globally, such a keyword applies to every dataset.

>>> settings = Settings({"rshmin": 0.7, "ntimes": 10})
>>> Settings.describe("feedback")
'feedback (.EPSR.inp, real, default 0.9): Confidence factor - should be < 1. [0.8]'

EPSR writes the settings back to its files at the end of every run, so a setting
**persists** into later runs.
"""

from __future__ import annotations

import difflib
import functools
from collections.abc import Iterator, Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from . import _core, _format

if TYPE_CHECKING:
    from .system import System

__all__ = ["Keyword", "Settings", "keywords", "suggest_rshmin", "MANAGED_KEYWORDS"]

#: Where each file kind of the table is saved, for messages.
FILE_NAMES = {"inp": ".EPSR.inp", "pcof": ".pcof", "data": ".EPSR.inp data block"}
VALUE_TYPES = {"i": "integer", "r": "real", "c": "text"}

#: Keywords the setup API sets itself, with the reason a user cannot set them.
MANAGED_KEYWORDS: dict[str, str] = {
    "fnameato": "the box file is named after the System",
    "fnamepcof": "the potential file is named after the System",
    "ndata": "the number of datasets is the number of samples",
    "datafile": "the data file is NeutronSample(data_file=...)",
    "wtsfile": "the weights file is made by nwts for each sample",
    "rho": "EPSR derives the density from the box",
    "ireset": "the scratch start sets it, and EPSR resets it after every run",
    "iinit": "the scratch start sets it, and EPSR resets it after every run",
    "nsumt": "EPSR counts the accumulated iterations itself",
    "niter": "give the number of iterations to run(), equilibrate() or refine()",
}


@dataclass(frozen=True)
class Keyword:
    """One entry of EPSR's keyword table."""

    name: str
    value_type: str  # "i" integer, "r" real, "c" text
    flag: int  # Fortran ``iflag``
    default: str
    description: str
    file: str  # "inp", "pcof" or "data"
    position: int  # 1-based position in the table

    def __str__(self) -> str:
        return (
            f"{self.name} ({FILE_NAMES[self.file]}, {VALUE_TYPES[self.value_type]}, "
            f"default {self.default}): {self.description}"
        )


@functools.lru_cache(maxsize=1)
def _table() -> dict[str, Keyword]:
    return {
        name: Keyword(name, value_type, flag, default, description, file, position)
        for position, (name, value_type, flag, default, description, file) in enumerate(
            _core.epsr_keywords()["symbols"], start=1
        )
    }


def keywords() -> dict[str, Keyword]:
    """EPSR's keyword table: name → :class:`Keyword`, in table order."""
    return dict(_table())


def _keyword(name: str) -> Keyword:
    table = _table()
    if name in table:
        return table[name]
    suggestions = difflib.get_close_matches(name, list(table), n=3)
    hint = (
        f"; did you mean {', '.join(repr(s) for s in suggestions)}?"
        if suggestions
        else ""
    )
    raise KeyError(f"{name!r} is not an EPSR keyword{hint}")


def _words(value: Any) -> list[Any]:
    if isinstance(value, str):
        return value.split()
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def value_text(name: str, value: Any) -> str:
    """The text of a keyword value as written to an EPSR file, checked against its type.

    Integers as integers; reals in their shortest exact form (``3`` is written
    ``3.0``); sequences (``(0.2, 5.0)``) and strings of several words joined by
    blanks.
    """
    keyword = _keyword(name)
    words = _words(value)
    if not words:
        raise ValueError(f"{name}: empty value")
    texts = []
    for word in words:
        if isinstance(word, bool):
            raise TypeError(f"{name}: {word!r} is not a valid value")
        if keyword.value_type == "c":
            text = str(word)
            if any(character.isspace() for character in text):
                raise ValueError(f"{name}: {word!r} contains a blank")
            texts.append(text)
            continue
        if isinstance(word, str):
            try:
                word = int(word) if keyword.value_type == "i" else float(word)
            except ValueError:
                raise ValueError(
                    f"{name} takes {VALUE_TYPES[keyword.value_type]} values, got {word!r}"
                ) from None
        if keyword.value_type == "i":
            if not isinstance(word, int) and not (
                isinstance(word, float) and word.is_integer()
            ):
                raise TypeError(f"{name} takes integer values, got {word!r}")
            texts.append(str(int(word)))
        else:
            if not isinstance(word, (int, float)):
                raise TypeError(f"{name} takes real values, got {word!r}")
            texts.append(repr(float(word)))
    return " ".join(texts)


class Settings(MutableMapping[str, Any]):
    """EPSR keyword values, checked against the engine's keyword table.

    A mapping from EPSR keyword names to values (a number, a string, or a tuple
    for keywords with several values such as ``qwidthqmax=(0.2, 5.0)``). Unknown
    names are refused with a suggestion. The key ``"data"`` holds per-dataset
    values: ``{"data": {"D2O": {"tweak": 1.05}}}``, keyed by sample name (or
    0-based dataset index). Keywords of a data block given globally apply to
    every dataset that does not set them itself.

    Keywords the setup API manages (:data:`MANAGED_KEYWORDS`: file names, the
    number of datasets, the run-state flags, ``niter``) are refused.

    Parameters
    ----------
    values:
        A mapping of keyword names to values.
    keyword_values:
        More values as keyword arguments (``Settings(ntimes=10)``).
    """

    def __init__(
        self, values: Mapping[str, Any] | None = None, /, **keyword_values: Any
    ) -> None:
        self._values: dict[str, Any] = {}
        self._datasets: dict[str | int, dict[str, Any]] = {}
        self.update(values or {})
        self.update(keyword_values)

    # ---------------------------------------------------------------- mapping
    def __setitem__(self, name: str, value: Any) -> None:
        if name == "data":
            self._set_datasets(value)
            return
        if name in MANAGED_KEYWORDS:
            raise KeyError(f"{name!r} cannot be set: {MANAGED_KEYWORDS[name]}")
        value_text(name, value)  # check
        self._values[name] = tuple(value) if isinstance(value, list) else value

    def _set_datasets(self, datasets: Mapping[str | int, Mapping[str, Any]]) -> None:
        if not isinstance(datasets, Mapping):
            raise TypeError(
                "'data' maps sample names (or 0-based indices) to keyword values"
            )
        for dataset, values in datasets.items():
            if isinstance(dataset, bool) or not isinstance(dataset, (str, int)):
                raise TypeError(
                    f"datasets are given by sample name or index, not {dataset!r}"
                )
            block = self._datasets.setdefault(dataset, {})
            for name, value in values.items():
                if name in MANAGED_KEYWORDS:
                    raise KeyError(f"{name!r} cannot be set: {MANAGED_KEYWORDS[name]}")
                if _keyword(name).file != "data":
                    raise KeyError(
                        f"{name!r} is not a per-dataset keyword; set it outside 'data'"
                    )
                value_text(name, value)
                block[name] = tuple(value) if isinstance(value, list) else value

    def __getitem__(self, name: str) -> Any:
        if name == "data":
            return {dataset: dict(values) for dataset, values in self._datasets.items()}
        return self._values[name]

    def __delitem__(self, name: str) -> None:
        if name == "data":
            self._datasets.clear()
        else:
            del self._values[name]

    def __iter__(self) -> Iterator[str]:
        yield from self._values
        if self._datasets:
            yield "data"

    def __len__(self) -> int:
        return len(self._values) + bool(self._datasets)

    def __repr__(self) -> str:
        items = dict(self._values)
        if self._datasets:
            items["data"] = self["data"]
        return f"Settings({items!r})"

    # ---------------------------------------------------------------- the table
    @staticmethod
    def describe(name: str) -> str:
        """The description of a keyword: file, type, default and EPSR's own text."""
        return str(_keyword(name))

    @staticmethod
    def file_of(name: str) -> str:
        """``"inp"``, ``"pcof"`` or ``"data"``: where EPSR saves the keyword."""
        return _keyword(name).file

    def part(self, file: str) -> dict[str, Any]:
        """The global values saved to one file kind (``"inp"``, ``"pcof"`` or ``"data"``)."""
        if file not in FILE_NAMES:
            raise ValueError(f"file must be one of {list(FILE_NAMES)}")
        return {
            name: value
            for name, value in self._values.items()
            if _keyword(name).file == file
        }

    def dataset_values(self, index: int, name: str | None = None) -> dict[str, Any]:
        """The per-dataset values for dataset ``index`` (0-based) of sample ``name``.

        Globally given data-block keywords first, then the values of ``data`` given
        for this sample's name or index.
        """
        values = self.part("data")
        for key in (index, name):
            if key is not None and key in self._datasets:
                values.update(self._datasets[key])
        return values

    def check_datasets(self, names: Sequence[str]) -> None:
        """Every dataset key of ``data`` must be one of ``names`` or a valid index."""
        for key in self._datasets:
            if isinstance(key, int):
                if not 0 <= key < len(names):
                    raise IndexError(
                        f"dataset index {key} out of range ({len(names)} datasets)"
                    )
            elif key not in names:
                raise KeyError(f"no dataset named {key!r}; datasets are {list(names)}")

    # ---------------------------------------------------------------- output
    def override_lines(self, names: Sequence[str] | None = None) -> list[str]:
        """The settings as ``.EPSR.inp`` lines for one run (``Simulation(overrides=...)``).

        Global ``.inp`` and ``.pcof`` keywords first, then one ``data n`` block per
        dataset (``names`` gives the sample names in dataset order, needed when
        data-block keywords are set).
        """
        lines = [
            f"{name} {value_text(name, value)}"
            for name, value in self._values.items()
            if _keyword(name).file != "data"
        ]
        has_data_values = bool(self.part("data") or self._datasets)
        if not has_data_values:
            return lines
        if names is None:
            raise ValueError("per-dataset settings need the dataset names")
        self.check_datasets(names)
        for index, name in enumerate(names):
            values = self.dataset_values(index, name)
            if values:
                lines.append(f"data {index + 1}")
                lines += [
                    f"{key} {value_text(key, value)}" for key, value in values.items()
                ]
        return lines

    def render_inp(
        self,
        run_name: str,
        box: str,
        datasets: Sequence[Mapping[str, Any]],
    ) -> str:
        """A complete ``.EPSR.inp`` for the scratch start (``ireset 1``, ``iinit 1``, ``nsumt -1``).

        Every keyword of the table in table order: the managed ones set by the
        setup API, the others from these settings or the table defaults. One
        ``data n`` block per entry of ``datasets`` (each with ``datafile``,
        ``wtsfile`` and the sample's data-block values). There is no ``.pcof`` yet;
        the scratch run creates it.
        """
        table = _table()
        fixed = {
            "fnameato": f"{box}.ato",
            "fnamepcof": f"{box}.pcof",
            "ireset": 1,
            "iinit": 1,
            "nsumt": -1,
            "niter": 1,
            "ndata": len(datasets),
        }
        lines = [f"{run_name}{_format.DESCRIPTION_GAP}{_format.TITLE_MARKER}"]
        for keyword in table.values():
            if keyword.file != "inp":
                continue
            if keyword.name in fixed:
                text = value_text(keyword.name, fixed[keyword.name])
            elif keyword.name in self._values:
                text = value_text(keyword.name, self._values[keyword.name])
            else:
                text = keyword.default
            lines.append(_format.keyword_line(keyword.name, text, keyword.description))
        for number, dataset in enumerate(datasets, start=1):
            lines += ["", f"data {number:3d}", ""]
            for keyword in table.values():
                if keyword.file != "data":
                    continue
                if keyword.name in dataset:
                    text = value_text(keyword.name, dataset[keyword.name])
                else:
                    text = keyword.default
                lines.append(
                    _format.keyword_line(keyword.name, text, keyword.description)
                )
        lines.append("q")
        return "".join(f"{line}\n" for line in lines)

    # ---------------------------------------------------------------- serialisation
    def to_dict(self) -> dict[str, Any]:
        """JSON-serialisable values (tuples become lists; dataset indices stay ints)."""
        values: dict[str, Any] = {
            name: list(value) if isinstance(value, tuple) else value
            for name, value in self._values.items()
        }
        if self._datasets:
            values["data"] = [
                [key, dict(block)] for key, block in self._datasets.items()
            ]
        return values

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> Settings:
        values = {name: value for name, value in data.items() if name != "data"}
        settings = cls(values)
        if "data" in data:
            settings["data"] = {key: block for key, block in data["data"]}
        return settings


#: ``suggest_rshmin`` values (Å): EPSR's default for molecular systems, and a value
#: below the shortest expected distance for systems of single atoms (silica Si-O 1.6 Å).
MOLECULAR_RSHMIN = 0.7
ATOMIC_RSHMIN = 1.2


def suggest_rshmin(system: System) -> float:
    """A starting value for ``rshmin``, the r below which g(r) counts as background.

    0.7 Å (EPSR's default) for molecular systems, less if a bond is shorter (0.2 Å
    below the shortest bond); 1.2 Å for systems of single atoms only. Never applied
    automatically: pass it yourself, ``settings={"rshmin": suggest_rshmin(system)}``.
    """
    bonds = [
        molecule.bond_length(*bond)
        for molecule in system.molecules
        for bond in molecule.bonds
    ]
    if not bonds:
        return ATOMIC_RSHMIN
    shortest = min(bonds)
    if shortest < MOLECULAR_RSHMIN:
        return round(shortest - 0.2, 3)
    return MOLECULAR_RSHMIN
