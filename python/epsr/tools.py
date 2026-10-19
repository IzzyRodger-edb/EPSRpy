"""The auxiliary EPSR programs (setup and analysis tools).

All 35 programs of the EPSR suite are available through :func:`run`, with the
same arguments as on the command line (``epsr-rs <folder> <command> <args>``).

* **Keyword-file programs** (``nwts``, ``xwts``, ``partials``, ``coord``,
  ``sharm``, ...) take the name of their ``.dat`` file as argument and rewrite
  it with defaults filled in. :func:`write_keyword_file` creates such a file
  from a dict.
* **Interactive programs** (``makeato``, ``mixato``, ``fmole``, ``changeato``,
  ``readxyz``...) read their answers from standard input; pass them as
  ``answers``. The ``*.stdin`` files of the ``test_problems/*/EPSR-test-ex*``
  cases show the expected order of answers for each program.

Outputs are files in the folder; read the column files with :func:`epsr.io.read_table`.

**Typed setup helpers.** The programs that build a simulation from scratch have
wrappers that take Python values instead of prompt answers and check that the
expected file was written: :func:`makemole`, :func:`mixato`, :func:`randomise`
(with an optional random seed, see :func:`set_rng_seed`), :func:`fmole` and
:func:`nwts` with ``components``. :class:`epsr.Project` ``.create`` uses them; they
can also be called on their own.

**Working directory.** Like the Fortran, several programs (the box builders
``makeato``, ``mixato``, ``fmole``, ``addato``, ``randomise``, ``makelattice``,
``makelatticeato``) read ``system_commands.txt`` and write ``dirlist.txt`` in the
*current working directory*, because the EPSR shell runs every program from
inside the simulation folder. :func:`run` therefore changes the working directory
to ``folder`` for the duration of the call and passes ``./`` as the folder, exactly
as ``epsr-rs`` is used by the EPSR verification scripts. Because the working
directory is process-wide, tool runs are serialised by a lock (one tool at a time
per process), and other threads should not rely on relative paths while a tool
runs. :class:`epsr.Simulation` does not depend on the working directory and is
not affected.
"""

from __future__ import annotations

import math
import os
import shutil
import tempfile
import threading
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any

from . import _core
from .errors import EpsrError

__all__ = [
    "commands",
    "run",
    "keywords",
    "ToolKeyword",
    "ToolKeywords",
    "write_keyword_file",
    "render_keyword_file",
    "makemole",
    "mixato",
    "randomise",
    "set_rng_seed",
    "fmole",
    "nwts",
    "render_nwts",
    "isotopes",
    "abundances_text",
    "xwts",
    "partials",
    "coord",
    "triangles",
    "torangles",
    "writexyz",
]

#: Serialises tool runs: they change the process-wide working directory.
_WORKING_DIRECTORY_LOCK = threading.RLock()

#: Column layout of the keyword files EPSR writes (Fortran format ``(a11,1x,a,15x,a)``).
KEYWORD_COLUMN_WIDTH = 11
DESCRIPTION_GAP = " " * 15


def commands() -> list[str]:
    """Names of the 35 auxiliary programs."""
    return list(_core.tool_commands())


@contextmanager
def _working_directory(folder: Path) -> Iterator[None]:
    """Temporarily make ``folder`` the working directory (one thread at a time)."""
    with _WORKING_DIRECTORY_LOCK:
        previous_directory = os.getcwd()
        os.chdir(folder)
        try:
            yield
        finally:
            os.chdir(previous_directory)


def _standard_input_text(answers: str | Sequence[str]) -> str:
    """The answers as the text a program reads from standard input (newline-terminated)."""
    if isinstance(answers, str):
        if answers and not answers.endswith("\n"):
            return answers + "\n"
        return answers
    return "".join(f"{answer}\n" for answer in answers)


def run(
    folder: str | PathLike[str],
    command: str,
    args: Sequence[str] = (),
    answers: str | Sequence[str] = "",
    *,
    verbose: bool = False,
) -> None:
    """Run the auxiliary program ``command`` in ``folder``.

    Parameters
    ----------
    folder:
        Working folder; files are read from and written to it in place.
    command:
        Program name, see :func:`commands`.
    args:
        Words after the command name on the EPSR command line (normally the
        input file name), e.g. ``["DMSODWaterD.NWTS.dat"]``.
    answers:
        Standard input for interactive programs: a string, or a sequence of
        lines (joined with newlines).

    Runs with the working directory set to ``folder`` (see the module notes).
    Raises :class:`epsr.EpsrError` subclasses on failure. Releases the GIL.

    Like the Fortran programs, keyword-file tools do nothing when their input
    file is missing; the typed helpers (:func:`nwts`, :func:`partials`...)
    raise ``FileNotFoundError`` instead.
    """
    folder = Path(folder).resolve()
    if not folder.is_dir():
        raise FileNotFoundError(f"folder {folder} does not exist")
    arguments = [str(argument) for argument in args]
    with _working_directory(folder):
        _core.run_tool("./", command, arguments, _standard_input_text(answers), verbose)


@dataclass(frozen=True)
class ToolKeyword:
    """One entry of a program's keyword table."""

    name: str
    value_type: str  # "c" text, "i" integer, "r" real
    flag: int  # Fortran ``iflag``
    default: str
    description: str


@dataclass(frozen=True)
class ToolKeywords:
    """The keyword table of a keyword-file program (from ``symtab_tools.rs``)."""

    name: str
    extinp: str  # input-file extension, e.g. ".NWTS.dat"
    dataprompt: str  # name of the repeated block ("component", "distribution"...) or ""
    dataprompttext: str
    symbols: tuple[ToolKeyword, ...]
    #: Fortran markers, e.g. ``nsymboluptodata``: the number of keywords before the first block.
    marks: Mapping[str, int]

    @property
    def n_global(self) -> int:
        """Number of keywords before the first repeated block."""
        if self.dataprompt and "nsymboluptodata" in self.marks:
            return self.marks["nsymboluptodata"]
        return len(self.symbols)

    @property
    def global_keywords(self) -> list[str]:
        return [keyword.name for keyword in self.symbols[: self.n_global]]

    @property
    def block_keywords(self) -> list[str]:
        end_of_block = self.marks.get("nsymboldata", len(self.symbols))
        return [keyword.name for keyword in self.symbols[self.n_global : end_of_block]]

    def defaults(self) -> dict[str, str]:
        return {keyword.name: keyword.default for keyword in self.symbols}

    def descriptions(self) -> dict[str, str]:
        return {keyword.name: keyword.description for keyword in self.symbols}


def keywords(command: str) -> ToolKeywords:
    """The keyword table of ``command`` (``nwts``, ``partials``, ``coord``...)."""
    table = _core.tool_keywords(command)
    symbols = tuple(
        ToolKeyword(
            name=str(name),
            value_type=str(value_type),
            flag=int(flag),
            default=str(default),
            description=str(description),
        )
        for name, value_type, flag, default, description in table["symbols"]
    )
    return ToolKeywords(
        name=table["name"],
        extinp=table["extinp"],
        dataprompt=table["dataprompt"],
        dataprompttext=table["dataprompttext"],
        symbols=symbols,
        marks=dict(table["marks"]),
    )


def _format_value(value: Any) -> str:
    """Write a Python value the way EPSR keyword files expect it."""
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, (list, tuple)):
        return " ".join(_format_value(item) for item in value)
    return str(value)


def _check_keywords(given: Mapping[str, Any], allowed: list[str], what: str) -> None:
    unknown = set(given) - set(allowed)
    if unknown:
        raise ValueError(f"unknown {what} {sorted(unknown)}; valid: {allowed}")


def render_keyword_file(
    command: str,
    values: Mapping[str, Any],
    blocks: Sequence[Mapping[str, Any]] | None = None,
    *,
    title: str,
) -> str:
    """The text of a keyword input file for a keyword-file program.

    ``values`` sets global keywords; keywords not given take the table
    defaults. ``blocks`` gives one mapping per repeated block (``component``
    for ``nwts``, ...); keywords missing in a block take the defaults. Without
    ``blocks`` no block section is written, and programs that derive their
    blocks from another file (``nwts`` from the ``.ato`` file) fill them in.
    The text uses the layout EPSR writes (title, keyword, value, description,
    ``q``). See :func:`write_keyword_file`.
    """
    table = keywords(command)
    blocks = list(blocks or [])
    _check_keywords(values, table.global_keywords, f"{command} keywords")
    if blocks and not table.dataprompt:
        raise ValueError(f"{command} input files have no repeated blocks")
    for block in blocks:
        _check_keywords(block, table.block_keywords, f"{command} block keywords")
    defaults = table.defaults()
    descriptions = table.descriptions()

    def keyword_line(keyword: str, value: Any) -> str:
        return (
            f"{keyword:<{KEYWORD_COLUMN_WIDTH}} {_format_value(value)}"
            f"{DESCRIPTION_GAP}{descriptions[keyword]}"
        )

    lines = [f"{title}{DESCRIPTION_GAP}Title of this file"]
    for keyword in table.global_keywords:
        lines.append(keyword_line(keyword, values.get(keyword, defaults[keyword])))
    for block_number, block in enumerate(blocks, start=1):
        lines += ["", f"{table.dataprompt} {block_number:>3}", ""]
        for keyword in table.block_keywords:
            lines.append(keyword_line(keyword, block.get(keyword, defaults[keyword])))
    lines.append("q")
    return "\n".join(lines) + "\n"


def write_keyword_file(
    path: str | PathLike[str],
    command: str,
    values: Mapping[str, Any],
    blocks: Sequence[Mapping[str, Any]] | None = None,
    *,
    title: str | None = None,
) -> Path:
    """Write a keyword input file for a keyword-file program.

    See :func:`render_keyword_file` for ``values`` and ``blocks``. The title
    defaults to the file name without the program's extension. The program
    rewrites the file on exit.
    """
    path = Path(path)
    if title is None:
        extension = keywords(command).extinp
        title = (
            path.name.removesuffix(extension)
            if path.name.endswith(extension)
            else path.stem
        )
    path.write_text(
        render_keyword_file(command, values, blocks, title=title), encoding="latin-1"
    )
    return path


def _run_keyword_tool(
    command: str, folder: str | PathLike[str], dat_file: str, verbose: bool
) -> None:
    """Run a keyword-file program, refusing a missing input file.

    The Fortran program (and :func:`run`) silently does nothing in that case.
    """
    input_file = Path(folder) / dat_file
    if not input_file.is_file():
        raise FileNotFoundError(f"{command} input file {input_file} does not exist")
    run(folder, command, [dat_file], verbose=verbose)


#: Most isotopes ``nwts`` accepts for one atom type (``MISOTOPE`` in ``tools/nwts.rs``).
MAX_ISOTOPES = 5
#: Tolerance on the sum of the isotope fractions of one atom type.
FRACTION_SUM_TOLERANCE = 1e-6
#: Valid random seeds: 1 to 2**31 - 2 (the modulus of ``ran1`` minus one).
MAX_SEED = 2**31 - 2


def _base_name(name: str, extension: str) -> str:
    """``name`` without ``extension``; refuses blanks and paths (EPSR reads one word)."""
    base = name.removesuffix(extension)
    if not base or any(character.isspace() for character in base) or "/" in base:
        raise ValueError(f"{name!r} is not a valid EPSR file name (one word, no path)")
    return base


def _modification_time(path: Path) -> int | None:
    return path.stat().st_mtime_ns if path.exists() else None


def _check_written(path: Path, time_before: int | None, command: str) -> Path:
    """Raise if ``command`` did not (re)write ``path``.

    Several programs report a failure on their log and return normally, as the
    Fortran does; run with ``verbose=True`` to see that log.
    """
    if not path.exists() or _modification_time(path) == time_before:
        raise EpsrError(
            f"{command} did not write {path.name}; run with verbose=True to see its log"
        )
    return path


def _number_text(value: float) -> str:
    """A number as EPSR input text: ``1`` for whole numbers, else the shortest exact form."""
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"cannot write {value!r} to an EPSR file")
    return str(int(value)) if value.is_integer() else repr(value)


def makemole(
    folder: str | PathLike[str], molecule: str, *, verbose: bool = False
) -> Path:
    """Build the one-molecule template ``<molecule>.ato`` from ``<molecule>.mol``.

    Returns the path of the template. ``makemole`` also rewrites the ``.mol``
    file in its own layout. Raises ``FileNotFoundError`` if the ``.mol`` file is
    missing and :class:`epsr.EpsrError` if no template was written.
    """
    folder = Path(folder)
    name = _base_name(molecule, ".mol")
    if not (folder / f"{name}.mol").is_file():
        raise FileNotFoundError(f"makemole: {folder / name}.mol does not exist")
    template = folder / f"{name}.ato"
    time_before = _modification_time(template)
    run(folder, "makemole", [name], verbose=verbose)
    return _check_written(template, time_before, "makemole")


def mixato(
    folder: str | PathLike[str],
    components: Sequence[tuple[str, int]],
    density: float,
    box: str,
    *,
    verbose: bool = False,
) -> Path:
    """Build the box ``<box>.ato`` from molecule templates, as the ``mixato`` program.

    Parameters
    ----------
    components:
        ``(template, count)`` per molecule type, e.g. ``[("Water", 1000)]``. Each
        template ``<template>.ato`` must exist in ``folder`` (see :func:`makemole`).
        **The order fixes the atom-type order** of the box, and with it the order
        of every pair-resolved quantity.
    density:
        Atomic number density of the box in atoms/Å³. The box is cubic.
    box:
        Name of the box (without ``.ato``).

    Returns the path of the box ``.ato`` file. The answers given to the program
    are the number of types, the template and count of each, the density and
    the box name.
    """
    folder = Path(folder)
    box = _base_name(box, ".ato")
    if not components:
        raise ValueError("mixato: no molecule templates given")
    if not density > 0:
        raise ValueError(f"mixato: density must be > 0 atoms/Å³, got {density!r}")
    answers = [str(len(components))]
    for template, count in components:
        template = _base_name(template, ".ato")
        if template == box:
            raise ValueError(
                f"mixato: the box name {box!r} is also a template name; "
                "the box would overwrite the template"
            )
        if not (folder / f"{template}.ato").is_file():
            raise FileNotFoundError(
                f"mixato: template {folder / template}.ato does not exist"
            )
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ValueError(f"mixato: count of {template} must be an int >= 1")
        answers += [template, str(count)]
    answers += [repr(float(density)), box]
    box_file = folder / f"{box}.ato"
    time_before = _modification_time(box_file)
    run(folder, "mixato", answers=answers, verbose=verbose)
    return _check_written(box_file, time_before, "mixato")


def set_rng_seed(ato_file: str | PathLike[str], seed: int) -> None:
    """Seed the random-number generator stored in a box ``.ato`` file.

    The generator state (Numerical Recipes ``ran1``) is set to ``idum = -seed``,
    which ``ran1`` turns into a fresh shuffle table on its first draw. The file is
    read and written by the engine itself, so nothing else in it changes. Valid
    seeds are 1 to 2**31 - 2. A box made by ``mixato`` holds an all-zero state,
    which behaves exactly like seed 1, so not seeding at all is the same as
    ``seed=1``.
    """
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError(f"seed must be an int, got {type(seed).__name__}")
    _core.set_rng_seed(str(ato_file), seed)


def randomise(
    folder: str | PathLike[str],
    box: str,
    *,
    seed: int | None = None,
    verbose: bool = False,
) -> Path:
    """Give every molecule of ``<box>.ato`` a random position and orientation.

    ``seed`` seeds the box's random-number generator first (:func:`set_rng_seed`);
    ``None`` keeps the stored state, which for a new box equals seed 1.
    """
    folder = Path(folder)
    box_file = folder / f"{_base_name(box, '.ato')}.ato"
    if not box_file.is_file():
        raise FileNotFoundError(f"randomise: {box_file} does not exist")
    if seed is not None:
        set_rng_seed(box_file, seed)
    time_before = _modification_time(box_file)
    run(folder, "randomise", [box_file.stem], verbose=verbose)
    return _check_written(box_file, time_before, "randomise")


def fmole(
    folder: str | PathLike[str],
    box: str,
    iterations: int,
    *,
    verbose: bool = False,
) -> Path:
    """Relax the intramolecular geometry of every molecule of ``<box>.ato``.

    ``iterations`` Monte Carlo shakes on the bond and angle restraints only. The
    Fortran ``fmole`` first re-runs ``makemole`` on every template that has a
    ``.mol`` file in the folder; the Rust port does not, and refuses to run while
    such a file exists. The ``.mol`` files are therefore moved aside for the run
    and put back afterwards, which gives the Fortran result for templates that
    are up to date with their ``.mol`` file.
    """
    folder = Path(folder)
    box_file = folder / f"{_base_name(box, '.ato')}.ato"
    if not box_file.is_file():
        raise FileNotFoundError(f"fmole: {box_file} does not exist")
    if isinstance(iterations, bool) or not isinstance(iterations, int):
        raise TypeError("fmole: iterations must be an int")
    if iterations < 1:
        raise ValueError(f"fmole: iterations must be >= 1, got {iterations}")
    aside = Path(tempfile.mkdtemp(prefix=".fmole-", dir=folder))
    moved = [
        Path(shutil.move(path, aside / path.name)) for path in folder.glob("*.mol")
    ]
    try:
        time_before = _modification_time(box_file)
        run(folder, "fmole", [box_file.stem, str(iterations)], verbose=verbose)
        return _check_written(box_file, time_before, "fmole")
    finally:
        for path in moved:
            shutil.move(path, folder / path.name)
        aside.rmdir()


def isotopes(element: str) -> dict[int, float]:
    """Mass numbers known to ``nwts`` for ``element`` and their scattering lengths.

    From the Sears (1991) table compiled into the engine; lengths in 10⁻¹² cm,
    mass number 0 is the natural element. Empty for an unknown element.
    """
    return {
        mass_number: length
        for symbol, mass_number, length in _core.isotope_table()
        if symbol == element
    }


def abundances_text(
    composition: Mapping[int, float], element: str | None = None
) -> str:
    """The ``abundances`` value of a ``.NWTS.dat`` component: mass numbers and fractions.

    ``{2: 1.0}`` gives ``"2 1"``, ``{1: 0.5, 2: 0.5}`` gives ``"1 0.5 2 0.5"``, in
    the given order. Mass number 0 is the natural element. Checks that there are
    at most five isotopes and that the fractions sum to 1; with ``element``, also
    that every mass number is in the engine's isotope table (``nwts`` would wait
    for typed input otherwise).
    """
    if not composition:
        raise ValueError("isotopic composition is empty")
    if len(composition) > MAX_ISOTOPES:
        raise ValueError(
            f"at most {MAX_ISOTOPES} isotopes per atom type, got {len(composition)}"
        )
    known = isotopes(element) if element is not None else None
    if known is not None and not known:
        raise ValueError(f"element {element!r} is not in the isotope table")
    words = []
    for mass_number, fraction in composition.items():
        if isinstance(mass_number, bool) or not isinstance(mass_number, int):
            raise TypeError(f"mass numbers are ints (0 = natural), got {mass_number!r}")
        if known is not None and mass_number not in known:
            raise ValueError(
                f"{element} has no isotope with mass number {mass_number} in the "
                f"isotope table; known: {sorted(known)} (0 = natural)"
            )
        if not 0 <= float(fraction) <= 1:
            raise ValueError(f"isotope fraction {fraction!r} is not between 0 and 1")
        words += [str(mass_number), _number_text(fraction)]
    total = sum(float(fraction) for fraction in composition.values())
    if abs(total - 1.0) > FRACTION_SUM_TOLERANCE:
        raise ValueError(f"isotope fractions sum to {total:g}, not 1")
    return " ".join(words)


def _nwts_component(component: Mapping[str, Any]) -> dict[str, Any]:
    """One ``component`` block, with ``abundances`` given as a mapping or as text."""
    block = dict(component)
    abundances = block.get("abundances")
    if isinstance(abundances, Mapping):
        block["abundances"] = abundances_text(abundances, block.get("atsymbol"))
    if isinstance(block.get("iexchange"), bool):
        block["iexchange"] = int(block["iexchange"])
    return block


def render_nwts(sample: str, box: str, components: Sequence[Mapping[str, Any]]) -> str:
    """The ``<sample>.NWTS.dat`` text for a one-sample neutron measurement.

    ``components`` as for :func:`nwts`: one mapping per atom type of the box
    ``<box>.ato``, in box order.
    """
    if not components:
        raise ValueError("nwts: no components")
    return render_keyword_file(
        "nwts",
        {
            "fnameato": f"{_base_name(box, '.ato')}.ato",
            "nsamples": 1,
            "ncomponent": len(components),
        },
        [_nwts_component(component) for component in components],
        title=f"{_base_name(sample, '.NWTS.dat')}.NWTS",
    )


def nwts(
    folder: str | PathLike[str],
    dat_file: str,
    *,
    box: str | None = None,
    components: Sequence[Mapping[str, Any]] | None = None,
    verbose: bool = False,
) -> Path:
    """Neutron weights of one sample, from a ``.NWTS.dat`` file.

    Without ``components`` the file ``dat_file`` must exist in ``folder``
    (``FileNotFoundError`` otherwise). With ``components``, the file is written
    first from them: one mapping per atom type of the box, in box order, with
    the ``nwts`` keywords ``atom``, ``atsymbol``, ``iexchange`` (int or bool) and
    ``abundances``, the last either as text or as ``{mass_number: fraction}``
    (see :func:`abundances_text`). ``box`` names the box ``.ato`` file.

    Returns the path of the weights file of a one-sample input,
    ``<sample>.NWTStot.wts`` for ``<sample>.NWTS.dat``: the ``wtsfile`` of an
    EPSR dataset. When the file is written from ``components``, raises
    :class:`epsr.EpsrError` if ``nwts`` did not write it.
    """
    folder = Path(folder)
    weights_file = folder / f"{Path(dat_file).stem}tot.wts"
    if components is None:
        if box is not None:
            raise ValueError("nwts: box= is only used together with components=")
        _run_keyword_tool("nwts", folder, dat_file, verbose)
        return weights_file

    extension = keywords("nwts").extinp
    if not dat_file.endswith(extension):
        raise ValueError(
            f"nwts: write the input as <sample>{extension}, got {dat_file!r}"
        )
    if box is None:
        raise ValueError("nwts: pass box= with components=")
    sample = _base_name(dat_file, extension)
    (folder / dat_file).write_text(
        render_nwts(sample, box, components), encoding="latin-1"
    )
    time_before = _modification_time(weights_file)
    _run_keyword_tool("nwts", folder, dat_file, verbose)
    return _check_written(weights_file, time_before, "nwts")


def xwts(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
    """X-ray weights from a ``.XWTS.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("xwts", folder, dat_file, verbose)


def partials(
    folder: str | PathLike[str], dat_file: str, *, verbose: bool = False
) -> None:
    """Partial g(r) from a ``.PARTIALS.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("partials", folder, dat_file, verbose)


def coord(folder: str | PathLike[str], dat_file: str, *, verbose: bool = False) -> None:
    """Coordination numbers from a ``.COORD.dat`` file.

    Raises FileNotFoundError if the file is missing.
    """
    _run_keyword_tool("coord", folder, dat_file, verbose)


def triangles(
    folder: str | PathLike[str], dat_file: str, *, verbose: bool = False
) -> None:
    """Bond angles from a ``.TRIANGLES.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("triangles", folder, dat_file, verbose)


def torangles(
    folder: str | PathLike[str], dat_file: str, *, verbose: bool = False
) -> None:
    """Torsion angles from a ``.TORANGLES.dat`` file. Raises FileNotFoundError if it is missing."""
    _run_keyword_tool("torangles", folder, dat_file, verbose)


def writexyz(
    folder: str | PathLike[str], ato_file: str, answers: str | Sequence[str] = ""
) -> None:
    """Run ``writexyz`` on ``ato_file`` (answers as in the interactive program)."""
    run(folder, "writexyz", [ato_file], answers)
