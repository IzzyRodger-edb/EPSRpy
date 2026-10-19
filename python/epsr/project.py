"""A simulation folder: building it from Python definitions, running and reading it.

:class:`Project` is an EPSR simulation folder. :meth:`Project.create` builds one
from a :class:`~epsr.system.System`, its :class:`~epsr.samples.NeutronSample`
datasets and :class:`~epsr.settings.Settings`, by writing the files EPSR's own
setup programs read and running those programs. The refinement then follows a
fixed vocabulary:

>>> project = epsr.Project.create("work/water", system, samples, settings={"rshmin": 0.7})
>>> project.equilibrate(iterations=100)        # scratch start + 99, reference potential
>>> history = project.refine(ereq=3.0, runs=2, iterations=100)
>>> result = project.results()                 # a Simulation over the finished folder

An existing folder is opened with ``Project(folder)``; :meth:`Project.from_template`
copies one first.
"""

from __future__ import annotations

import json
import shutil
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ._elements import ELEMENT_MASSES
from .io import KeywordFile, Table, read_erg, read_keyword_file, read_table_series
from .simulation import PROJECT_FILE, Overrides, Simulation, Summary
from .simulation import dataset_names as _folder_dataset_names
from .simulation import find_inp

if TYPE_CHECKING:
    from .refine import History
    from .samples import NeutronSample
    from .settings import Settings
    from .system import System

__all__ = ["Project", "ProjectStatus", "ProjectStateError", "copy_folder"]

#: Most datasets EPSR fits at once (``MFILE`` in ``state.rs``).
MAX_DATASETS = 30
#: Version of the ``epsr_project.json`` layout.
PROJECT_FILE_FORMAT = 1


class ProjectStateError(RuntimeError):
    """An action that does not fit the project's state (e.g. ``refine`` before ``equilibrate``)."""


def copy_folder(
    source: str | PathLike[str],
    destination: str | PathLike[str],
    *,
    overwrite: bool = False,
) -> Path:
    """Copy the regular files of ``source`` (not subfolders) into a new folder ``destination``.

    Refuses to touch an existing non-empty ``destination`` unless ``overwrite``
    is true, in which case ``destination`` is deleted first.
    """
    source, destination = Path(source), Path(destination)
    if not source.is_dir():
        raise FileNotFoundError(f"template folder {source} does not exist")
    _prepare_empty_folder(destination, overwrite)
    for source_file in source.iterdir():
        if source_file.is_file():
            shutil.copy2(source_file, destination / source_file.name)
    return destination


def _prepare_empty_folder(folder: Path, overwrite: bool) -> None:
    if folder.exists() and any(folder.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"{folder} exists and is not empty; pass overwrite=True to replace it"
            )
        shutil.rmtree(folder)
    folder.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class ProjectStatus:
    """Where a project stands (:meth:`Project.status`)."""

    folder: Path
    has_run: bool
    iterations: int  # iterations run so far (lines of the .erg file)
    ereq: float  # requested empirical-potential amplitude
    datasets: list[str]
    n_molecules: int | None
    seed: int | None
    r_factor: float | None  # raw R-factor of the last iteration (.erg)
    pending: Mapping[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        datasets = f"{len(self.datasets)} datasets ({', '.join(self.datasets)})"
        molecules = f"{self.n_molecules} molecules; " if self.n_molecules else ""
        seed = "default" if self.seed is None else str(self.seed)
        if not self.has_run:
            pending = (
                f"; pending .pcof settings {dict(self.pending)}" if self.pending else ""
            )
            return f"not run yet; {molecules}{datasets}; seed {seed}{pending}"
        r_factor = "" if self.r_factor is None else f"; R-factor {self.r_factor:.3e}"
        return (
            f"{self.iterations} iterations; ereq {self.ereq:g}{r_factor}; "
            f"{molecules}{datasets}; seed {seed}"
        )


class Project:
    """An EPSR simulation folder.

    >>> project = Project.from_template("EPSR-test", "work/dmso")  # copy, then use the copy
    >>> simulation = project.simulation(overrides={"ntimes": 5})
    >>> simulation.run()

    Every EPSR run rewrites ``.ato``, ``.pcof``, ``.EPSR.inp`` and ``.EPSR.inpa``
    in the folder and appends to ``.erg``. Always work on a copy of your
    original folder; :meth:`snapshot` keeps intermediate states.

    A folder built by :meth:`create` also holds ``epsr_project.json``, the Python
    definitions it was built from (:attr:`system`, :attr:`samples`); EPSR does not
    read it.
    """

    def __init__(
        self, folder: str | PathLike[str], inp_name: str | None = None
    ) -> None:
        self.folder = Path(folder)
        if not self.folder.is_dir():
            raise FileNotFoundError(f"simulation folder {self.folder} does not exist")
        self.inp_name = inp_name if inp_name is not None else find_inp(self.folder)
        self._recorder: Any = None
        #: Summary of the last run made through this object (see :meth:`results`).
        self.last_summary: Summary | None = None

    def __repr__(self) -> str:
        return f"<epsr.Project {self.folder / self.inp_name}>"

    @classmethod
    def from_template(
        cls,
        template: str | PathLike[str],
        work: str | PathLike[str],
        *,
        inp_name: str | None = None,
        overwrite: bool = False,
    ) -> Project:
        """Copy the files of ``template`` into ``work`` and open the copy."""
        return cls(copy_folder(template, work, overwrite=overwrite), inp_name)

    # ---------------------------------------------------------------- building a folder
    @classmethod
    def create(
        cls,
        folder: str | PathLike[str],
        system: System,
        samples: Sequence[NeutronSample],
        settings: Mapping[str, Any] | None = None,
        *,
        relax_iterations: int | None = None,
        randomise: bool | None = None,
        seed: int | None = None,
        overwrite: bool = False,
        verbose: bool = False,
    ) -> Project:
        """Build a simulation folder from Python definitions. EPSR does not run yet.

        The steps are those of EPSR's own setup, with the same programs. For a
        system of molecules:

        1. ``makemole`` for each molecule (its ``.mol`` file is written from the
           :class:`~epsr.system.Molecule` and kept in the folder);
        2. ``mixato``: the box ``<system.name>.ato``; then ``changeato`` if
           :meth:`System.monte_carlo <epsr.system.System.monte_carlo>` set values;
        3. the random seed, if given (see ``seed``);
        4. ``randomise``: random positions and orientations;
        5. ``fmole``: ``relax_iterations`` shakes of the molecular geometry (skipped
           when every molecule is a single atom, or with ``relax_iterations=0``).

        For a :class:`~epsr.system.Framework` (``System(structure=...)``):

        1. ``readpdb``: the PDB file (copied as ``<system.name>.pdb``, with EPSR's
           ``vanderWaalsRadii.txt``) becomes the box ``<system.name>.ato``;
        2. ``changeato``: the temperature, the Monte Carlo parameters and the
           ε, σ, charge (and mass, where it is not the element's) of every type;
        3. the random seed, if given. The crystal positions are kept: no
           ``randomise`` and no ``fmole`` unless asked for.

        Then, for both:

        6. ``nwts`` for each sample: ``<sample>.NWTS.dat`` → ``<sample>.NWTStot.wts``;
        7. the ``<system.name>.EPSR.inp`` run-control file, set for a scratch start.

        There is no ``.pcof`` file yet: the first :meth:`equilibrate` creates it.

        Parameters
        ----------
        folder:
            The simulation folder; created, and replaced if ``overwrite``.
        system, samples:
            What to simulate and the data to fit (at most 30 datasets).
        settings:
            EPSR keyword values (:class:`~epsr.settings.Settings`). ``.inp`` and
            dataset keywords go into the ``.EPSR.inp``; ``.pcof`` keywords
            (``rminpt``, ``psigma2``...) are kept and applied by the scratch start.
            ``ereq`` must be 0 here: the scratch start always has no empirical
            potential; switch it on with :meth:`refine`.
        relax_iterations:
            ``fmole`` iterations (0: none). Default: 100, or 0 for a framework.
        randomise:
            Randomise the box (``False`` keeps ``mixato``'s lattice). Default: true,
            or false for a framework.
        seed:
            Seed of EPSR's random-number generator (1 to 2**31 - 2), applied
            before the first random draw. ``None`` keeps EPSR's default, which is
            the same as ``seed=1``, so the folder is identical to one built by
            hand. A run is always reproducible: the generator state is saved in
            the box file after every program.
        overwrite:
            Replace a non-empty ``folder``.

        Notes
        -----
        The setup programs work in the folder as their working directory (as EPSR
        does), so one process builds one folder at a time; build several in
        parallel with separate processes.
        """
        from . import tools
        from .samples import NeutronSample, check_sample_names
        from .settings import Settings, value_text
        from .system import System

        if not isinstance(system, System):
            raise TypeError("Project.create needs an epsr.System")
        framework = system.structure
        if not system.components and framework is None:
            raise ValueError(f"system {system.name} has no molecules; use System.add")
        if randomise is None:
            randomise = framework is None
        if relax_iterations is None:
            relax_iterations = 0 if framework is not None else 100
        samples = list(samples)
        if not samples:
            raise ValueError("Project.create needs at least one sample to fit")
        if len(samples) > MAX_DATASETS:
            raise ValueError(
                f"{len(samples)} samples; EPSR fits at most {MAX_DATASETS}"
            )
        for sample in samples:
            if not isinstance(sample, NeutronSample):
                raise TypeError(
                    f"samples must be NeutronSample objects, got {sample!r}"
                )
            sample.check(system)
        check_sample_names(samples)
        names = [sample.name for sample in samples]
        settings = Settings(settings or {})
        settings.check_datasets(names)
        if "ereq" in settings and float(value_text("ereq", settings["ereq"])) != 0:
            raise ValueError(
                "ereq must be 0 when the project is created: the scratch start has no "
                "empirical potential. Switch it on with project.refine(ereq=...)."
            )
        if seed is not None:
            if isinstance(seed, bool) or not isinstance(seed, int):
                raise TypeError(f"seed must be an int, got {type(seed).__name__}")
            if not 1 <= seed <= tools.MAX_SEED:
                raise ValueError(
                    f"seed must be between 1 and {tools.MAX_SEED}, got {seed}"
                )
        if isinstance(relax_iterations, bool) or not isinstance(relax_iterations, int):
            raise TypeError("relax_iterations must be an int")
        if relax_iterations < 0:
            raise ValueError("relax_iterations must be >= 0")

        folder = Path(folder)
        _prepare_empty_folder(folder, overwrite)
        for sample in samples:
            shutil.copy2(sample.data_file, folder / sample.data_file.name)
        box = system.name
        if framework is not None:
            shutil.copy2(framework.path, folder / f"{box}.pdb")
            tools.copy_vdw_radii(folder)
            tools.readpdb(
                folder,
                f"{box}.pdb",
                framework.vdw_factor,
                dihedrals=framework.dihedrals,
                verbose=verbose,
            )
            commands = tools.changeato_commands(
                temperature=system.temperature,
                parameters=system.monte_carlo_parameters,
                atom_types=[atom_type.to_dict() for atom_type in system.atom_types],
                default_masses=ELEMENT_MASSES,
            )
            tools.changeato(folder, box, commands, verbose=verbose)
        else:
            for molecule in system.molecules:
                mol_file = folder / f"{molecule.name}.mol"
                mol_file.write_text(
                    molecule.render_mol(system.temperature), encoding="latin-1"
                )
                tools.makemole(folder, molecule.name, verbose=verbose)
            tools.mixato(
                folder,
                [(molecule.name, count) for molecule, count in system.components],
                system.density,
                box,
                verbose=verbose,
            )
            if system.monte_carlo_parameters:
                commands = tools.changeato_commands(
                    parameters=system.monte_carlo_parameters
                )
                tools.changeato(folder, box, commands, verbose=verbose)
        if randomise:
            tools.randomise(folder, box, seed=seed, verbose=verbose)
        elif seed is not None:  # the first random draw is then in the scratch run
            tools.set_rng_seed(folder / f"{box}.ato", seed)
        has_molecules = framework is not None or any(
            m.n_atoms > 1 for m in system.molecules
        )
        if relax_iterations and has_molecules:
            tools.fmole(folder, box, relax_iterations, verbose=verbose)
        for sample in samples:
            tools.nwts(
                folder,
                sample.nwts_file,
                box=box,
                components=sample.nwts_components(system),
                verbose=verbose,
            )
        datasets = [
            sample.dataset_keywords(
                settings.part("data"), settings.dataset_overrides(index, sample.name)
            )
            for index, sample in enumerate(samples)
        ]
        inp_name = f"{box}.EPSR.inp"
        (folder / inp_name).write_text(
            settings.render_inp(f"{box}.EPSR", box, datasets), encoding="latin-1"
        )
        pending = {k: v for k, v in settings.part("pcof").items() if k != "ereq"}
        record = {
            "format": PROJECT_FILE_FORMAT,
            "system": system.to_dict(),
            "samples": [
                {
                    **sample.to_dict(),
                    "data_file": sample.data_file.name,
                    "source": str(sample.data_file),
                }
                for sample in samples
            ],
            "settings": settings.to_dict(),
            "pending_pcof": Settings(pending).to_dict(),
            "seed": seed,
            "randomise": randomise,
            "relax_iterations": relax_iterations,
        }
        project = cls(folder, inp_name)
        project._write_record(record)
        return project

    # ---------------------------------------------------------------- the project file
    @property
    def project_file(self) -> Path:
        """``epsr_project.json``: how :meth:`create` built the folder (not read by EPSR)."""
        return self.folder / PROJECT_FILE

    def _record(self) -> dict[str, Any] | None:
        if not self.project_file.is_file():
            return None
        record: dict[str, Any] = json.loads(
            self.project_file.read_text(encoding="utf-8")
        )
        return record

    def _write_record(self, record: Mapping[str, Any]) -> None:
        self.project_file.write_text(
            json.dumps(record, indent=2) + "\n", encoding="utf-8"
        )

    def _require_record(self, what: str) -> dict[str, Any]:
        record = self._record()
        if record is None:
            raise ProjectStateError(
                f"{what}: {self.folder} was not made by Project.create (no {PROJECT_FILE})"
            )
        return record

    @property
    def system(self) -> System:
        """The :class:`~epsr.system.System` the folder was built from."""
        from .system import System

        data = self._require_record("system")["system"]
        structure_path = self.folder / f"{data['name']}.pdb"
        return System.from_dict(
            data, structure_path=structure_path if data.get("structure") else None
        )

    @property
    def samples(self) -> list[NeutronSample]:
        """The samples the folder was built from, with their data files in the folder."""
        from .samples import NeutronSample

        return [
            NeutronSample.from_dict(entry, data_file=self.folder / entry["data_file"])
            for entry in self._require_record("samples")["samples"]
        ]

    @property
    def seed(self) -> int | None:
        """The random seed given to :meth:`create` (``None``: EPSR's default)."""
        record = self._record()
        return None if record is None else record.get("seed")

    @property
    def dataset_names(self) -> list[str]:
        """Dataset names in ``data n`` order (sample names, else data-file stems)."""
        return _folder_dataset_names(self.folder, self.inp_name)

    def _pending_pcof(self) -> Settings:
        from .settings import Settings

        record = self._record()
        return Settings.from_dict(record.get("pending_pcof", {}) if record else {})

    def _set_pending_pcof(self, pending: Settings) -> None:
        record = self._require_record("pending .pcof settings")
        record["pending_pcof"] = pending.to_dict()
        self._write_record(record)

    # ---------------------------------------------------------------- naming
    @property
    def run_name(self) -> str:
        """Prefix of the output files, e.g. ``200DMSO400Water.EPSR``."""
        return self.inp_name.removesuffix(".inp")

    def simulation(
        self,
        overrides: Overrides = None,
        *,
        settings: Mapping[str, Any] | None = None,
        verbose: bool = False,
    ) -> Simulation:
        """Read the folder into a new :class:`Simulation` (one EPSR run).

        ``settings`` (checked) or ``overrides`` (raw) change keywords for this run;
        EPSR saves them back to its files, so they persist.
        """
        return Simulation(
            self.folder,
            self.inp_name,
            overrides=overrides,
            settings=settings,
            verbose=verbose,
        )

    def snapshot(self, label: str, *, overwrite: bool = False) -> Path:
        """Copy the folder's files to ``<folder>.snapshots/<label>/`` and return that path."""
        snapshot_folder = self.folder.parent / f"{self.folder.name}.snapshots" / label
        return copy_folder(self.folder, snapshot_folder, overwrite=overwrite)

    def output(self, extension: str) -> Path:
        """Path of an output file: ``output("g01")`` is ``<folder>/<run_name>.g01``."""
        return self.folder / f"{self.run_name}.{extension.lstrip('.')}"

    def read_table(self, extension: str) -> Table:
        """Read a result table (``"g01"``, ``"u01"``...).

        Continuation files (``02``, ``03``...) are merged in.
        """
        return read_table_series(self.output(extension))

    def history(self) -> History:
        """The ``.erg`` energy history as a :class:`epsr.History`."""
        from .refine import History

        return History.from_erg(self.output("erg"))

    def input_file(self) -> KeywordFile:
        """The parsed ``.EPSR.inp`` file."""
        return read_keyword_file(self.folder / self.inp_name)

    def files(self) -> dict[str, Path | list[Path]]:
        """The main input files named by the ``.EPSR.inp`` file."""
        input_file = self.input_file()
        paths: dict[str, Path | list[Path]] = {"inp": self.folder / self.inp_name}
        for key, keyword in (("ato", "fnameato"), ("pcof", "fnamepcof")):
            file_name = input_file.words(keyword)
            if file_name:
                paths[key] = self.folder / file_name[0]
        paths["data"] = [
            self.folder / block["datafile"].split()[0]
            for block in input_file.blocks
            if block.get("datafile")
        ]
        paths["weights"] = [
            self.folder / block["wtsfile"].split()[0]
            for block in input_file.blocks
            if block.get("wtsfile")
        ]
        paths["erg"] = self.output("erg")
        paths["out"] = self.output("out")
        return paths

    def _pcof_path(self) -> Path:
        words = self.input_file().words("fnamepcof")
        if not words:
            raise ProjectStateError(f"{self.inp_name} names no .pcof file (fnamepcof)")
        return self.folder / words[0]

    @property
    def has_run(self) -> bool:
        """True once EPSR has run: the ``.pcof`` file exists (the engine's own test)."""
        return self._pcof_path().is_file()

    # ---------------------------------------------------------------- running
    def _settings_lines(self, settings: Mapping[str, Any]) -> list[str]:
        from .settings import Settings

        return Settings(settings).override_lines(self.dataset_names)

    def _recorder_for(self, record: Any) -> Any:
        """The recorder to use for ``record`` (True, False or a TrajectoryRecorder)."""
        if record is False or record is None:
            return None
        from . import vis

        if not vis.available():
            raise ImportError("record=True needs the 'plot' extra (matplotlib, scipy)")
        if record is not True:
            if not isinstance(record, vis.TrajectoryRecorder):
                raise TypeError(
                    "record is True, False or an epsr.vis.TrajectoryRecorder"
                )
            self._recorder = record
            return record
        if self._recorder is None:
            self._recorder = vis.TrajectoryRecorder(case=self._vis_case())
        return self._recorder

    def _vis_case(self) -> str:
        """The ``epsr.vis`` case whose atom types are this system's (figure styling)."""
        from .vis.systems import CASES

        labels = set(self._atom_type_labels())
        for name, case in CASES.items():
            if set(case.element) == labels:
                return name
        return ""

    def _atom_type_labels(self) -> list[str]:
        """The atom-type labels of the system (from ``epsr_project.json``)."""
        record = self._record()
        if record is not None:
            return [atom_type["label"] for atom_type in record["system"]["atom_types"]]
        return []

    def _run_once(self, lines: list[str], record: Any, verbose: bool) -> Summary:
        return self._run_simulation(lines, record, verbose).summary

    def _run_simulation(
        self, lines: list[str], record: Any, verbose: bool
    ) -> Simulation:
        recorder = self._recorder_for(record)
        simulation = self.simulation(lines, verbose=verbose)
        if recorder is None:
            simulation.run()
        else:
            from . import vis

            vis.run_recorded(simulation, recorder)
        self.last_summary = simulation.summary
        return simulation

    def _require_run(self, what: str) -> None:
        if not self.has_run:
            raise ProjectStateError(
                f"{what}: EPSR has not run on this project yet. Start with "
                "project.equilibrate(iterations=...), which makes the scratch start "
                "(it creates the .pcof file) and equilibrates on the reference potential."
            )

    @staticmethod
    def _check_count(value: Any, what: str, minimum: int = 1) -> int:
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{what} must be an int, got {value!r}")
        if value < minimum:
            raise ValueError(f"{what} must be >= {minimum}, got {value}")
        return value

    def run(
        self,
        iterations: int | None = None,
        *,
        record: Any = False,
        verbose: bool = False,
        **settings: Any,
    ) -> Summary:
        """One EPSR run (``iterations`` refinements; default: the ``niter`` of the files).

        ``settings`` are EPSR keywords for this run (and later ones: EPSR saves
        them). ``record=True`` records every iteration (:meth:`trajectory`). On a
        project that has never run, raises :class:`ProjectStateError`: start with
        :meth:`equilibrate`.
        """
        self._require_run("run")
        lines = self._settings_lines(settings)
        if iterations is not None:
            lines.append(f"niter {self._check_count(iterations, 'iterations')}")
        return self._run_once(lines, record, verbose)

    def equilibrate(
        self,
        iterations: int,
        *,
        scratch_iterations: int = 1,
        record: Any = False,
        callback: Callable[[int, Simulation], bool | None] | None = None,
        verbose: bool = False,
        **settings: Any,
    ) -> History:
        """Equilibrate on the reference potential only (``ereq = 0``).

        **On a new project** (no ``.pcof`` file yet) the first ``scratch_iterations``
        iterations are the *scratch start*: a run with ``ireset 1``/``iinit 1`` and
        no ``.pcof``, which creates the ``.pcof`` file and applies any pending
        ``.pcof`` settings of :meth:`create`. Then one run of
        ``iterations - scratch_iterations``. ``equilibrate(iterations=100)`` is the
        1 + 99 of the validated EPSR schedule. The returned
        :class:`~epsr.History` then has two rows, the scratch run first.

        **Later calls** make one run of ``iterations``. If the empirical potential
        is on (``ereq > 0``) it is switched off (``ereq 0``), which sets it to zero.

        ``callback(run_number, simulation)`` is called after each run (1: the
        scratch start); returning true after the scratch start skips the rest.
        """
        from .refine import History

        self._check_count(iterations, "iterations")
        history = History()
        if not self.has_run:
            self._check_count(scratch_iterations, "scratch_iterations")
            if iterations <= scratch_iterations:
                raise ValueError(
                    f"iterations ({iterations}) must be more than the "
                    f"scratch_iterations ({scratch_iterations}) on a new project"
                )
            pending = self._pending_pcof()
            lines = pending.override_lines() + self._settings_lines(settings)
            if scratch_iterations != 1:  # the .inp of a new project has niter 1
                lines.append(f"niter {scratch_iterations}")
            scratch = self._run_simulation(lines, record, verbose)
            history.append(scratch.summary)
            if pending:
                self._set_pending_pcof(type(pending)())
            if callback is not None and callback(1, scratch):
                return history
            lines = [f"niter {iterations - scratch_iterations}"]
        else:
            lines = self._settings_lines(settings)
            if self._current_ereq() != 0:
                lines.append("ereq 0.0")
            lines.append(f"niter {iterations}")
        simulation = self._run_simulation(lines, record, verbose)
        history.append(simulation.summary)
        if callback is not None:
            callback(len(history), simulation)
        return history

    def refine(
        self,
        ereq: float | None = None,
        runs: int = 1,
        iterations: int | None = None,
        *,
        record: Any = False,
        until: Callable[[History], bool] | None = None,
        callback: Callable[[int, Simulation], bool | None] | None = None,
        snapshot_every: int | None = None,
        verbose: bool = False,
        progress: bool | None = False,
        **settings: Any,
    ) -> History:
        """Refine the empirical potential: ``runs`` EPSR runs of ``iterations`` each.

        Parameters
        ----------
        ereq:
            Requested EP amplitude (kJ/mol), set before the first run; it persists.
            ``None`` keeps the current value.
        runs, iterations:
            Number of runs, and iterations per run (default: the files' ``niter``).
        record:
            ``True`` (or a ``TrajectoryRecorder``) records every iteration.
        until, callback, snapshot_every, progress:
            As for :func:`epsr.refine`: stop when ``until(history)`` is true or
            ``callback(run, simulation)`` returns true; copy the folder every
            ``snapshot_every`` runs; show a progress bar.
        settings:
            EPSR keywords for these runs (they persist).

        Raises :class:`ProjectStateError` on a project that has never run: the
        reference-potential stage (:meth:`equilibrate`) comes first.
        """
        from .refine import refine
        from .settings import value_text

        self._require_run("refine")
        self._check_count(runs, "runs")
        lines = self._settings_lines(settings)
        if iterations is not None:
            lines.append(f"niter {self._check_count(iterations, 'iterations')}")
        initial = None
        if ereq is not None:
            if float(ereq) < 0:
                raise ValueError(f"ereq must be >= 0, got {ereq}")
            initial = [f"ereq {value_text('ereq', ereq)}"]
        recorder = self._recorder_for(record)
        runner = None
        if recorder is not None:
            from . import vis

            def runner(simulation: Simulation) -> Summary:
                return vis.run_recorded(simulation, recorder)

        history = refine(
            self,
            runs,
            overrides=lines,
            initial_overrides=initial,
            callback=callback,
            until=until,
            snapshot_every=snapshot_every,
            verbose=verbose,
            progress=progress,
            runner=runner,
        )
        if history.summaries:
            self.last_summary = history.summaries[-1]
        return history

    def trajectory(self) -> Any:
        """Every iteration recorded with ``record=True`` on this project, as one ``Trajectory``.

        Equilibration and refinement appear together, in the order they ran.
        Needs the ``plot`` extra. Recording is in memory: save it with
        ``trajectory.save(folder)``.
        """
        if self._recorder is None:
            raise ProjectStateError(
                "nothing recorded: pass record=True to equilibrate/refine/run"
            )
        return self._recorder.trajectory()

    # ---------------------------------------------------------------- settings
    def _current_ereq(self) -> float:
        pcof = self._pcof_path()
        if not pcof.is_file():
            return 0.0
        words = read_keyword_file(pcof).words("ereq")
        return float(words[0]) if words else 0.0

    def set(self, **settings: Any) -> None:
        """Change EPSR keywords in the files without running (they apply from the next run).

        ``.EPSR.inp`` keywords and dataset keywords (globally, or per dataset with
        ``data={"D2O": {...}}``) are edited in the ``.EPSR.inp`` file; ``.pcof``
        keywords in the ``.pcof`` file, or, before the first run, kept for the
        scratch start. The edit is the one a user makes in a text editor between
        runs (``validation/set_keyword.py``).
        """
        from .settings import Settings, value_text

        values = Settings(settings)
        names = self.dataset_names
        values.check_datasets(names)
        inp = self.folder / self.inp_name
        for name, value in values.part("inp").items():
            _edit_keyword(inp, name, value_text(name, value))
        for index, dataset in enumerate(names):
            for name, value in values.dataset_values(index, dataset).items():
                _edit_keyword(inp, name, value_text(name, value), block=index + 1)
        pcof_values = values.part("pcof")
        if not pcof_values:
            return
        if self.has_run:
            for name, value in pcof_values.items():
                _edit_keyword(self._pcof_path(), name, value_text(name, value))
            return
        if (
            "ereq" in pcof_values
            and float(value_text("ereq", pcof_values["ereq"])) != 0
        ):
            raise ValueError(
                "ereq can only be set after the scratch start; use refine(ereq=...)"
            )
        pending = self._pending_pcof()
        pending.update({k: v for k, v in pcof_values.items() if k != "ereq"})
        self._set_pending_pcof(pending)

    @property
    def settings(self) -> Settings:
        """The current EPSR keyword values of the files (and pending ``.pcof`` ones).

        Values are the text EPSR wrote; datasets under ``"data"`` by name.
        """
        from .settings import MANAGED_KEYWORDS, Settings, keywords

        table = keywords()
        current = Settings()
        files = [read_keyword_file(self.folder / self.inp_name)]
        if self.has_run:
            files.append(read_keyword_file(self._pcof_path()))
        for keyword_file in files:
            for name, value in keyword_file.values.items():
                if name in table and name not in MANAGED_KEYWORDS and value:
                    current._values[name] = value
        if not self.has_run:
            current._values.update(self._pending_pcof()._values)
        blocks = files[0].blocks
        for name, block in zip(self.dataset_names, blocks, strict=False):
            current._datasets[name] = {
                key: value
                for key, value in block.items()
                if key in table and key not in MANAGED_KEYWORDS
            }
        return current

    def status(self) -> ProjectStatus:
        """Whether EPSR has run, iterations accumulated, ``ereq``, last R-factor, seed."""
        record = self._record()
        has_run = self.has_run
        iterations = 0
        r_factor = None
        erg = self.output("erg")
        if has_run and erg.is_file():
            energies = read_erg(erg)  # one line per iteration ever run
            iterations = len(energies)
            if iterations:
                r_factor = float(energies[-1, 2])
        n_molecules = None
        if record is not None:
            n_molecules = sum(entry["count"] for entry in record["system"]["molecules"])
        pending = self._pending_pcof() if not has_run else {}
        return ProjectStatus(
            folder=self.folder,
            has_run=has_run,
            iterations=iterations,
            ereq=self._current_ereq(),
            datasets=self.dataset_names,
            n_molecules=n_molecules,
            seed=None if record is None else record.get("seed"),
            r_factor=r_factor,
            pending=dict(pending),
        )

    # ---------------------------------------------------------------- results
    def results(self, *, verbose: bool = False) -> Simulation:
        """The finished folder as a :class:`Simulation`, for reading results.

        No file is changed: reading a folder makes EPSR rewrite its ``.qdr`` trend
        file, which is put back as it was. Do not :meth:`~Simulation.run` the
        returned object; use :meth:`run` or :meth:`refine` to continue.

        The arrays (fits, g(r), potentials) are those of the folder. Its
        ``summary`` is that of the last run made through this :class:`Project`
        object (:attr:`last_summary`); reading a folder runs no iteration, so a
        project opened afresh has no summary until it runs again (use
        :meth:`history` for the ``.erg`` record).
        """
        trend = self.output("qdr")
        saved = trend.read_bytes() if trend.is_file() else None
        try:
            simulation = Simulation(self.folder, self.inp_name, verbose=verbose)
            simulation._results_summary = self.last_summary
            return simulation
        finally:
            if saved is None:
                trend.unlink(missing_ok=True)
            else:
                trend.write_bytes(saved)


def _edit_keyword(
    path: Path, keyword: str, text: str, block: int | None = None
) -> None:
    """Set ``keyword`` to ``text`` in an EPSR keyword file, keeping its description.

    ``block`` selects a ``data n`` block of an ``.EPSR.inp`` file (1-based);
    otherwise the keywords before the first block are searched.
    """
    lines = path.read_text(encoding="latin-1").split("\n")
    in_block = block is None
    for index, line in enumerate(lines):
        words = line.split()
        if words == ["q"]:
            break
        if len(words) == 2 and words[0] == "data" and words[1].isdigit():
            if in_block and block is None:
                break
            in_block = int(words[1]) == block
            continue
        if in_block and words and words[0] == keyword:
            gap = line.find(" " * 15, len(keyword))
            description = line[gap:] if gap >= 0 else ""
            lines[index] = f"{keyword:<11} {text}{description}"
            path.write_text("\n".join(lines), encoding="latin-1")
            return
    where = f"data block {block}" if block is not None else "the keyword section"
    raise KeyError(f"{keyword!r} not found in {where} of {path.name}")
