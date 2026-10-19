"""Analysis and visualisation of EPSR refinements, iteration by iteration (``epsr.vis``).

Needs the ``plot`` extra (matplotlib, scipy); the 3D viewer
(:mod:`epsr.vis.viewer`) also needs the ``viewer`` extra (pyvista, imageio).
``epsr.vis`` itself can always be imported: :func:`available` says whether the
extra is installed, and any other name raises ``ImportError`` without it.

EPSR writes its result files only after the last iteration of a run.  To follow
a refinement, record it: :func:`run_recorded` / :func:`refine_recorded` replace
``Simulation.run()`` / :func:`epsr.refine` and store a snapshot after every
iteration in a :class:`TrajectoryRecorder`; the simulation files are unchanged.

    from epsr import vis
    recorder = vis.TrajectoryRecorder(case="silica_bulk")
    vis.run_recorded(project.simulation(), recorder)
    vis.refine_recorded(project, 2, recorder, overrides={"niter": 100},
                        initial_overrides={"ereq": 50.0})
    trajectory = recorder.trajectory()        # or Trajectory.load(saved folder)
    vis.make_report(trajectory, "figures/")   # 9 figures + summary tables

Main entry points
-----------------
:class:`Trajectory`            per-iteration record (in memory, or loaded from disk)
:class:`TrajectoryRecorder`    collects snapshots; ``save()`` writes the folder format
:func:`run_recorded`           ``Simulation.run()``, recorded
:func:`refine_recorded`        :func:`epsr.refine`, recorded
:func:`make_report`            all figures + ``summary.csv``/``summary.md``
:mod:`epsr.vis.metrics`        R-factor, R², χ², coordination numbers, bond angles, ...
:mod:`epsr.vis.plots`          the individual matplotlib figures
:mod:`epsr.vis.viewer`         pyvista 3D configuration + partial g(r) (``viewer`` extra)

This package is the visualisation toolkit of
``visualisation/python_visualisation_toolkit`` (standalone package ``epsrvis``);
apart from this file and :mod:`epsr.vis.wrapper` the modules are identical.
"""

from __future__ import annotations

from typing import Any

try:
    from . import metrics, plots, style
except ImportError as error:  # matplotlib or scipy missing
    _missing = ImportError(
        "epsr.vis needs the 'plot' extra (matplotlib, scipy): "
        "uv sync --extra plot, or pip install 'epsr-py[plot]'"
    )
    _missing.__cause__ = error
    _MISSING_EXTRA: ImportError | None = _missing
else:
    _MISSING_EXTRA = None
    from .outputs import (
        OutputFolder,
        load_reference_partials,
        read_epsr_table,
        read_erg,
        read_isis_gr,
    )
    from .recorder import TrajectoryRecorder, record_run
    from .report import make_report
    from .trajectory import Trajectory
    from .wrapper import refine_recorded, run_recorded


def available() -> bool:
    """True when the ``plot`` extra (matplotlib, scipy) is installed and ``epsr.vis`` works.

    ``epsr.vis`` can always be imported; without the extra, every name other
    than :func:`available` raises ``ImportError`` when it is used.
    """
    return _MISSING_EXTRA is None


def __getattr__(name: str) -> Any:
    if _MISSING_EXTRA is not None and not name.startswith("__"):
        raise _MISSING_EXTRA
    raise AttributeError(f"module 'epsr.vis' has no attribute {name!r}")


__all__ = [
    "available",
    "OutputFolder",
    "Trajectory",
    "TrajectoryRecorder",
    "load_reference_partials",
    "make_report",
    "metrics",
    "plots",
    "read_epsr_table",
    "read_erg",
    "read_isis_gr",
    "record_run",
    "refine_recorded",
    "run_recorded",
    "style",
]
