"""Repeated EPSR runs with convergence monitoring, as the EPSR shell loop does.

Run from the ``python/`` folder:

    uv run python examples/02_refinement_loop/refinement_loop.py        # 10 runs, ~30 s
    uv run python examples/02_refinement_loop/refinement_loop.py --runs 3

Continues the refinement of a copy of the DMSO/water test case (``template/``, from ``EPSR-test/``) for a
number of EPSR runs. Each run is reported as it finishes, the loop stops early once the
fit is good enough, and the folder is copied every five runs so that an earlier state
can be restored. Everything goes to ``output/`` next to this script.
"""

import argparse
import os
from pathlib import Path

import epsr

EXAMPLE_FOLDER = Path(__file__).resolve().parent
TEMPLATE = (
    EXAMPLE_FOLDER / "template"
)  # the run files of EPSR-test/: 200 DMSO + 400 water
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"
TARGET_R_FACTOR = 1.0e-4

parser = argparse.ArgumentParser()
parser.add_argument("--runs", type=int, default=10, help="number of EPSR runs")
options = parser.parse_args()

project = epsr.Project.from_template(TEMPLATE, OUTPUT / "epsr_run", overwrite=True)


def report_progress(run_number: int, simulation: epsr.Simulation) -> None:
    """Called after every run; returning True would stop the loop."""
    summary = simulation.summary
    print(
        f"run {run_number:3d}  R {summary.r_factor:.4e}  E {summary.energy:8.3f}  "
        f"P {summary.pressure:6.3f}  EP {summary.ep_amplitude:.2f}"
    )


def target_reached(history: epsr.History) -> bool:
    """Stop early once the fit is good enough."""
    return history.r_factor[-1] < TARGET_R_FACTOR


# ereq is left as it is in the files (10 kJ/mol); each run has the files' niter (1)
# iterations of ntimes = 5 Monte Carlo cycles.
history = project.refine(
    runs=options.runs,
    ntimes=5,  # EPSR keyword: it persists in the files
    callback=report_progress,
    until=target_reached,
    snapshot_every=5,  # copies to output/epsr_run.snapshots/
)

print(
    f"{len(history)} runs; R-factor {history.r_factor[0]:.4e} -> "
    f"{history.r_factor[-1]:.4e}"
)
print("columns:", ", ".join(history.columns))
print("history written to", history.to_csv(OUTPUT / "history.csv"))
print("whole .erg record:", len(project.history()), "iterations")  # every run ever made

try:
    import matplotlib

    matplotlib.use("Agg")  # write the figure to a file, no display needed
    axes = epsr.plot.history(history)
    plot_file = OUTPUT / "history.png"
    axes[0].figure.savefig(plot_file, dpi=120)
    print("plot written to", plot_file)
except ImportError:
    print("(install the 'plot' extra for a convergence plot)")
