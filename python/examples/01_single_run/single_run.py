"""One EPSR run on a copy of the DMSO/water test case, then read the results.

Run from the ``python/`` folder:

    uv run python examples/01_single_run/single_run.py    # a few seconds

``template/`` holds the run files of the DMSO/water test case ``EPSR-test/`` (200
DMSO + 400 water, three SANDALS datasets). It is a finished EPSR folder: its box,
potential, accumulated results and weights already exist, so EPSR can continue from it
straight away. The folder is copied first, because every EPSR run rewrites its input
files in place. Everything goes to ``output/`` next to this script.
"""

import os
from pathlib import Path

import numpy as np

import epsr

EXAMPLE_FOLDER = Path(__file__).resolve().parent
TEMPLATE = (
    EXAMPLE_FOLDER / "template"
)  # the run files of EPSR-test/: 200 DMSO + 400 water
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"

# --- 1. Open a copy of the folder ------------------------------------------------------
project = epsr.Project.from_template(TEMPLATE, OUTPUT / "epsr_run", overwrite=True)
print(project.status())  # already run; iterations accumulated; ereq; 3 datasets ...

# --- 2. One EPSR run ------------------------------------------------------------------
# EPSR keywords are given by their own names and checked against EPSR's keyword table;
# EPSR saves them back to its files, so they also apply to later runs.
summary = project.run(iterations=1, ntimes=5)  # releases the GIL while the engine runs
print(f"R-factor {summary.r_factor:.4e}")
print(
    f"energy   {summary.energy:.3f} kJ/mol per molecule, "
    f"pressure {summary.pressure:.3f} kbar"
)

# --- 3. Results -----------------------------------------------------------------------
result = project.results()  # read-only Simulation of the folder
print(f"Q grid: {result.q.size} points up to {result.q[-1]:.2f} 1/Å")
print(f"r grid: {result.r.size} points up to {result.r[-1]:.2f} Å")
for name in project.dataset_names:  # the data-file names of the test case
    fit = result.fit(name)  # select a dataset by name (or index)
    rms = float(np.sqrt(np.mean(fit.residual**2)))
    print(
        f"  {name:12s} Q {fit.q[0]:.2f}-{fit.q[-1]:.2f} 1/Å, rms(data - fit) = {rms:.4f}, "
        f"R = {summary.r_factor_for(name):.2e}"
    )

# Atom types: S, C, O, H belong to DMSO; O1, H1 to water.
print("pairs:", result.pair_labels[:6], "...")
for a, b, window, meaning in (
    ("O1", "H1", (1.5, 2.3), "water-water hydrogen bond"),
    ("O", "H1", (1.5, 2.3), "DMSO O ... water H"),
):
    peak = epsr.analysis.first_peak(result.r, result.gr_inter(a, b), window)
    print(f"{a}-{b} intermolecular g(r): first peak at {peak:.2f} Å ({meaning})")

# Partial g(r) and the potentials are (n_pairs, nr) float32 arrays.
empirical = result.potential.empirical()
print(f"EP range: {empirical.min():.3f} .. {empirical.max():.3f} kJ/mol")
print(
    f"configuration: {result.positions.shape} positions (Å), "
    f"box {result.box_size:.3f} Å"
)
