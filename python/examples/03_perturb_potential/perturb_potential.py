"""Perturbing the EPSR potential from Python.

Run from the ``python/`` folder:

    uv run python examples/03_perturb_potential/perturb_potential.py    # ~10 s

Three ways to perturb the potential of a copy of the DMSO/water test case:

(a) the standard EPSR control keywords, set in the files with ``project.set``;
(b) direct, in-memory edits of the empirical-potential (EP) coefficients, then a run;
(c) a fixed EP with refinement switched off (Monte Carlo only).

Everything goes to ``output/`` next to this script.
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
WATER_HYDROGEN_BOND = ("O1", "H1")  # water oxygen ... water hydrogen

project = epsr.Project.from_template(TEMPLATE, OUTPUT / "epsr_run", overwrite=True)

# --- (a) Standard knobs: EP amplitude limit, feedback factor, refinement factor --------
# Written into the .pcof / .EPSR.inp files, as a user would edit them between runs.
project.set(ereq=15.0, feedback=0.8, potfac=(1.0, 1.0))
print("ereq now", project.settings["ereq"], "kJ/mol")

# --- (b) Direct manipulation, in memory, before running ----------------------------------
simulation = project.simulation()  # one EPSR run, read from the files
potential = simulation.potential
print(potential)
print("parameters:", potential.params)

original = potential.ep_coefficients()  # (n_pairs, ncoeffp - 1)
pair = potential.pair_index(*WATER_HYDROGEN_BOND)
print(
    f"O1-H1 EP minimum before: "
    f"{potential.empirical(*WATER_HYDROGEN_BOND).min():.3f} kJ/mol"
)
damped = original.copy()
damped[pair] *= 0.5  # damp the water hydrogen-bond EP
potential.set_ep_coefficients(damped)  # writes pcoefp and rebuilds the EP
print(
    f"O1-H1 EP minimum after damping: "
    f"{potential.empirical(*WATER_HYDROGEN_BOND).min():.3f} kJ/mol"
)
potential.perturb_ep(WATER_HYDROGEN_BOND, 0.002)  # add to every coefficient of a pair
potential.scale_ep(0.9)  # scale the whole EP (never beyond ereq)
sum_matches = np.allclose(
    potential.total(), potential.reference() + potential.empirical(), atol=1e-3
)
print("total = reference + empirical:", sum_matches)

summary = simulation.run()  # the modified EP is refined further and saved to the .pcof
print(
    f"after run: R-factor {summary.r_factor:.4e}, "
    f"EP amplitude {summary.ep_amplitude:.3f}"
)

# --- (c) Keep the EP fixed and only run Monte Carlo: potfac 0 disables refinement ------
before = project.results().potential.ep_coefficients()
project.run(potfac=(0.0, 0.0))  # persists: set potfac back to refine again
after = project.results().potential.ep_coefficients()
print("EP unchanged by a run with potfac 0:", np.allclose(after, before))

# Reference-potential parameters are read-only (they come from the .ato file).
for atom_type in potential.type_params:
    print(
        f"  {atom_type.label:3s} epsilon {atom_type.epsilon:.4f} kJ/mol  "
        f"sigma {atom_type.sigma:.3f} Å  charge {atom_type.charge:+.4f} e"
    )
