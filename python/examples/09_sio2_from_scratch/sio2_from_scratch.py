"""Refine amorphous silica against neutron data, with every input defined in Python.

Run from the ``python/`` folder:

    uv run python examples/09_sio2_from_scratch/sio2_from_scratch.py            # ~12 minutes
    uv run python examples/09_sio2_from_scratch/sio2_from_scratch.py --quick    # ~1 minute
    uv run python examples/09_sio2_from_scratch/sio2_from_scratch.py --seed 42  # another starting box

or, with ``epsr-py`` installed by pip (see the top-level README), from anywhere in the
repository clone: ``python python/examples/09_sio2_from_scratch/sio2_from_scratch.py``.

1500 Si and 3000 O atoms (each a one-atom "molecule") are placed at random and fitted
to one NIMROD measurement of SiO2 glass. The schedule is that of
``validation/silica_bulk/steps.txt``: equilibrate for 100 iterations on the reference
potential, then two runs of 100 iterations with the empirical potential (EReq 50
kJ/mol). No EPSR input file is written by hand: the library writes them from the
definitions below and runs EPSR's own setup programs. Everything goes to ``output/``
next to this script.
"""

import argparse
import os
from pathlib import Path

import numpy as np

import epsr
from epsr.samples import NeutronSample
from epsr.system import AtomType, Molecule, System

EXAMPLE_FOLDER = Path(__file__).resolve().parent
REPOSITORY = EXAMPLE_FOLDER.parents[2]
DATA = REPOSITORY / "experimental_data" / "silica_bulk" / "data"
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"
EREQ = 50.0  # kJ/mol, as in the Dissolve silica example
SI_O_BOND_CUTOFF = 2.0  # Å, the first minimum of the Si-O g(r)

parser = argparse.ArgumentParser()
parser.add_argument("--quick", action="store_true", help="10 iterations per run")
parser.add_argument(
    "--seed", type=int, default=None, help="random seed (default: EPSR's own, = 1)"
)
parser.add_argument("--no-3d", action="store_true", help="no 3D configuration image")
options = parser.parse_args()
EQUILIBRATION = 11 if options.quick else 100  # includes the 1-iteration scratch start
REFINEMENT = 10 if options.quick else 100

summary_lines: list[str] = []


def report(text: object = "") -> None:
    """Print, and keep the line for ``output/summary.txt``."""
    print(text)
    summary_lines.extend(str(text).splitlines() or [""])


# --- 1. The atoms -------------------------------------------------------------------
# Reference potential of the Dissolve silica example: Lennard-Jones epsilon (kJ/mol),
# sigma (Å) and formal charges (e). Si has a small sigma, so Si-O is held apart by the
# O-O repulsion and the Si-O Coulomb attraction; system.pair_parameters() shows the pairs.
SI = AtomType(
    "Si",
    element="Si",
    epsilon=0.175,
    sigma=1.03,
    charge=2.0,
    coarse_grain={"rmax": 1.11},  # as in the validated .mol files
)
OX = AtomType(
    "O",
    element="O",
    epsilon=0.165,
    sigma=3.5,
    charge=-1.0,
    coarse_grain={"rmax": 1.55},
)
# A network glass has no molecules: every atom is its own one-atom "molecule", so EPSR
# moves atoms one by one and there are no intramolecular restraints (and no fmole).
silicon = Molecule.atom(SI)  # named after the type: Si.mol
oxygen = Molecule.atom(OX)

# --- 2. The box ---------------------------------------------------------------------
system = System("silica", temperature=298.0, density=0.0664)  # atoms/Å³ (2.2 g/cm³)
system.add(silicon, 1500).add(oxygen, 3000)  # the add order fixes the type order
report(system)  # 4500 atoms, box 40.78 Å, types Si, O
report(system.pair_parameters())

# --- 3. The data --------------------------------------------------------------------
samples = [NeutronSample("SiO2", DATA / "SiO2_NIMROD.mint01")]  # natural isotopes

# --- 4. Build the EPSR project --------------------------------------------------------
project = epsr.Project.create(  # does not run EPSR yet
    OUTPUT / "epsr_run",
    system,
    samples,
    settings={
        "qmin": (0.1, 0.0),  # Q below 0.1 1/Å left out of the fits
        "rshmin": 1.2,  # background subtraction below 1.2 Å (the Si-O bond is 1.6 Å)
    },
    seed=options.seed,  # None: EPSR's default (= seed 1), byte-identical
    overwrite=True,
)
report(project.status())  # not run yet; 4500 molecules; 1 dataset (SiO2); ...

# --- 5. Refine ------------------------------------------------------------------------
record = epsr.vis.available()  # True when the 'plot' extra is installed
equilibration = project.equilibrate(iterations=EQUILIBRATION, record=record)
#   first call: scratch start (1 iteration, creates the .pcof) + 99 iterations, ereq = 0
history = project.refine(ereq=EREQ, runs=2, iterations=REFINEMENT, record=record)
#   ereq switched on for the first run; the second run continues from the files

# --- 6. Results -----------------------------------------------------------------------
result = project.results()  # read-only Simulation of the finished folder
final = history.summaries[-1]
report()
report(f"R-factor: {equilibration.r_factor[-1]:.3e} on the reference potential")
report(f"          {final.r_factor:.3e} with the empirical potential")
report(
    f"Energy {final.energy:.1f} kJ/mol per atom, pressure {final.pressure:.1f} kbar, "
    f"EP amplitude {final.ep_amplitude:.1f} of {EREQ:.0f} kJ/mol"
)
fit = result.fit("SiO2")
report(f"  SiO2: fitted over Q = {fit.q[0]:.2f}-{fit.q[-1]:.2f} 1/Å")
for a, b, window, meaning in (
    ("Si", "O", (1.2, 2.0), "the Si-O bond"),
    ("O", "O", (2.2, 3.0), "edge of the SiO4 tetrahedron"),
    ("Si", "Si", (2.6, 3.6), "corner-sharing tetrahedra"),
):
    peak = epsr.analysis.first_peak(result.r, result.gr_inter(a, b), window)
    report(f"{a}-{b} first peak at {peak:.2f} Å ({meaning})")
oxygen_density = result.density * 2 / 3  # O atoms per Å³
silicon_oxygen = epsr.analysis.coordination_number(
    result.r, result.gr_inter("Si", "O"), oxygen_density, SI_O_BOND_CUTOFF
)
report(
    f"O neighbours of Si within {SI_O_BOND_CUTOFF} Å: {silicon_oxygen:.2f} (4 in SiO2)"
)
# Is the refinement finished? The EP amplitude should have grown to ereq.
if final.ep_amplitude < 0.95 * EREQ:
    report(
        f"  not converged: the EP amplitude ({final.ep_amplitude:.1f} kJ/mol) is still "
        f"below ereq ({EREQ:.0f}); run more iterations with project.refine(...)"
    )


def write_summary() -> None:
    (OUTPUT / "summary.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )


if not record:
    report("(install the 'plot' extra for the figures)")
    write_summary()
    raise SystemExit

# --- 7. Figures (existing epsr.plot / epsr.vis API) -------------------------------------
import matplotlib  # noqa: E402  (available: epsr.vis.available() is True)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from epsr import vis  # noqa: E402

figures = OUTPUT / "figures"
trajectory = project.trajectory()  # every iteration recorded (both stages)
trajectory.save(OUTPUT / "trajectory")  # reload later: vis.Trajectory.load(...)
report(f"\n{trajectory.nsteps} EPSR iterations recorded; figures in {figures}")
vis.make_report(trajectory, figures, verbose=False)  # 01-09 + summary tables
fit_quality = vis.metrics.fit_metrics(trajectory)
report(
    f"  SiO2: R² of the fit {fit_quality['r2'][0, 0]:7.3f} (first iteration)"
    f" -> {fit_quality['r2'][-1, 0]:.3f} (last)"
)

axes = epsr.plot.fits(result)
axes.figure.savefig(figures / "fits.png", dpi=150)
plt.close(axes.figure)

figure, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
epsr.plot.potentials(result, which="reference", ax=left)
left.set_xlim(0.0, 12.5)  # the reference potential is truncated at rmaxpt = 12 Å
left.set_ylim(-1500.0, 1500.0)  # formal charges: Coulomb terms of ~1000 kJ/mol at 2 Å
left.set_title("reference potential (LJ + Coulomb)")
epsr.plot.potentials(result, which="empirical", ax=right)
right.set_title("empirical potential after refinement")
figure.savefig(figures / "12_reference_and_empirical_potential.png", dpi=150)
plt.close(figure)

# Formation of the tetrahedral network, every 5th iteration: the fraction of Si with
# four O neighbours, of O bridging two Si, and the mean tetrahedral order of Si.
labels = trajectory.atom_labels
sampled = np.unique(np.r_[np.arange(0, trajectory.nsteps, 5), trajectory.nsteps - 1])
fourfold_si, bridging_o, order_si = [], [], []
for step in sampled:
    positions, box = trajectory.positions[step], trajectory.box(step)
    _, fraction = vis.metrics.coordination_distribution(
        positions, box, labels, "Si", "O", SI_O_BOND_CUTOFF
    )
    fourfold_si.append(fraction[4] if len(fraction) > 4 else 0.0)
    _, fraction = vis.metrics.coordination_distribution(
        positions, box, labels, "O", "Si", SI_O_BOND_CUTOFF
    )
    bridging_o.append(fraction[2] if len(fraction) > 2 else 0.0)
    order_si.append(
        vis.metrics.tetrahedral_order(positions, box, labels, "Si", "O").mean()
    )
figure, axes = plt.subplots(1, 3, figsize=(13, 3.6))
iterations = trajectory.steps[sampled]
for axis, values, title in (
    (
        axes[0],
        fourfold_si,
        f"Fraction of 4-coordinated Si (O within {SI_O_BOND_CUTOFF} Å)",
    ),
    (axes[1], bridging_o, f"Fraction of bridging O (2 Si within {SI_O_BOND_CUTOFF} Å)"),
    (axes[2], order_si, "Mean tetrahedral order q of Si"),
):
    axis.plot(iterations, values, marker="o", markersize=3, color=vis.style.series(0))
    for first_step, _, _ in trajectory.run_starts[1:]:  # the start of each EPSR run
        axis.axvline(trajectory.steps[first_step] - 0.5, color="#8a8985", ls=":")
    axis.set_title(title, fontsize=9)
    axis.set_xlabel("EPSR iteration")
figure.suptitle("Silica: formation of the tetrahedral network", x=0.01, ha="left")
figure.tight_layout()
figure.savefig(figures / "10_network.png", dpi=150)
plt.close(figure)
report(
    f"  4-coordinated Si {100 * fourfold_si[-1]:.1f} %, bridging O "
    f"{100 * bridging_o[-1]:.1f} %, tetrahedral order <q> {order_si[-1]:.3f}"
)

if not options.no_3d:
    try:
        from epsr.vis.viewer import ConfigurationViewer
    except ImportError:
        report("  (install the 'viewer' extra for the 3D configuration image)")
    else:
        ConfigurationViewer(
            trajectory, slab=12.0, color_by="coordination", off_screen=True
        ).screenshot(figures / "11_configuration.png", step=-1)
        report("  3D configuration (Si coloured by O neighbours): 11_configuration.png")

write_summary()
