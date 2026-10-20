"""Refine liquid benzene against neutron data, with every input defined in Python.

Run from the ``python/`` folder:

    uv run python examples/10_benzene_from_scratch/benzene_from_scratch.py            # ~10 minutes
    uv run python examples/10_benzene_from_scratch/benzene_from_scratch.py --quick    # ~1 minute
    uv run python examples/10_benzene_from_scratch/benzene_from_scratch.py --seed 42  # another starting box

or, with ``epsr-py`` installed by pip (see the top-level README), from anywhere in the
repository clone: ``python python/examples/10_benzene_from_scratch/benzene_from_scratch.py``.

200 benzene molecules, read from an XYZ file, are fitted to three NIMROD measurements
(C6H6, C6D6 and a 50:50 mixture of the two). The force field, density and EReq are
those of the Dissolve benzene example. The schedule follows the other liquids (examples 07 and 09):
equilibrate for 100 iterations on the reference potential, then two runs of 100
iterations with the empirical potential (EReq 20 kJ/mol). No EPSR input file is written
by hand: the library writes them from the definitions below and runs EPSR's own setup
programs. Everything goes to ``output/`` next to this script.
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
DATA = REPOSITORY / "experimental_data" / "benzene" / "data"
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"
AVOGADRO = 6.02214076e23

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


# --- 1. The molecule ----------------------------------------------------------------
# OPLS-AA aromatic carbon and hydrogen, as in the Dissolve benzene example: Lennard-
# Jones epsilon (kJ/mol), sigma (Å) and charge (e). Unlike water, H has its own LJ site.
CA = AtomType("CA", element="C", epsilon=0.29288, sigma=3.55, charge=-0.115)
HA = AtomType("HA", element="H", epsilon=0.12552, sigma=2.42, charge=0.115)

benzene = Molecule.from_xyz(  # atoms named element + running number: C1..C6, H1..H6
    EXAMPLE_FOLDER / "benzene.xyz",  # ideal D6h ring: C-C 1.397 Å, C-H 1.084 Å
    types={"C": CA, "H": HA},  # first column of the XYZ file -> atom type
    name="Benzene",
)  # bonds="auto": atoms closer than their covalent radii + 0.4 Å
report(benzene)  # 12 atoms, 12 bonds, 18 angles (every pair of bonds at an atom)
# The 18 angle restraints (120°) keep the ring planar and the H atoms in its plane: a
# hexagon with all six angles at 120° cannot pucker. No rotations: nothing to rotate.

# --- 2. The box ---------------------------------------------------------------------
# EPSR's density unit is atoms/Å³; convert from the mass density (0.876 g/cm³).
molar_mass = 6 * 12.011 + 6 * 1.008  # g/mol
density = 0.876 / molar_mass * AVOGADRO * 1e-24 * benzene.n_atoms  # atoms/Å³
system = System("benzene", temperature=300.0, density=density)
system.add(benzene, 200)
report(system)  # 200 × Benzene, 2400 atoms, box 30.94 Å, types CA, HA
report(system.pair_parameters())

# --- 3. The data --------------------------------------------------------------------
# isotopes: {atom type: {mass number: fraction}}; only the hydrogen changes. Benzene H
# does not exchange (AtomType exchangeable=False, the default).
samples = [
    NeutronSample("C6H6", DATA / "C6H6.mint01"),  # natural abundance
    NeutronSample("C6D6", DATA / "C6D6.mint01", isotopes={"HA": {2: 1.0}}),
    NeutronSample("HD5050", DATA / "5050.mint01", isotopes={"HA": {1: 0.5, 2: 0.5}}),
]
# The 50:50 sample is a mixture of whole C6H6 and C6D6 molecules; the weights treat H
# and D as mixed at random over every H site, which changes only the intramolecular
# H-H weights. (Two species with H types HA and DA would be exact, but the EP then
# refines HA and DA pairs separately and the two identical species partly demix.)

# --- 4. Build the EPSR project --------------------------------------------------------
project = epsr.Project.create(  # does not run EPSR yet
    OUTPUT / "epsr_run",
    system,
    samples,
    settings={"rshmin": 0.7},  # background subtraction below 0.7 Å (C-H is 1.08 Å)
    relax_iterations=100,  # fmole
    seed=options.seed,  # None: EPSR's default (= seed 1)
    overwrite=True,
)
report(project.status())  # not run yet; 200 molecules; 3 datasets (C6H6, C6D6, HD5050)

# --- 5. Refine ------------------------------------------------------------------------
record = epsr.vis.available()  # True when the 'plot' extra is installed
equilibration = project.equilibrate(iterations=EQUILIBRATION, record=record)
#   first call: scratch start (1 iteration, creates the .pcof) + 99 iterations, ereq = 0
history = project.refine(ereq=20.0, runs=2, iterations=REFINEMENT, record=record)
#   ereq switched on for the first run; the second run continues from the files

# --- 6. Results -----------------------------------------------------------------------
result = project.results()  # read-only Simulation of the finished folder
report()
report(f"R-factor: {equilibration.r_factor[-1]:.3e} on the reference potential")
report(f"          {history.r_factor[-1]:.3e} with the empirical potential")
report(
    f"Energy {history.summaries[-1].energy:.1f} kJ/mol per molecule (see the check below)"
)
for sample in project.samples:
    fit = result.fit(sample.name)
    report(
        f"  {sample.name}: fitted over Q = {fit.q[0]:.2f}-{fit.q[-1]:.2f} 1/Å, "
        f"R = {result.summary.r_factor_for(sample.name):.2e}"
    )
ring = epsr.analysis.first_peak(result.r, result.gr_intra("CA", "CA"), (1.0, 1.8))
report(f"CA-CA intramolecular first peak at {ring:.2f} Å (the ring, 1.397 Å target)")
for a, b, window in (("CA", "CA", (3.0, 7.0)), ("CA", "HA", (2.5, 6.0))):
    peak = epsr.analysis.first_peak(result.r, result.gr_inter(a, b), window)
    report(f"{a}-{b} intermolecular first peak at {peak:.2f} Å")

# Ring centres and ring normals of the final configuration, from the carbon positions.
box = result.box_size
positions = np.asarray(result.positions, dtype=float)
is_carbon = np.asarray(result.atom_labels) == "CA"
centres, normals = [], []
for atoms in result.molecule_atoms:
    carbons = positions[np.asarray(atoms)[is_carbon[atoms.start : atoms.stop]]]
    carbons -= box * np.round((carbons - carbons[0]) / box)  # whole molecule
    centre = carbons.mean(axis=0)
    centres.append(centre)
    normals.append(np.linalg.svd(carbons - centre)[2][-1])  # least-variance axis
centres, normals = np.array(centres), np.array(normals)
first, second = np.triu_indices(len(centres), k=1)
separation = centres[second] - centres[first]
separation -= box * np.round(separation / box)  # minimum image
distance = np.linalg.norm(separation, axis=1)
# The angle between the two ring planes (0°: parallel, 90°: perpendicular / T-shaped).
plane_angle = np.degrees(np.arccos(np.abs((normals[first] * normals[second]).sum(1))))
edges = np.arange(0.0, box / 2, 0.2)
counts, _ = np.histogram(distance, edges)
r_centre = 0.5 * (edges[1:] + edges[:-1])
shell = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
g_centre = counts / (shell * len(centres) / 2 * len(centres) / box**3)
inside = (r_centre > 4.0) & (r_centre < 8.0)
report(
    f"ring centre-centre first peak at {r_centre[inside][g_centre[inside].argmax()]:.1f} Å"
)
# Interlocked rings: the random start puts molecules through each other, and the Monte
# Carlo cannot always pull a ring off another (an H through the ring of its neighbour).
# Such pairs keep a huge repulsive energy, so the energy above is not meaningful while
# any are left; they are left out of the orientation statistics below.
interlocked = distance < 3.0
if interlocked.any():
    report(
        f"  check: {interlocked.sum()} pairs of rings interlocked (centres < 3 Å apart), "
        "left from the random start"
    )
for low, high in ((3.0, 5.0), (5.0, 7.5)):
    close = (distance >= low) & (distance < high)
    if close.any():
        report(
            f"  rings {low:.0f}-{high} Å apart: {close.sum()} pairs, "
            f"{100 * (plane_angle[close] < 30).mean():.0f} % near-parallel (< 30°), "
            f"{100 * (plane_angle[close] > 60).mean():.0f} % near-perpendicular (> 60°)"
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
for index, name in enumerate(project.dataset_names):
    report(
        f"  {name}: R² of the fit {fit_quality['r2'][0, index]:7.3f} (first iteration)"
        f" -> {fit_quality['r2'][-1, index]:.3f} (last)"
    )

axes = epsr.plot.fits(result)
axes.figure.savefig(figures / "fits.png", dpi=150)
plt.close(axes.figure)

figure, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
epsr.plot.potentials(result, which="reference", ax=left)
left.set_xlim(0.0, 12.5)  # the reference potential is truncated at rmaxpt = 12 Å
left.set_ylim(-5.0, 10.0)  # small charges: the LJ wells (~0.3 kJ/mol) and Coulomb tails
left.set_title("reference potential (LJ + Coulomb)")
epsr.plot.potentials(result, which="empirical", ax=right)
right.set_title("empirical potential after refinement")
figure.savefig(figures / "12_reference_and_empirical_potential.png", dpi=150)
plt.close(figure)

# Ring-ring structure: centre-centre g(r), and the angle between the ring planes of
# every pair of molecules against their separation.
figure, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
left.plot(r_centre, g_centre, color=vis.style.series(0))
left.axhline(1.0, color="#8a8985", ls=":")
left.set_xlabel("ring centre-centre distance (Å)")
left.set_ylabel("g(r) of ring centres")
left.set_title("benzene: ring centres (final configuration)")
near = distance < 9.0
image = right.hist2d(
    distance[near],
    plane_angle[near],
    bins=(36, 18),
    range=((3.0, 9.0), (0.0, 90.0)),
    cmap="viridis",
)
figure.colorbar(image[3], ax=right, label="pairs")
right.set_xlabel("ring centre-centre distance (Å)")
right.set_ylabel("angle between ring planes (°)")
right.set_title("0°: parallel (stacked), 90°: perpendicular (T-shaped)")
figure.savefig(figures / "10_ring_orientation.png", dpi=150)
plt.close(figure)

if not options.no_3d:
    try:
        from epsr.vis.viewer import ConfigurationViewer
    except ImportError:
        report("  (install the 'viewer' extra for the 3D configuration image)")
    else:
        ConfigurationViewer(
            trajectory, gr="total", slab=10.0, off_screen=True
        ).screenshot(figures / "11_configuration.png", step=-1)
        report("  3D configuration + total g(r): 11_configuration.png")

write_summary()
