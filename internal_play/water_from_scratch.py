"""Refine liquid water against neutron data, with every input defined in Python.

A copy of ``python/examples/07_water_from_scratch`` for experimenting with. Run it from
this folder, with ``epsr-py`` installed (see README.md here):

    python water_from_scratch.py            # ~4 minutes
    python water_from_scratch.py --quick    # ~30 s
    python water_from_scratch.py --seed 42  # another starting box

1000 water molecules are built from scratch and fitted to three SANDALS measurements
(H2O, D2O and HDO 50:50), read from ``data/``. The schedule is that of
``validation/water/steps.txt``: equilibrate for 100 iterations on the reference
potential, then two runs of 100 iterations with the empirical potential (EReq 3 kJ/mol).
No EPSR input file is written by hand: the library writes them from the definitions
below and runs EPSR's own setup programs. Everything goes to ``output/`` next to this
script.
"""

import argparse
from pathlib import Path

import epsr
from epsr.samples import NeutronSample
from epsr.system import AtomType, Molecule, System

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"  # the three measured .mint01 files
REFERENCE_PARTIALS = HERE / "reference_partials"  # Soper (2000), for figure 10
OUTPUT = HERE / "output"

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
# Reference potential per atom type: Lennard-Jones epsilon (kJ/mol), sigma (Å) and charge (e).
# EPSR always derives each pair from its two types: eps_ij = sqrt(eps_i eps_j),
# sigma_ij = (sigma_i + sigma_j)/2, plus Coulomb 1390 q_i q_j / r. Per-pair LJ is not possible.
# HW has epsilon = 0, i.e. no LJ site (SPC convention), so OW-HW and HW-HW are Coulomb-only;
# the hydrogen bond comes from the O-H Coulomb attraction (and, once ereq > 0, the EP).
# system.pair_parameters() below prints this, with a "terms" column per pair.
OW = AtomType(
    "OW",
    element="O",
    epsilon=0.6503,
    sigma=3.165492,
    charge=-0.82,
    coarse_grain={"rmax": 1.55},  # as in the validated .ato
)
HW = AtomType(
    "HW",
    element="H",
    epsilon=0.0,
    sigma=0.0,
    charge=0.41,
    mass=2.0,
    exchangeable=True,  # labile H: exchanges between molecules in every sample
)

water = Molecule.from_geometry(  # builds the coordinates from the internal geometry
    "Water",
    centre=("O", OW),  # atoms are named within the molecule
    arms=[("H1", HW, 0.976), ("H2", HW, 0.976)],  # bond lengths (Å)
    angle=107.134,  # H1-O-H2 (degrees)
    vibtemp=65,
    angtemp=1,
    dihtemp=10,  # restraint stiffness ("temperatures")
)
# Intramolecular: bonds and angles are harmonic *restraints* to this geometry (angles as
# 1-3 distances). There are no per-bond force constants: stiffness is vibtemp (x angtemp
# for angles), scaled by sqrt(reduced mass)/r0 and the temperature.
# The same molecule written atom by atom (atoms named; bonds/angles by name or 0-based
# index); from_geometry above is exact, so it uses the precise H2 position:
#   Molecule("Water",
#            atoms=[("O", OW, (0, 0, 0)), ("H1", HW, (0.976, 0, 0)),
#                   ("H2", HW, (0.976 * cos(107.134°), 0.976 * sin(107.134°), 0))],
#            bonds=[("O", "H1"), ("O", "H2")], angles=[("H1", "O", "H2")])
# or from a file: Molecule.from_xyz("water.xyz", types={"O": OW, "H": HW})
#   (atoms then named element + running number: "O1", "H1", "H2")

# --- 2. The box ---------------------------------------------------------------------
system = System("waterbox", temperature=300.0, density=0.1)  # atoms/Å³, cubic box
system.add(water, 1000)
report(system)  # 1000 × Water, 3000 atoms, box 31.07 Å, types OW, HW
report(system.pair_parameters())  # one row per pair, with a "terms" column

# --- 3. The data --------------------------------------------------------------------
# isotopes: {atom type: {mass number: fraction}}, exactly EPSR's `abundances` pairs.
# Mass number 0 = the natural element; types left out are natural ({0: 1.0}).
# Fractions must sum to 1, at most 5 isotopes per type, and each mass number is checked
# against the engine's Sears (1991) scattering-length table before anything runs.
# Exchange comes from AtomType(exchangeable=True).
samples = [
    NeutronSample("D2O", DATA / "SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}),
    NeutronSample(
        "HDO", DATA / "SLS18500-HDO5050.mint01", isotopes={"HW": {1: 0.5, 2: 0.5}}
    ),
    NeutronSample("H2O", DATA / "SLS18498-H2O.mint01"),  # natural abundance
]

# --- 4. Build the EPSR project --------------------------------------------------------
project = epsr.Project.create(  # does not run EPSR yet
    OUTPUT / "epsr_run",
    system,
    samples,
    settings={"rshmin": 0.7},  # EPSR keyword names; checked against the keyword table
    relax_iterations=100,  # fmole
    seed=options.seed,  # None: EPSR's default (= seed 1), byte-identical
    overwrite=True,
)
report(project.status())  # not run yet; 1000 molecules; 3 datasets (D2O, HDO, H2O); ...

# --- 5. Refine ------------------------------------------------------------------------
record = epsr.vis.available()  # True when the 'plot' extra is installed
equilibration = project.equilibrate(iterations=EQUILIBRATION, record=record)
#   first call: scratch start (1 iteration, creates the .pcof) + 99 iterations, ereq = 0
history = project.refine(ereq=3.0, runs=2, iterations=REFINEMENT, record=record)
#   ereq switched on for the first run; the second run continues from the files

# --- 6. Results -----------------------------------------------------------------------
result = project.results()  # read-only Simulation of the finished folder
report()
report(f"R-factor: {equilibration.r_factor[-1]:.3e} on the reference potential")
report(f"          {history.r_factor[-1]:.3e} with the empirical potential")
for sample in project.samples:  # names from the project definition
    fit = result.fit(sample.name)  # select a dataset by name
    report(
        f"  {sample.name}: fitted over Q = {fit.q[0]:.2f}-{fit.q[-1]:.2f} 1/Å, "
        f"R = {result.summary.r_factor_for(sample.name):.2e}"
    )
for a, b, window, meaning in (
    ("OW", "OW", (2.0, 4.0), "nearest-neighbour waters"),
    ("OW", "HW", (1.5, 2.3), "hydrogen bond"),
):
    peak = epsr.analysis.first_peak(result.r, result.gr_inter(a, b), window)
    report(f"{a}-{b} first peak at {peak:.2f} Å ({meaning})")
report(result.potential.pair_params())  # mixed parameters as EPSR uses them


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
left.set_ylim(-400.0, 600.0)  # the Coulomb terms are hundreds of kJ/mol at 2-3 Å
left.set_title("reference potential (LJ + Coulomb)")
epsr.plot.potentials(result, which="empirical", ax=right)
right.set_title("empirical potential after refinement")
figure.savefig(figures / "12_reference_and_empirical_potential.png", dpi=150)
plt.close(figure)

reference = vis.load_reference_partials(
    {
        "OW-OW": REFERENCE_PARTIALS / "h2o298aoogor.txt",
        "OW-HW": REFERENCE_PARTIALS / "h2o298aohgor.txt",
        "HW-HW": REFERENCE_PARTIALS / "h2o298ahhgor.txt",
    }
)
vis.plots.partials_vs_reference(
    trajectory, reference, source="Soper (2000), 298 K"
).savefig(figures / "10_partials_vs_experiment.png", dpi=150)
plt.close("all")

# Hydrogen bonds per molecule in the final configuration (O···H < 2.45 Å, other molecule).
oxygen, hydrogen, _ = vis.metrics.neighbour_pairs(
    trajectory.positions[-1],
    trajectory.box(-1),
    trajectory.atom_labels,
    "OW",
    "HW",
    2.45,
)
molecule = trajectory.molecule_of_atom
hydrogen_bonds = (molecule[oxygen] != molecule[hydrogen]).sum()
report(
    f"  {2 * hydrogen_bonds / trajectory.meta['nmol']:.2f} hydrogen bonds per molecule"
)

if not options.no_3d:
    try:
        from epsr.vis.viewer import ConfigurationViewer
    except ImportError:
        report("  (install the 'viewer' extra for the 3D configuration image)")
    else:
        ConfigurationViewer(
            trajectory, gr="total", slab=10.0, reference=reference, off_screen=True
        ).screenshot(figures / "11_configuration.png", step=-1)
        report("  3D configuration + partial g(r): 11_configuration.png")

write_summary()
