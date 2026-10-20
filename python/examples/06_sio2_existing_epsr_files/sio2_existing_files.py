"""Refine amorphous silica from existing EPSR input files, run within the Python framework.

Run from the ``python/`` folder:

    uv run python examples/06_sio2_existing_epsr_files/sio2_existing_files.py           # ~12 minutes
    uv run python examples/06_sio2_existing_epsr_files/sio2_existing_files.py --quick   # ~1 minute

The starting point is a set of EPSR input files written for the Fortran EPSR, as an
existing project would have them: the one-atom templates ``Si.mol`` and ``O.mol``, the
weights file ``SiO2.NWTS.dat`` and the run-control file ``silica.EPSR.inp`` (in
``template/``, copied from ``validation/silica_bulk/``), plus the NIMROD measurement. EPSR's own setup programs
build the box (1500 Si + 3000 O at 0.0664 atoms/Å³, 298 K), called from Python, and
:class:`epsr.Project` then runs the schedule of ``validation/silica_bulk/steps.txt``:
equilibrate for 100 iterations on the reference potential, then two runs of 100
iterations with the empirical potential (EReq 50 kJ/mol).

This is the same calculation as example 09, which defines every input in Python
instead; the two folders are byte-identical. Everything goes to ``output/`` next to
this script.
"""

import argparse
import os
import shutil
from pathlib import Path

import epsr
from epsr import tools

EXAMPLE_FOLDER = Path(__file__).resolve().parent
REPOSITORY = EXAMPLE_FOLDER.parents[2]
EPSR_INPUTS = (
    EXAMPLE_FOLDER / "template"
)  # written for the Fortran EPSR (silica_bulk of validation/)
DATA = REPOSITORY / "experimental_data" / "silica_bulk" / "data"
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"

parser = argparse.ArgumentParser()
parser.add_argument("--quick", action="store_true", help="10 iterations per run")
options = parser.parse_args()
EQUILIBRATION = 11 if options.quick else 100  # includes the 1-iteration scratch start
REFINEMENT = 10 if options.quick else 100

# --- 1. Gather the existing EPSR files in a fresh folder --------------------------------
folder = OUTPUT / "epsr_run"
if folder.exists():
    shutil.rmtree(folder)
folder.mkdir(parents=True)
for name in ("Si.mol", "O.mol", "SiO2.NWTS.dat", "silica.EPSR.inp"):
    shutil.copy2(EPSR_INPUTS / name, folder / name)
shutil.copy2(DATA / "SiO2_NIMROD.mint01", folder)  # the data file named in the .inp

# --- 2. Build the box with EPSR's setup programs (validation/silica_bulk/steps.txt) -----
tools.makemole(folder, "Si")  # one-atom "molecules": Si.mol -> Si.ato
tools.makemole(folder, "O")
tools.mixato(folder, [("Si", 1500), ("O", 3000)], density=0.0664, box="silica")
tools.randomise(folder, "silica")  # random positions; no fmole: nothing to relax
tools.nwts(folder, "SiO2.NWTS.dat")  # natural isotopes -> SiO2.NWTStot.wts

# --- 3. Open the folder as a Project ----------------------------------------------------
project = epsr.Project(folder)  # finds silica.EPSR.inp
print(project.status())  # not run yet; 1 dataset (named after the data file)
settings = project.settings  # the keyword values in the files, by EPSR name
print({name: settings[name] for name in ("ntimes", "qmin")})
print(settings["data"]["SiO2_NIMROD"])  # rshmin 1.2 Å: below the 1.6 Å Si-O bond

# --- 4. Refine --------------------------------------------------------------------------
# A folder not made by Project.create does not know its system, so name the epsr.vis
# case (atom colours and sizes for the figures) when recording.
record = (
    epsr.vis.TrajectoryRecorder(case="silica_bulk") if epsr.vis.available() else False
)
equilibration = project.equilibrate(iterations=EQUILIBRATION, record=record)
#   the .inp asks for a scratch start (ireset 1, iinit 1, niter 1): 1 + 99 iterations
history = project.refine(ereq=50.0, runs=2, iterations=REFINEMENT, record=record)

# --- 5. Results -------------------------------------------------------------------------
result = project.results()
final = history.summaries[-1]
print()
print(f"R-factor: {equilibration.r_factor[-1]:.3e} on the reference potential")
print(f"          {final.r_factor:.3e} with the empirical potential")
print(
    f"Energy {final.energy:.1f} kJ/mol per atom, pressure {final.pressure:.1f} kbar, "
    f"EP amplitude {final.ep_amplitude:.1f} of 50 kJ/mol"
)
fit = result.fit("SiO2_NIMROD")
print(f"  {fit.name} fitted over Q = {fit.q[0]:.2f}-{fit.q[-1]:.2f} 1/Å")
for a, b, window in (
    ("Si", "O", (1.2, 2.0)),
    ("O", "O", (2.2, 3.0)),
    ("Si", "Si", (2.6, 3.6)),
):
    peak = epsr.analysis.first_peak(result.r, result.gr_inter(a, b), window)
    print(f"{a}-{b} first peak at {peak:.2f} Å")

if not record:
    print("(install the 'plot' extra for the figures)")
    raise SystemExit

# --- 6. Figures (epsr.vis) --------------------------------------------------------------
import matplotlib  # noqa: E402  (available: epsr.vis.available() is True)

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from epsr import vis  # noqa: E402

figures = OUTPUT / "figures"
trajectory = project.trajectory()  # every iteration recorded (both stages)
trajectory.save(OUTPUT / "trajectory")  # reload later: vis.Trajectory.load(...)
print(f"\n{trajectory.nsteps} EPSR iterations recorded; figures in {figures}")
vis.make_report(trajectory, figures, verbose=False)  # 01-09 + summary tables
axes = epsr.plot.fits(result)
axes.figure.savefig(figures / "fits.png", dpi=150)
plt.close("all")
