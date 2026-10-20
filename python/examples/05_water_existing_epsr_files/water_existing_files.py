"""Refine liquid water from existing EPSR input files, run within the Python framework.

Run from the ``python/`` folder:

    uv run python examples/05_water_existing_epsr_files/water_existing_files.py           # ~4 minutes
    uv run python examples/05_water_existing_epsr_files/water_existing_files.py --quick   # ~30 s

The starting point is a set of EPSR input files written for the Fortran EPSR, as an
existing project would have them: the molecule template ``Water.mol``, one
``.NWTS.dat`` weights file per sample and the run-control file ``waterbox.EPSR.inp``
(in ``template/``, copied from ``validation/water/``), plus the three SANDALS
measurements (H2O, D2O and HDO 50:50). EPSR's own setup programs build the box from them, called from Python, and
:class:`epsr.Project` then runs the schedule of ``validation/water/steps.txt``:
equilibrate for 100 iterations on the reference potential, then two runs of 100
iterations with the empirical potential (EReq 3 kJ/mol).

This is the same calculation as example 07, which defines every input in Python
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
)  # written for the Fortran EPSR (water of validation/)
DATA = REPOSITORY / "experimental_data" / "water" / "data"
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"
SAMPLES = ("D2O", "HDO", "H2O")  # the order of the data blocks in waterbox.EPSR.inp

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
for name in ["Water.mol", "waterbox.EPSR.inp", *(f"{s}.NWTS.dat" for s in SAMPLES)]:
    shutil.copy2(EPSR_INPUTS / name, folder / name)
for data_file in DATA.glob("*.mint01"):  # the data files named in waterbox.EPSR.inp
    shutil.copy2(data_file, folder / data_file.name)

# --- 2. Build the box with EPSR's setup programs (validation/water/steps.txt) -----------
tools.makemole(folder, "Water")  # Water.mol -> the one-molecule template Water.ato
tools.mixato(folder, [("Water", 1000)], density=0.1, box="waterbox")  # atoms/Å³
tools.randomise(folder, "waterbox")  # random positions and orientations
tools.fmole(folder, "waterbox", 100)  # relax the molecular geometry
for sample in SAMPLES:
    tools.nwts(folder, f"{sample}.NWTS.dat")  # -> <sample>.NWTStot.wts

# --- 3. Open the folder as a Project ----------------------------------------------------
project = epsr.Project(folder)  # finds waterbox.EPSR.inp
print(project.status())  # not run yet; 3 datasets (named after the data files)
settings = project.settings  # the keyword values in the files, by EPSR name
print({name: settings[name] for name in ("ntimes", "qmin", "feedback")})
print(settings["data"]["SLS18502-D2O"])  # the keywords of one data block

# --- 4. Refine --------------------------------------------------------------------------
# A folder not made by Project.create does not know its system, so name the epsr.vis
# case (atom colours and sizes for the figures) when recording.
record = epsr.vis.TrajectoryRecorder(case="water") if epsr.vis.available() else False
equilibration = project.equilibrate(iterations=EQUILIBRATION, record=record)
#   the .inp asks for a scratch start (ireset 1, iinit 1, niter 1): 1 + 99 iterations
history = project.refine(ereq=3.0, runs=2, iterations=REFINEMENT, record=record)

# --- 5. Results -------------------------------------------------------------------------
result = project.results()
print()
print(f"R-factor: {equilibration.r_factor[-1]:.3e} on the reference potential")
print(f"          {history.r_factor[-1]:.3e} with the empirical potential")
for name in project.dataset_names:
    fit = result.fit(name)
    print(
        f"  {name:18s} fitted over Q = {fit.q[0]:.2f}-{fit.q[-1]:.2f} 1/Å, "
        f"R = {result.summary.r_factor_for(name):.2e}"
    )
for a, b, window, meaning in (
    ("OW", "OW", (2.0, 4.0), "nearest-neighbour waters"),
    ("OW", "HW", (1.5, 2.3), "hydrogen bond"),
):
    peak = epsr.analysis.first_peak(result.r, result.gr_inter(a, b), window)
    print(f"{a}-{b} first peak at {peak:.2f} Å ({meaning})")

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
