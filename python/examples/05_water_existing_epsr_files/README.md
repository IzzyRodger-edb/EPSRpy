# Water from existing EPSR files

The same refinement as [example 07](../07_water_from_scratch/README.md) (liquid water
against SANDALS H2O, D2O and HDO 50:50 data), but starting from **existing EPSR input
files**, as a project set up for the Fortran EPSR would have them. Nothing is defined
in Python: [`water_existing_files.py`](water_existing_files.py) copies the files, runs
EPSR's own setup programs on them and then refines with `epsr.Project`. Use this route
to bring an existing EPSR project into the Python framework.

| existing file (in `template/`, copied from `validation/water/`) | what it holds | turned into | by |
|---|---|---|---|
| `Water.mol` | geometry, restraints and reference potential of one molecule | `Water.ato` (the template) | `tools.makemole(folder, "Water")` |
| (none: the answers of `mixato`) | 1000 molecules at 0.1 atoms/Å³ | `waterbox.ato` (the box) | `tools.mixato(folder, [("Water", 1000)], density=0.1, box="waterbox")`, `tools.randomise`, `tools.fmole(..., 100)` |
| `D2O.NWTS.dat`, `HDO.NWTS.dat`, `H2O.NWTS.dat` | isotopes of each sample | `<sample>.NWTStot.wts` | `tools.nwts(folder, "<sample>.NWTS.dat")` |
| `waterbox.EPSR.inp` | the run control, set for a scratch start; names the data files | read by every EPSR run | `epsr.Project(folder)` |
| `SLS18498-H2O.mint01`... (`experimental_data/water/data/`) | the measured data | | |

Then the same calls as example 07: `project.equilibrate(iterations=100)` (scratch start +
99 iterations on the reference potential) and
`project.refine(ereq=3.0, runs=2, iterations=100)`. `project.status()` and
`project.settings` show what the files hold. The datasets are named after their data
files (`SLS18502-D2O`...), because the files carry no sample names. A folder that
`Project.create` did not build does not know its system, so the `epsr.vis` recorder is
given its case (`TrajectoryRecorder(case="water")`) for the atom colours in the figures.

Run from the `python/` folder:

```bash
uv run python examples/05_water_existing_epsr_files/water_existing_files.py           # ~4 minutes
uv run python examples/05_water_existing_epsr_files/water_existing_files.py --quick   # 10 iterations per run
```

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the EPSR folder, with the copied input files and everything EPSR wrote.
- `trajectory/`: every iteration, recorded by `epsr.vis` (reload with `vis.Trajectory.load`).
- `figures/`: `01`–`09` and `summary.csv`/`summary.md` from `vis.make_report`; `fits.png`
  the final fit to each dataset. Example 07 adds the comparison with the experimental
  partials, the hydrogen-bond count and the 3D configuration.

The R-factor falls from 6.1×10⁻⁴ (reference potential) to 2.2×10⁻⁴ (empirical
potential); the O–O peak is at 2.79 Å. Every file of `output/epsr_run/` is
byte-identical to the committed output of example 07 (and so to the validated Fortran
run), except the date line of `waterbox.EPSR.out` and `waterbox.EPSR.qdr`: both routes
build the same folder.
