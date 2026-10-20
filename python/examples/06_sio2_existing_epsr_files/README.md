# Silica from existing EPSR files

The same refinement as [example 09](../09_sio2_from_scratch/README.md) (amorphous silica
against NIMROD data), but starting from **existing EPSR input files**, as a project set
up for the Fortran EPSR would have them. Nothing is defined in Python:
[`sio2_existing_files.py`](sio2_existing_files.py) copies the files, runs EPSR's own
setup programs on them and then refines with `epsr.Project`. Use this route to bring an
existing EPSR project into the Python framework.

| existing file (in `template/`, copied from `validation/silica_bulk/`) | what it holds | turned into | by |
|---|---|---|---|
| `Si.mol`, `O.mol` | one-atom "molecules" and their reference potential | `Si.ato`, `O.ato` (the templates) | `tools.makemole(folder, "Si")`, `tools.makemole(folder, "O")` |
| (none: the answers of `mixato`) | 1500 Si + 3000 O at 0.0664 atoms/Å³ | `silica.ato` (the box) | `tools.mixato(folder, [("Si", 1500), ("O", 3000)], density=0.0664, box="silica")`, `tools.randomise` |
| `SiO2.NWTS.dat` | isotopes of the sample (natural) | `SiO2.NWTStot.wts` | `tools.nwts(folder, "SiO2.NWTS.dat")` |
| `silica.EPSR.inp` | the run control (`qmin 0.1`, `rshmin 1.2`), set for a scratch start; names the data file | read by every EPSR run | `epsr.Project(folder)` |
| `SiO2_NIMROD.mint01` (`experimental_data/silica_bulk/data/`) | the measured data | | |

Then the same calls as example 09: `project.equilibrate(iterations=100)` (scratch start +
99 iterations on the reference potential) and
`project.refine(ereq=50.0, runs=2, iterations=100)`. `project.status()` and
`project.settings` show what the files hold. The dataset is named after its data file
(`SiO2_NIMROD`). A folder that `Project.create` did not build does not know its system,
so the `epsr.vis` recorder is given its case (`TrajectoryRecorder(case="silica_bulk")`)
for the atom colours in the figures.

Run from the `python/` folder:

```bash
uv run python examples/06_sio2_existing_epsr_files/sio2_existing_files.py           # ~12 minutes
uv run python examples/06_sio2_existing_epsr_files/sio2_existing_files.py --quick   # 10 iterations per run
```

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the EPSR folder, with the copied input files and everything EPSR wrote.
- `trajectory/`: every iteration, recorded by `epsr.vis` (reload with `vis.Trajectory.load`).
- `figures/`: `01`–`09` and `summary.csv`/`summary.md` from `vis.make_report`; `fits.png`
  the final fit. Example 09 adds the formation of the tetrahedral network, the
  potentials and the 3D configuration.

The R-factor falls from 1.0×10⁻³ (reference potential) to 2.3×10⁻⁴ (empirical
potential). The first peaks are Si–O 1.59 Å, O–O 2.55 Å and Si–Si 3.09 Å. Every file of
`output/epsr_run/` is byte-identical to example 09's, except the date line of
`silica.EPSR.out`: both routes build the same folder, which is in turn byte-identical
to the original Fortran EPSR's for the same steps (`validation/silica_bulk/steps.txt`).
