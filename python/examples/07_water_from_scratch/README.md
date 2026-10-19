# Water from scratch

EPSR refinement of liquid water against neutron data (SANDALS: H2O, D2O and HDO 50:50).
The only input files are the three measured data files
(`experimental_data/water/data/*.mint01`). [`water_from_scratch.py`](water_from_scratch.py)
defines everything else in Python with the setup API, and the library writes the EPSR
files from it and runs EPSR's own setup programs:

| input | defined in the script as | written by |
|---|---|---|
| reference potential (ε, σ, charge, mass per atom type) | `AtomType("OW", ...)`, `AtomType("HW", ...)` | the `potential` lines of `Water.mol` |
| molecule geometry and restraints | `Molecule.from_geometry("Water", ...)` | `Water.mol` → `makemole` |
| box (1000 molecules, 0.1 atoms/Å³, 300 K) | `System("waterbox", ...).add(water, 1000)` | `mixato`, `randomise`, `fmole` |
| neutron weights | `NeutronSample(..., isotopes={"HW": {2: 1.0}})` per sample | `<sample>.NWTS.dat` → `nwts` |
| EPSR control file | `settings={"rshmin": 0.7}` (EPSR keyword names) | `waterbox.EPSR.inp`, every other keyword at EPSR's default |

Run from the `python/` folder (the `plot` extra gives the figures, `viewer` the 3D image):

```bash
uv run python examples/07_water_from_scratch/water_from_scratch.py            # ~4 minutes
uv run python examples/07_water_from_scratch/water_from_scratch.py --quick    # 10 iterations per run
uv run python examples/07_water_from_scratch/water_from_scratch.py --seed 42  # another starting box
```

The schedule is that of `validation/water/steps.txt`: `project.equilibrate(iterations=100)`
is the scratch start (1 iteration, which creates the `.pcof` file) plus 99 iterations on
the reference potential, and `project.refine(ereq=3.0, runs=2, iterations=100)` two runs
of 100 iterations with the empirical potential (EReq 3 kJ/mol).

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the EPSR folder (`.ato`, `.pcof`, `.EPSR.inp`, `.erg`, every `.EPSR.*`
  result), plus `Water.mol` and `epsr_project.json` (the Python definitions), which EPSR
  does not read.
- `trajectory/`: every iteration, recorded by `epsr.vis` (reload with `vis.Trajectory.load`).
- `figures/`: `01`–`09` and `summary.csv`/`summary.md` from `vis.make_report`; `10` the partial
  g(r) against the experimental partials of Soper (2000); `11` the 3D configuration; `12` the
  reference and empirical potentials; `fits.png` the final fit to each dataset.
- `summary.txt`: the console summary (system, pair parameters, R-factors, peaks).

The full run reproduces the validated water calculation (Fortran, Rust CLI and wrapper,
`validation/water`): every file of the committed `output/epsr_run/` is byte-identical to
the run's, except the date line of `waterbox.EPSR.out` and `waterbox.EPSR.qdr`. The slow
test `uv run pytest -m slow -k example_07` checks this. `--seed` gives a different
starting box and therefore a different (equally valid) run.

`mixato` prints `sh: 1: ./Water.ato: Permission denied`. This is harmless: without a
`system_commands.txt`, its file-listing command is empty (the Fortran does the same).
