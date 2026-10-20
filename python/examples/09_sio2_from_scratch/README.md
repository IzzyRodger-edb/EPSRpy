# Silica from scratch

EPSR refinement of amorphous silica (SiO₂ glass) against the NIMROD data
(`experimental_data/silica_bulk/data/SiO2_NIMROD.mint01`). The only input file is the
measured data. [`sio2_from_scratch.py`](sio2_from_scratch.py) defines everything else in
Python with the setup API, and the library writes the EPSR files from it and runs
EPSR's own setup programs. Unlike water, a network glass has no molecules: each atom is
its own one-atom "molecule", moved on its own by the Monte Carlo.

| input | defined in the script as | written by |
|---|---|---|
| reference potential (ε, σ, formal charge per atom type) | `AtomType("Si", ...)`, `AtomType("O", ...)` | the `potential` lines of `Si.mol`, `O.mol` |
| one-atom "molecules" | `Molecule.atom(SI)`, `Molecule.atom(OX)` | `Si.mol`, `O.mol` → `makemole` |
| box (1500 Si + 3000 O, 0.0664 atoms/Å³, 298 K) | `System("silica", ...).add(silicon, 1500).add(oxygen, 3000)` | `mixato`, `randomise` (no `fmole`: nothing to relax) |
| neutron weights | `NeutronSample("SiO2", ...)` (natural isotopes) | `SiO2.NWTS.dat` → `nwts` |
| EPSR control file | `settings={"qmin": (0.1, 0.0), "rshmin": 1.2}` | `silica.EPSR.inp`, every other keyword at EPSR's default |

Run from the `python/` folder (the `plot` extra gives the figures, `viewer` the 3D image):

```bash
uv run python examples/09_sio2_from_scratch/sio2_from_scratch.py            # ~12 minutes
uv run python examples/09_sio2_from_scratch/sio2_from_scratch.py --quick    # 10 iterations per run
uv run python examples/09_sio2_from_scratch/sio2_from_scratch.py --seed 42  # another starting box
```

With `epsr-py` installed by pip instead (`pip install "./python[plot]"` from the
repository root, see the [top-level README](../../../README.md#1-installation)), run
`python python/examples/09_sio2_from_scratch/sio2_from_scratch.py` from the root.

The schedule is that of `validation/silica_bulk/steps.txt`: `project.equilibrate(iterations=100)`
is the scratch start (1 iteration, which creates the `.pcof` file) plus 99 iterations on
the reference potential, and `project.refine(ereq=50.0, runs=2, iterations=100)` two runs
of 100 iterations with the empirical potential (EReq 50 kJ/mol, as in the Dissolve example).

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the EPSR folder (`.ato`, `.pcof`, `.EPSR.inp`, `.erg`, every `.EPSR.*`
  result), plus `Si.mol`, `O.mol` and `epsr_project.json` (the Python definitions),
  which EPSR does not read.
- `trajectory/`: every iteration, recorded by `epsr.vis` (reload with `vis.Trajectory.load`).
- `figures/`: `01`–`09` and `summary.csv`/`summary.md` from `vis.make_report`; `10` the
  formation of the tetrahedral network (4-coordinated Si, bridging O and the tetrahedral
  order parameter, every 5th iteration); `11` the 3D configuration with Si coloured by
  its number of O neighbours; `12` the reference and empirical potentials; `fits.png`
  the final fit.
- `summary.txt`: the console summary (system, pair parameters, R-factors, peaks, network).

The R-factor falls from 1.0×10⁻³ (reference potential) to 2.3×10⁻⁴ (empirical
potential). The first peaks are Si–O 1.59 Å, O–O 2.55 Å and Si–Si 3.09 Å; at the end
81 % of Si have four O neighbours within 2.0 Å and 89 % of O bridge two Si. The EP
amplitude reaches 24 of the requested 50 kJ/mol, so the script reports that the
refinement has not converged: continue with `project.refine(...)`.

The full run reproduces the validated silica calculation (`validation/silica_bulk`):
every file is byte-identical to the original Fortran EPSR's folder for the same steps
(see the top-level README), and to example 06, which runs it from the existing EPSR input
files. `--seed` gives a different starting box and therefore a different run.
