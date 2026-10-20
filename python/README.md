# epsr-py: Python interface to the Rust EPSR engine

`epsr` runs and inspects EPSR (Empirical Potential Structure Refinement, A. K. Soper)
from Python. The engine is the Rust translation of the original Fortran in
[`../rust`](../rust), used unchanged through its public API; it reproduces the
Fortran output byte for byte.

## Install

With pip, from the repository root (you need a Rust toolchain, stable 1.85 or later, and
Python 3.10 or later; the first build compiles the engine, 5–10 minutes):

```bash
pip install "./python[plot]"
# or without a clone:
pip install "epsr-py[plot] @ git+https://github.com/IzzyRodger-edb/EPSRpy#subdirectory=python"
```

The extras are `plot` (matplotlib, scipy), `viewer` (pyvista, for 3D views and
movies), `pandas`, `progress` (tqdm) and `all`. Pre-built wheels, which need no Rust,
are described in the [top-level README](../README.md#1-installation).

For development, use [uv](https://docs.astral.sh/uv/):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh        # if uv is not installed
export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
cd python
uv sync --all-extras        # creates .venv, compiles the Rust engine, installs epsr
```

`uv sync` rebuilds the extension whenever the Rust engine or the bindings change.

## Run

One EPSR run on a copy of a simulation folder:

```python
import epsr

project = epsr.Project.from_template("../EPSR-test", "work/dmso")  # copy: EPSR rewrites its inputs
simulation = project.simulation(overrides={"ntimes": 5})           # read .EPSR.inp, .ato, .pcof, data
summary = simulation.run()                                         # one EPSR run (Rust, GIL released)

summary.r_factor                           # fit quality
simulation.gr_inter("O1", "H1")            # partial g(r) as a float32 numpy array
simulation.potential.empirical()           # (n_pairs, nr) empirical potential, kJ/mol
history = epsr.refine(project, n_runs=20)  # repeated runs with a convergence history
```

Every run overwrites `.ato`, `.pcof`, `.EPSR.inp` and `.EPSR.inpa` and appends to
`.erg`; that is how EPSR carries its state from run to run. The 35 setup and
analysis programs run as `epsr.tools.run(folder, "nwts", ["sample.NWTS.dat"])`.

## Set up a system from Python

A new simulation is defined in Python and built by EPSR's own setup programs; no
EPSR input file is written by hand
([examples/07_water_from_scratch](examples/07_water_from_scratch/water_from_scratch.py)):

```python
from epsr import AtomType, Molecule, NeutronSample, Project, System

OW = AtomType("OW", element="O", epsilon=0.6503, sigma=3.165492, charge=-0.82,
              coarse_grain={"rmax": 1.55})
HW = AtomType("HW", element="H", epsilon=0.0, sigma=0.0, charge=0.41, mass=2.0,
              exchangeable=True)
water = Molecule.from_geometry("Water", centre=("O", OW),
                               arms=[("H1", HW, 0.976), ("H2", HW, 0.976)], angle=107.134)
system = System("waterbox", temperature=300.0, density=0.1).add(water, 1000)
samples = [NeutronSample("D2O", "SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}),
           NeutronSample("H2O", "SLS18498-H2O.mint01")]

project = Project.create("work/water", system, samples, settings={"rshmin": 0.7})
project.equilibrate(iterations=100)              # scratch start + 99 on the reference potential
project.refine(ereq=3.0, runs=2, iterations=100) # empirical potential on
result = project.results()                       # result.fit("D2O"), result.gr_inter("OW", "HW")
```

With the three SANDALS samples, this is byte for byte the validated water
calculation (`uv run pytest -m slow -k example_07`).

| command | purpose |
|---|---|
| `uv run python examples/07_water_from_scratch/water_from_scratch.py --quick` | refine water against real data (about 1 minute) |
| `uv run python examples/09_sio2_from_scratch/sio2_from_scratch.py --quick` | refine amorphous silica against real data (about 1.5 minutes) |
| `uv run pytest -n auto` | the whole test suite (365 tests, including the slow refinements on real data) |
| `uv run pytest -m "not slow" -n 4` | the quick tests (222 tests, about 3 minutes on 4 cores) |
| `uv run isort . && uv run black . && uv run ruff check . && uv run mypy` | format, lint, type-check |

## Examples on real experimental data

Four systems from [`../experimental_data`](../experimental_data) are built from scratch
in Python and refined against measured neutron data. Each example is a folder with one
script and a README; see [examples/README.md](examples/README.md) for the full list.

| example | system | data (ISIS) | how the system is defined |
|---|---|---|---|
| [07_water_from_scratch](examples/07_water_from_scratch) | 1000 H₂O, 300 K | SANDALS: H₂O, D₂O, HDO | `Molecule.from_geometry` |
| [08_cu_btc_from_scratch](examples/08_cu_btc_from_scratch) | Cu-BTC (HKUST-1) unit cell, Bragg peaks | NIMROD | `Framework.from_pdb` |
| [09_sio2_from_scratch](examples/09_sio2_from_scratch) | 1500 Si + 3000 O, amorphous | NIMROD | `Molecule.atom` (one-atom "molecules") |
| [10_benzene_from_scratch](examples/10_benzene_from_scratch) | 200 C₆H₆, 300 K | NIMROD: C₆H₆, C₆D₆, 50:50 | `Molecule.from_xyz` |

Examples [05](examples/05_water_existing_epsr_files) and
[06](examples/06_sio2_existing_epsr_files) run the water and silica calculations
from existing EPSR input files (`.mol`, `.NWTS.dat`, `.EPSR.inp` written for the
Fortran EPSR) inside the Python framework, and give the same folders as 07 and 09.

Water, silica and Cu-BTC have been run by the original Fortran EPSR, by the Rust CLI
and by this wrapper, and every output file of the three is byte-identical; the plots
below overlay their partial g(r). Benzene has no Fortran reference yet.

```bash
uv run python examples/07_water_from_scratch/water_from_scratch.py            # full schedule, about 4 minutes
uv run python examples/07_water_from_scratch/water_from_scratch.py --quick    # 10 iterations per run
```

After the full water schedule the R-factor falls from 6.1×10⁻⁴ (reference potential)
to 2.3×10⁻⁴, with the O···H hydrogen-bond peak at 1.80 Å and the O–O peak at 2.79 Å.

![Water: partial g(r) from the Fortran, the Rust CLI and the Python wrapper](../validation/water/gr_agreement.png)

![Silica: partial g(r) from the Fortran, the Rust CLI and the Python wrapper](../validation/silica_bulk/gr_agreement.png)

![Cu-BTC: partial g(r) from the Fortran, the Rust CLI and the Python wrapper](../validation/cu-btc/gr_agreement.png)

### Figures of the whole refinement

The examples on real data record every EPSR iteration with `epsr.vis` and end with
figures of the refinement in `examples/<example>/output/figures/`:
- convergence (energy, pressure, R-factor and R² per data set, EP amplitude, MC acceptance);
- F(Q) fits and residual maps, partials, the empirical potential, coordination numbers and bond angles;
- the 3D configuration next to the partial g(r);
- for water, the partial g(r) against the experimental partials of Soper (2000);
- for silica, the formation of the tetrahedral network;
- for benzene, the relative orientation of neighbouring rings.

Add `--no-3d` to skip the 3D image, e.g. on a machine without OpenGL.

## More

* [docs/user_guide.md](docs/user_guide.md): setting up a system from Python, the workflow, overrides and settings, perturbing the potential, tools, plotting and `epsr.vis`, parallel runs.
* [docs/api.md](docs/api.md): every public class and function.
* [docs/architecture.md](docs/architecture.md) and [docs/testing.md](docs/testing.md): how the wrapper is built and verified.
* [examples/](examples/README.md): `01` to `04` show a single run, a refinement loop, potential perturbation and the tools on the DMSO/water test case; `05` and `06` run existing EPSR input files; `07` to `10` set up water, Cu-BTC, silica and benzene entirely in Python.
