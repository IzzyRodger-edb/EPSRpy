# internal_play: trying out EPSRpy

A scratch folder for installing EPSRpy and experimenting with it. It holds a copy of
example 07 (liquid water built from scratch and refined against neutron data) that
reads its data from this folder and writes everything to `internal_play/output/`, so
you can edit it freely without touching the real examples.

```
internal_play/
├── README.md                  this file
├── water_from_scratch.py      the script to experiment with (copy of example 07)
├── data/                      SANDALS neutron data: H2O, D2O, HDO 50:50 (.mint01)
├── reference_partials/        Soper (2000) partial g(r) of water, for comparison figures
└── output/                    created by the script (ignored by git)
```

These instructions are for Linux (including WSL2 on Windows). The top-level
[README](../README.md#1-installation) has the details for other distributions.

## 1. System packages (once)

Installing from a git clone compiles the Rust engine, so you need Python ≥ 3.10, git,
a C compiler and a Rust toolchain ≥ 1.85.

On Debian / Ubuntu / WSL2:

```bash
sudo apt update
sudo apt install python3 python3-venv python3-pip git build-essential curl
```

On Fedora / RHEL: `sudo dnf install python3 python3-pip git gcc curl`.

Install Rust with rustup (no `sudo`; distribution packages are often too old):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
```

Check: `python3 --version` should say 3.10 or later, and `rustc --version` 1.85 or later.

## 2. Get the code with git

```bash
git clone https://github.com/IzzyRodger-edb/EPSRpy
cd EPSRpy
```

To pick up later changes: `git pull` from inside `EPSRpy/`, then repeat the `pip install`
line in step 3 so the engine is rebuilt.

## 3. Install EPSRpy into a virtual environment

From the `EPSRpy/` folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install "./python[plot]"
```

The first install compiles the engine and takes about 5–10 minutes. `[plot]` adds
matplotlib and SciPy for the figures. Use `"./python[plot,viewer]"` if you also want the
3D configuration image. That needs the OpenGL libraries
(`sudo apt install libgl1 libegl1 libx11-6 libxcursor1`).

Check the installation:

```bash
python -c "import epsr; print(epsr.__version__, epsr.VERSION_DATE)"
```

This should print `0.1.0 26: 2017-03-10`.

In every new terminal, activate the environment again before using EPSRpy:

```bash
source ~/path/to/EPSRpy/.venv/bin/activate
```

## 4. Run the water example

```bash
cd internal_play
python water_from_scratch.py --quick    # 10 iterations per run, about 1–2 minutes
python water_from_scratch.py            # the full refinement, about 4 minutes or more
python water_from_scratch.py --seed 42  # a different random starting box
```

Add `--no-3d` to skip the 3D image (it is skipped anyway without the `viewer` extra).
`mixato` prints `sh: 1: ./Water.ato: Permission denied` during setup. This is harmless.

The script prints the system, the pair potential parameters, the R-factors and the
first-peak positions, and writes:

- `output/epsr_run/`: the EPSR simulation folder (all the usual `.ato`, `.pcof`,
  `.EPSR.inp`, `.EPSR.*` files)
- `output/figures/`: convergence, F(Q) fits, partial g(r), potentials, comparison with
  Soper's partials, etc.
- `output/trajectory/`: every iteration, for later analysis
- `output/summary.txt`: the console summary

Each run overwrites `output/`. Rename it first if you want to keep a run.

## 5. Importing the libraries yourself

With the environment active, everything is in the `epsr` package. The script uses:

```python
import epsr                                           # Project, analysis, plot, ...
from epsr.system import AtomType, Molecule, System    # build the box
from epsr.samples import NeutronSample                # the measured datasets
from epsr import vis                                  # refinement analysis and figures ([plot])
from epsr.vis.viewer import ConfigurationViewer       # 3D views ([viewer] only)
```

The script is organised in numbered sections, which are good places to experiment:

1. **The molecule**: `AtomType` (ε, σ, charge, isotopes) and `Molecule.from_geometry`
2. **The box**: `System(..., temperature=, density=)` and `system.add(water, 1000)`
3. **The data**: one `NeutronSample` per measurement, with its isotope substitution
4. **The EPSR project**: `epsr.Project.create(...)`, with EPSR keywords in `settings=`
5. **Refinement**: `project.equilibrate(...)` then `project.refine(ereq=..., runs=..., iterations=...)`
6. **Results**: `project.results()`, fits, partial g(r), peak positions
7. **Figures**

To look at a finished run interactively (e.g. in `python` or Jupyter, started in
`internal_play/`):

```python
import epsr

project = epsr.Project("output/epsr_run")
print(project.status())                     # iterations, ereq, R-factor, datasets
result = project.results()

for name in project.dataset_names:
    fit = result.fit(name)                  # fit.q and the data / model F(Q)
    print(name, fit.q[0], fit.q[-1])

g_oo = result.gr_inter("OW", "OW")          # intermolecular partial g(r) on result.r
print(epsr.analysis.first_peak(result.r, g_oo, (2.0, 4.0)))   # ≈ 2.8 Å

ax = epsr.plot.fits(result)                 # matplotlib figure of every fit
ax.figure.savefig("fits.png")
```

Use `help(epsr.Project)`, `help(epsr.system.Molecule)`, etc. for the full documentation
of each class. The other examples are in [`python/examples/`](../python/examples/), and
`epsrvis --help` lists the command-line tool for figures and movies.

## Troubleshooting

| problem | fix |
|---|---|
| `error: externally-managed-environment` | the virtual environment is not active: `source .venv/bin/activate` |
| `ModuleNotFoundError: No module named 'epsr'` | same: activate the environment the package was installed into |
| `cargo` / `rustc` not found, or "requires rustc 1.85" | install or update Rust with rustup (step 1), then open a new terminal |
| `linker 'cc' not found` | install the C compiler: `sudo apt install build-essential` |
| `ensurepip is not available` | `sudo apt install python3-venv` |
| no figures, "install the 'plot' extra" | `pip install "./python[plot]"` from `EPSRpy/` |
