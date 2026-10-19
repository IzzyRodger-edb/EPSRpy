# EPSRpy

**EPSR** (Empirical Potential Structure Refinement, A. K. Soper) turns neutron and X-ray diffraction data into atomistic models of liquids, glasses and disordered materials. It is a standard analysis route for total-scattering data at the ISIS Neutron and Muon Source and other facilities. For 30 years it has been a Fortran program run through an interactive shell, which makes it hard to install, script and reproduce.

This repository contains:

1. **`epsr`, a Python package** (the main focus of development) that drives the Rust engine in-process. It returns NumPy arrays, runs refinement loops and draws figures.
2. **A Rust translation of the whole Fortran EPSR suite**: the `epsr` engine and all 35 setup and analysis programs. It is a drop-in replacement for the Fortran, verified **byte for byte** against the original.
3. **The original Fortran source**, its test case, 115 feature test problems and measured neutron data used for validation.

> **Status**
>
> - ✅ Byte-identical to the Fortran on the supplied test case, 115 feature tests and 3 real-data refinements, at every step
> - ✅ Identical output at any thread count
> - ✅ 5 bugs in the original Fortran found and fixed, including one that silently corrupts three-body energies

Everything is MIT-licensed (see [LICENSE](LICENSE)).

---

## 1. Quick start

Tested on Linux (including WSL2 on Windows).

| tool | why | install (Ubuntu / WSL) |
|---|---|---|
| Rust toolchain (stable, ≥ 1.85) | builds the engine | `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \| sh` |
| C compiler and linker | needed by the Rust and Python builds | `sudo apt install build-essential` |
| Python ≥ 3.10 | the Python package | usually present; `sudo apt install python3` |
| [uv](https://docs.astral.sh/uv/) | Python package management | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| `gfortran` (optional) | only to rebuild the original Fortran and re-run the comparison | `sudo apt install gfortran` |

```bash
# the Python package: refine water against real SANDALS data
cd python && uv sync --all-extras          # compiles the Rust engine on first sync
uv run python examples/05_experimental_water.py --quick
uv run pytest

# the engine on its own, as a drop-in replacement for the Fortran `epsr`
cd ../rust && cargo build --release && cargo test --release
```

## 2. The Python package

[`python/`](python/README.md) is a PyO3 extension built with `uv sync`. It calls the Rust engine's public API unchanged, so its output is still byte-identical to the Fortran:

```python
import epsr
project = epsr.Project.from_template("../EPSR-test", "work/dmso")
sim = project.simulation(overrides={"ntimes": 5})
summary = sim.run()                       # one EPSR run in Rust, GIL released
sim.gr_inter("O1", "H1")                  # partial g(r) as a NumPy array
history = epsr.refine(project, n_runs=20) # refinement loop with convergence history
epsr.tools.run(folder, "nwts", ["sample.NWTS.dat"])   # any of the 35 tools
```

A new system is defined in Python (atom types, molecules, the box, neutron samples with their isotopes, EPSR settings) and built by EPSR's own setup programs, with no input file written by hand:

```python
system = epsr.System("waterbox", temperature=300.0, density=0.1).add(water, 1000)
project = epsr.Project.create("work/water", system, samples, settings={"rshmin": 0.7})
project.equilibrate(iterations=100)                  # scratch start + 99
project.refine(ereq=3.0, runs=2, iterations=100)     # byte-identical to validation/water
```

`epsr.vis` records every iteration and draws convergence plots, F(Q) fits and residual maps, partials, empirical potentials, coordination numbers, bond angles and 3D configurations. Examples 05–07 refine water and silica against real data.

Documentation: [user guide](python/docs/user_guide.md), [API reference](python/docs/api.md), [architecture](python/docs/architecture.md), [testing](python/docs/testing.md).

## 3. The Rust engine

[`rust/`](rust/README.md) translates EPSR version 26 (`fortran_src/epsr_standalone_rev3.f` and everything it calls) into about 35,000 lines of Rust whose only dependency is `rayon` for parallel loops. It covers every optional engine feature: X-ray data, non-cubic boxes, tethered molecules and dihedral groups, container forces, three-body terms, coordination constraints, coarse graining, pressure control, Bragg and diffuse scattering, and lattice refinement. All 35 setup and analysis programs (`makemole`, `mixato`, `fmole`, `nwts`, `partials`, `sdf`, `sharm`, …) are ported too, behind the single `epsr-rs` command. How the engine was parallelised without losing byte-identity is in [PARALLELISM_PLAN.md](PARALLELISM_PLAN.md).

`epsr-rs` takes **the same command-line arguments, reads the same input files and writes the same output files** as the Fortran programs:

```bash
# Fortran:   epsr  <folder> epsr 200DMSO400Water.EPSR.inp
epsr-rs EPSR-test/ epsr 200DMSO400Water.EPSR.inp
epsr-rs EPSR-test/ nwts DMSODWaterD.NWTS.dat      # the setup and analysis programs too
```

Because the output is byte-identical, an existing EPSR project folder (including one behind a published paper) can be continued, re-run or audited with the Rust engine and give exactly the same numbers.

**Fortran bugs found along the way.** Translating line by line exposed real defects in the original ([FORTRAN_BUGS.md](FORTRAN_BUGS.md)). The worst one: after every save, the Fortran overwrites its three-body angles (stored as cosines) with radians, so from the second iteration of any multi-iteration run that uses three-body terms the energy is wrong. In one test case the Fortran reports +1954 kJ/mol per molecule where the correct value is −39. Five such bugs are fixed in the Rust as documented, deliberate divergences. The rest are reproduced on purpose, so that the equivalence still holds, and each has a proposed fix. Inputs on which the Fortran crashes or reads uninitialised memory raise a clear "not supported by the Rust port" error instead.

---

## 4. Validation against the Fortran

**The standard is exact equivalence, not "close agreement".** Monte Carlo trajectories are chaotic: once any accept/reject decision comes out differently, the two runs drift apart for good. Byte-identical output after hundreds of coupled iterations is therefore a very strict test.

![Energy of the water refinement from the Fortran and the Rust over 300 steps, and the divergence caused by moving one molecule by a millionth of an ångström](docs/figures/rust_matches_fortran.png)

On the left, the Fortran and the Rust run the same 300-step water refinement and produce the same energy at every step. Every one of the 39 output files is identical too. On the right, both programs continue for 30 more steps. The Rust stays identical to the Fortran. A copy whose input differs only by moving one of the 1000 molecules by 10⁻⁶ Å matches for 12 steps and then diverges. Nudges of 10⁻⁵ to 10⁻² Å diverge within one or two steps. The run data and the script that draws these figures are in [docs/figures/](docs/figures/make_readme_figures.py).

The original Fortran was rebuilt as a reference (gfortran `-O2 -fno-tree-vectorize`), and the two codes were compared file by file at three levels:

| level | what was compared | result |
|---|---|---|
| **Supplied test case** ([`EPSR-test/`](EPSR-test), 200 DMSO + 400 water, three neutron data sets) | every output file over 12 chained EPSR runs, plus input variants; FNV hashes of every file are pinned in `cargo test` | byte-identical |
| **115 feature tests** ([`test_problems/`](test_problems/README.md)) | one case per feature and per tool; 232 chained passes, each with an archived Fortran reference | byte-identical (bar the 5 deliberate bug fixes) |
| **Real experimental data** ([`validation/`](validation/README.md), see [section 5](#5-refinement-of-measured-data)) | three systems built from scratch and refined for 300 iterations; Fortran, Rust CLI and Python wrapper run side by side, with **every file compared after every step** | all three cases: every step identical |

The parallel Rust also gives identical output at any thread count (`tests/thread_invariance.rs`), which fixes a data race present in the Fortran's OpenMP code.

To reproduce: [rust/docs/VALIDATION.md](rust/docs/VALIDATION.md) and [validation/README.md § How to reproduce](validation/README.md#how-to-reproduce).

---

## 5. Refinement of measured data

The port was developed against the supplied DMSO/water test case. To show that it works on data it was never tuned on, three new systems were built from scratch from the [Dissolve](https://www.projectdissolve.com/) examples and refined against measured neutron data from the ISIS Neutron and Muon Source. The Fortran, the Rust CLI and the Python wrapper were run side by side.

| case | system | data (ISIS instrument) | final R-factor |
|---|---|---|---|
| water | 1000 H₂O, 300 K | SANDALS: H₂O, D₂O, HDO | 2.3 × 10⁻⁴ (from 2.1 × 10⁻² at a random start) |
| silica | 1500 Si + 3000 O, amorphous | NIMROD | 2.3 × 10⁻⁴ |
| Cu-BTC (HKUST-1) | 624-atom MOF unit cell, with Bragg peaks | NIMROD | 1.6 × 10⁻² |

All three codes agree at every step. The overlay plots for [water](validation/water/gr_agreement.png), [silica](validation/silica_bulk/gr_agreement.png) and [Cu-BTC](validation/cu-btc/gr_agreement.png) show zero difference for every pair. The measured data are in [`experimental_data/`](experimental_data/README.md).

The figures below come from the water case: 1000 water molecules fitted to three SANDALS measurements, starting from a random box.

![Measured and simulated scattering signal of heavy water at the start, after 100 steps and after 300 steps](docs/figures/fit_before_after.png)

The first 100 steps use only the reference potential, which gets the general shape right. Switching on the empirical potential closes most of the remaining gap, most visibly around the main peak near Q = 2 Å⁻¹.

![Mismatch between simulation and experiment for the three water samples over 300 refinement steps](docs/figures/mismatch_history.png)

One model has to fit all three isotopic samples at once, and the mismatch for all three falls to about 2 × 10⁻⁴.

![Pair distribution functions of the refined model against the published reference partials of water](docs/figures/structure_vs_published.png)

The refined model shows the O–O peak at 2.79 Å and the O···H hydrogen bond at 1.80 Å, and its pair distribution functions closely follow the published reference partials of Soper (2000). The reference is itself derived from water diffraction data by an EPSR-type refinement, so this is a consistency check rather than an independent measurement. The 300-step schedule is also short, which is why the O–O curve is still noisy between 3 and 5 Å.

### Limitations

* **Some Fortran bugs are reproduced on purpose.** Fixing them would change every output and break byte-identity with the original. They are listed, with proposed fixes, in [FORTRAN_BUGS.md](FORTRAN_BUGS.md).
* **Not everything around the engine is ported.** The PGPLOT/gnuplot plotting front-ends are replaced by `epsr.vis`, and `system_commands.txt` (shell configuration for the interactive EPSRshell) is not read.
* **Installation still needs a compiler.** There are no prebuilt wheels or binaries yet. Testing so far has been on Linux and WSL2.

---

## 6. The theory, in brief

A diffraction experiment on a disordered material measures the total structure factor F(Q). It is a weighted sum of the partial structure factors of every pair of atom types α, β:

$$F(Q) = \sum_{\alpha \le \beta} (2-\delta_{\alpha\beta})\, c_\alpha c_\beta\, b_\alpha b_\beta \,\big[S_{\alpha\beta}(Q) - 1\big]$$

Here c are the atomic fractions and b the scattering lengths (form factors for X-rays). Each S_αβ(Q) is the Fourier transform of a partial pair distribution function g_αβ(r). A single measurement mixes all the partials together. Isotope substitution helps: H₂O, D₂O and HDO have the same structure but different weights b, so measuring several of them gives several independent combinations of the same partials.

EPSR inverts these data with a model:

* **Monte Carlo on a reference potential.** A periodic box of molecules is simulated with a physically sensible potential: Lennard-Jones plus Coulomb, with intramolecular geometry held by harmonic constraints.
* **An empirical potential (EP).** Every few cycles, the simulated F(Q) of each data set is compared with the measurement. The residual is inverted through the weights matrix into a correction for each atom pair and added to the potential as a sum of Poisson functions in r, which pushes the simulation towards the data.
* **A cap on the EP.** The EP's total amplitude is limited by `ereq`, so the data refine the model rather than overwrite the physics.
* **Ensemble averages.** Once the fit stabilises, EPSR accumulates partial g(r), coordination numbers, bond angles and spatial density functions over many configurations.

Fit quality is given by the R-factor, the mean squared difference between data and simulation over the fitted Q range. For water refined against good neutron data it falls to a few × 10⁻⁴.

Key references: Soper, *Chem. Phys.* **202**, 295 (1996); Soper, *Phys. Rev. B* **72**, 104204 (2005); Soper, *EPSRshell user's guide*, RAL-TR-2011-012.

---

## 7. Repository map

| path | contents |
|---|---|
| [`python/`](python/README.md) | the `epsr` Python package (PyO3): bindings, `epsr.vis`, examples, tests, user guide and API reference |
| [`rust/`](rust/README.md) | the Rust engine (`epsr-core`), the drop-in CLI `epsr-rs`, the Fortran reference build and comparison tools, and the architecture, Fortran-map and validation docs |
| [`fortran_src/`](fortran_src) | the original Fortran, untouched |
| [`EPSR-test/`](EPSR-test) | the original 200 DMSO / 400 water test case, the base of the test problems |
| [`test_problems/`](test_problems/README.md) | the 115 `EPSR-test-ex*` verification cases with archived Fortran references |
| [`experimental_data/`](experimental_data/README.md) | measured neutron data (water, silica, Cu-BTC, benzene), their Dissolve inputs and reference partials |
| [`validation/`](validation/README.md) | water, silica and Cu-BTC refined against measured data by the Fortran, the Rust CLI and the Python wrapper |
| [`visualisation/`](visualisation/README.md) | the per-iteration Rust recorder and the standalone `epsrvis` toolkit that `epsr.vis` grew from |
| [`docs/figures/`](docs/figures/make_readme_figures.py) | the figures in this README, the script that draws them and the run data they use |
| [`archive/`](archive/README.md) | retired components: [Mellor](archive/mellor/ARCHIVED.md), an LLM agent for EPSR refinements, kept so it can be reinitialised |
| [FORTRAN_BUGS.md](FORTRAN_BUGS.md) | bugs found in the Fortran and their status in the port |
| [PARALLELISM_PLAN.md](PARALLELISM_PLAN.md) | how the engine was parallelised without losing byte-identity |

### Acknowledgements

EPSR is the work of Alan Soper (ISIS Neutron and Muon Source, STFC). The example systems and data come from the [Dissolve](https://www.projectdissolve.com/) project (T. G. A. Youngs). The project started at Reinvent the Wheel 2026, hosted by Science Works and Inherent Laboratories.
