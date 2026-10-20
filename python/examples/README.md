# Examples

Each example is a folder with one script and a README that describes its inputs,
the calls it makes and its output. Run them from the `python/` folder with
`uv run python examples/<folder>/<script>.py`, or with `epsr-py` installed by pip,
from the repository root with `python python/examples/<folder>/<script>.py` (the
examples read data from the repository, so they need a clone). Each script writes to
`output/` in its own folder, or to `$EPSR_EXAMPLES_WORK/output/` if that variable is set.

## Working with an existing EPSR folder: the DMSO/water test case

`EPSR-test/` (200 DMSO + 400 water, three SANDALS datasets) is a finished EPSR folder.
Each of these examples has the run files of it in its own `template/` folder, copies
them and continues from them. Each takes seconds.

| example | shows |
|---|---|
| [01_single_run](01_single_run) | `Project.from_template`, `project.run`, `project.results`: one run and its fits, g(r) and potentials |
| [02_refinement_loop](02_refinement_loop) | `project.refine` with a progress callback, early stopping and snapshots; the `History` table |
| [03_perturb_potential](03_perturb_potential) | `project.set` for EPSR keywords, in-memory edits of the empirical potential, a run with refinement off |
| [04_tools](04_tools) | `epsr.tools`: the 35 setup and analysis programs, by keyword file, arguments or interactive answers |

## Existing EPSR input files, run within the Python framework

The validated water and silica inputs (`.mol`, `.NWTS.dat`, `.EPSR.inp`, written for
the Fortran EPSR, copied into each example's `template/` folder) are built with EPSR's
setup programs and refined with `Project`.
Each gives the same folder as its from-scratch twin below.

| example | system | full run |
|---|---|---|
| [05_water_existing_epsr_files](05_water_existing_epsr_files) | 1000 H₂O, SANDALS H₂O/D₂O/HDO (same as 07) | ~4 min |
| [06_sio2_existing_epsr_files](06_sio2_existing_epsr_files) | 1500 Si + 3000 O, NIMROD (same as 09) | ~12 min |

## Systems defined from scratch in Python

Every input is defined with the setup API (`AtomType`, `Molecule` or `Framework`,
`System`, `NeutronSample`, `Project.create`); no EPSR file is written by hand. Each
example defines its molecules a different way. Add `--quick` for 10 iterations per run.

| example | system | data (ISIS) | molecules from | full run |
|---|---|---|---|---|
| [07_water_from_scratch](07_water_from_scratch) | 1000 H₂O, 300 K | SANDALS: H₂O, D₂O, HDO 50:50 | `Molecule.from_geometry` | ~4 min |
| [08_cu_btc_from_scratch](08_cu_btc_from_scratch) | Cu-BTC (HKUST-1) unit cell, Bragg peaks | NIMROD | `Framework.from_pdb` | ~4 min |
| [09_sio2_from_scratch](09_sio2_from_scratch) | 1500 Si + 3000 O, amorphous | NIMROD | `Molecule.atom` | ~12 min |
| [10_benzene_from_scratch](10_benzene_from_scratch) | 200 C₆H₆, 300 K | NIMROD: C₆H₆, C₆D₆, 50:50 | `Molecule.from_xyz` | ~10 min |

Water, Cu-BTC and silica are byte-identical to the validated Fortran calculations
(`validation/`); benzene has no Fortran reference.

Tests: `tests/test_examples.py` runs every script but 07 in its short mode (08–10 are
marked `slow`); `tests/test_example_07.py` checks example 07 against its committed output.
