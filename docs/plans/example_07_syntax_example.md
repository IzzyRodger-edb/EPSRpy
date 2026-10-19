# Wireframe: example 07 in the proposed syntax

**Status:** design wireframe, 2026-10-19. **None of the new API below exists yet.** This
is the target that [python_setup_api.md](python_setup_api.md) is implemented towards.
When it is built, this script replaces
`python/examples/07_water_from_scratch/water_from_scratch.py`.

**Key**
- 🆕 marks a call that is proposed and doesn't exist yet. The § number points to the plan.
- Unmarked calls exist today (`epsr.plot`, `epsr.vis`, `Simulation` accessors).

**Decisions reflected here** (plan, top table):

| | decision | in the script |
|---|---|---|
| D1 | EPSR keyword names for settings | `settings={"rshmin": 0.7}`, `ereq=3.0` |
| D2 | scratch start inside the first `equilibrate` | `project.equilibrate(iterations=100)` = 1 + 99 |
| D3 | `seed=` option; `None` = EPSR's default (seed 1) | `Project.create(..., seed=options.seed)` |
| D4 | atoms named within a molecule; bonds/angles by name or 0-based index | `("O", OW, …)`, `bonds=[("O", "H1"), …]` |
| D5 | `exchangeable` on the `AtomType` | `AtomType("HW", …, exchangeable=True)` |
| D6 | number density in atoms/Å³ only | `System(..., density=0.1)` |
| D7 | isotopes as mass numbers only, `{mass_number: fraction}`; 0 = natural | `isotopes={"HW": {1: 0.5, 2: 0.5}}` |

**Acceptance.** The full run (without `--quick`) must write an `epsr_run/` folder
byte-identical to the committed `output/epsr_run/`, which is the validated
`validation/water` calculation. As today, the date line of `.EPSR.out` and the `.qdr`
file are excepted.

---

## 1. Shape of the script

| step | today (`water_from_scratch.py`) | proposed |
|---|---|---|
| define the molecule | dicts + `mol_file_text()` renders the `.mol` by hand | `AtomType`, `Molecule` 🆕 |
| define the box | `N_MOLECULES`, `ATOMIC_DENSITY` constants | `System(...).add(water, 1000)` 🆕 |
| define the data | `DATASETS` tuples with abundance strings (`"1 0.5 2 0.5"`) | `NeutronSample(isotopes={"HW": {1: 0.5, 2: 0.5}})` 🆕 |
| build the folder | copy data, write `Water.mol`, `makemole`, `mixato` (prompt answers), `randomise`, delete `Water.mol`, `fmole`, `write_keyword_file` + `nwts` ×3, write `.EPSR.inp` by hand | `Project.create(...)` 🆕 |
| scratch start + equilibrate | `project.simulation()` + `run()`, then `overrides={"niter": 99}` | `project.equilibrate(iterations=100)` 🆕 (1 scratch + 99) |
| refine with the EP | `refine(..., initial_overrides={"ereq": 3.0})`, or `vis.refine_recorded` when plotting is installed | `project.refine(ereq=3.0, runs=2, iterations=100, record=...)` 🆕 |
| results | `project.simulation()`; Lorentz–Berthelot and peak windows by hand | `project.results()` 🆕, `pair_params()` 🆕, `fit("D2O")` 🆕, `first_peak` 🆕 |
| figures | `epsr.plot` + `epsr.vis` | unchanged, fed from `project.trajectory()` 🆕 |

The setup part goes from about 230 lines to about 60, and no EPSR file text appears in
the script.

---

## 2. The script

```python
"""Refine liquid water against neutron data, with every input defined in Python.

Run from the ``python/`` folder:

    uv run python examples/07_water_from_scratch/water_from_scratch.py            # ~4 minutes
    uv run python examples/07_water_from_scratch/water_from_scratch.py --quick    # ~30 s
    uv run python examples/07_water_from_scratch/water_from_scratch.py --seed 42  # another starting box

1000 water molecules are built from scratch and fitted to three SANDALS measurements
(H2O, D2O and HDO 50:50). The schedule is that of ``validation/water/steps.txt``:
equilibrate for 100 iterations on the reference potential, then two runs of 100
iterations with the empirical potential (EReq 3 kJ/mol).
"""

import argparse
import os
from pathlib import Path

import epsr
from epsr.samples import NeutronSample                     # 🆕 §4.3
from epsr.system import AtomType, Molecule, System         # 🆕 §4.3

EXAMPLE_FOLDER = Path(__file__).resolve().parent
REPOSITORY = EXAMPLE_FOLDER.parents[2]
DATA = REPOSITORY / "experimental_data" / "water" / "data"
REFERENCE_PARTIALS = REPOSITORY / "experimental_data" / "water" / "reference_partials"
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"

parser = argparse.ArgumentParser()
parser.add_argument("--quick", action="store_true", help="10 iterations per run")
parser.add_argument("--seed", type=int, default=None, help="random seed (default: EPSR's own, = 1)")
parser.add_argument("--no-3d", action="store_true")
options = parser.parse_args()
EQUILIBRATION = 11 if options.quick else 100   # includes the 1-iteration scratch start
REFINEMENT = 10 if options.quick else 100

# --- 1. The molecule ----------------------------------------------------------------
# Reference potential per atom type: Lennard-Jones epsilon (kJ/mol), sigma (Å) and charge (e).
# EPSR always derives each pair from its two types: eps_ij = sqrt(eps_i eps_j),
# sigma_ij = (sigma_i + sigma_j)/2, plus Coulomb 1390 q_i q_j / r. Per-pair LJ is not possible.
# HW has epsilon = 0, i.e. no LJ site (SPC convention), so OW-HW and HW-HW are Coulomb-only;
# the hydrogen bond comes from the O-H Coulomb attraction (and, once ereq > 0, the EP).
# system.pair_parameters() below prints this, with a "terms" column per pair.
OW = AtomType("OW", element="O", epsilon=0.6503, sigma=3.165492, charge=-0.82,
              coarse_grain={"rmax": 1.55})          # as in the validated .ato
HW = AtomType("HW", element="H", epsilon=0.0, sigma=0.0, charge=0.41, mass=2.0,
              exchangeable=True)                    # labile H: exchanges between molecules in every sample

water = Molecule.from_geometry(                     # 🆕 builds coordinates from internal geometry
    "Water",
    centre=("O", OW),                               # atoms are named within the molecule
    arms=[("H1", HW, 0.976), ("H2", HW, 0.976)],    # bond lengths (Å)
    angle=107.134,                                  # H1-O-H2 (degrees)
    vibtemp=65, angtemp=1, dihtemp=10,              # restraint stiffness ("temperatures")
)
# Intramolecular: bonds and angles are harmonic *restraints* to this geometry (angles as
# 1-3 distances). There are no per-bond force constants: stiffness is vibtemp (x angtemp
# for angles), scaled by sqrt(reduced mass)/r0 and the temperature. See plan §1.8.
# Equivalent explicit form (atoms named; bonds/angles by name or 0-based index):
#   Molecule("Water",
#            atoms=[("O", OW, (0, 0, 0)), ("H1", HW, (0.976, 0, 0)), ("H2", HW, (-0.28754, 0.93268, 0))],
#            bonds=[("O", "H1"), ("O", "H2")], angles=[("H1", "O", "H2")])
# or from a file: Molecule.from_xyz("water.xyz", types={"O": OW, "H": HW})
#   (atoms then named element + running number: "O1", "H1", "H2")

# --- 2. The box ---------------------------------------------------------------------
system = System("waterbox", temperature=300.0, density=0.1)    # atoms/Å³, cubic box
system.add(water, 1000)
print(system)                       # 🆕 1000 × Water, 3000 atoms, box 31.07 Å, types OW, HW
print(system.pair_parameters())     # 🆕 one row per pair, with a "terms" column:
#   pair    eps_ij  sigma_ij  1390*qi*qj  terms
#   OW-OW   0.6503  3.165     +934.6      LJ + Coulomb
#   OW-HW   0       1.583     -467.3      Coulomb only
#   HW-HW   0       0         +233.7      Coulomb only

# --- 3. The data --------------------------------------------------------------------
# isotopes: {atom type: {mass number: fraction}}, exactly EPSR's `abundances` pairs (D7).
# Mass number 0 = the natural element; types left out are natural ({0: 1.0}).
# Fractions must sum to 1, at most 5 isotopes per type, and each mass number is checked
# against the engine's Sears (1991) scattering-length table when the sample is created.
# Exchange comes from AtomType(exchangeable=True) (D5).
samples = [
    NeutronSample("D2O", DATA / "SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}),
    NeutronSample("HDO", DATA / "SLS18500-HDO5050.mint01", isotopes={"HW": {1: 0.5, 2: 0.5}}),
    NeutronSample("H2O", DATA / "SLS18498-H2O.mint01"),        # natural abundance
]

# --- 4. Build the EPSR project --------------------------------------------------------
project = epsr.Project.create(                     # 🆕 §4.3; does not run EPSR yet
    OUTPUT / "epsr_run",
    system,
    samples,
    settings={"rshmin": 0.7},                      # EPSR keyword names; checked against the keyword table
    relax_iterations=100,                          # fmole
    seed=options.seed,                             # None: EPSR's default (= seed 1), byte-identical
    overwrite=True,
)
print(project.status())    # 🆕 "not run yet; 1000 molecules; 3 datasets (D2O, HDO, H2O); seed default"

# --- 5. Refine ------------------------------------------------------------------------
record = epsr.vis.available()                      # 🆕 True when the 'plot' extra is installed
equilibration = project.equilibrate(iterations=EQUILIBRATION, record=record)
#   🆕 first call: scratch start (1 iteration, creates the .pcof) + 99 iterations, ereq = 0
history = project.refine(ereq=3.0, runs=2, iterations=REFINEMENT, record=record)
#   🆕 ereq switched on for the first run; the second run continues from the files

# --- 6. Results -----------------------------------------------------------------------
result = project.results()                         # 🆕 read-only Simulation of the finished folder
print(f"R-factor: {equilibration.r_factor[-1]:.3e} on the reference potential")
print(f"          {history.r_factor[-1]:.3e} with the empirical potential")
for sample in project.samples:                     # 🆕 names from the project definition
    fit = result.fit(sample.name)                  # 🆕 select a dataset by name
    print(f"  {sample.name}: fitted over Q = {fit.q[0]:.2f}-{fit.q[-1]:.2f} 1/Å, "
          f"R = {result.summary.r_factor_for(sample.name):.2e}")   # 🆕
for a, b, window, meaning in (("OW", "OW", (2.0, 4.0), "nearest-neighbour waters"),
                              ("OW", "HW", (1.5, 2.3), "hydrogen bond")):
    peak = epsr.analysis.first_peak(result.r, result.gr_inter(a, b), window)   # 🆕
    print(f"{a}-{b} first peak at {peak:.2f} Å ({meaning})")
print(result.potential.pair_params())              # 🆕 mixed parameters as EPSR uses them

if not record:
    raise SystemExit("(install the 'plot' extra for the figures)")

# --- 7. Figures (existing epsr.plot / epsr.vis API) -------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from epsr import vis

figures = OUTPUT / "figures"
trajectory = project.trajectory()                  # 🆕 every iteration recorded with record=True (both stages)
trajectory.save(OUTPUT / "trajectory")
vis.make_report(trajectory, figures, verbose=False)            # 01-09 + summary tables

axes = epsr.plot.fits(result)
axes.figure.savefig(figures / "fits.png", dpi=150)

figure, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
epsr.plot.potentials(result, which="reference", ax=left)
epsr.plot.potentials(result, which="empirical", ax=right)
figure.savefig(figures / "12_reference_and_empirical_potential.png", dpi=150)

reference = vis.load_reference_partials({
    "OW-OW": REFERENCE_PARTIALS / "h2o298aoogor.txt",
    "OW-HW": REFERENCE_PARTIALS / "h2o298aohgor.txt",
    "HW-HW": REFERENCE_PARTIALS / "h2o298ahhgor.txt",
})
vis.plots.partials_vs_reference(trajectory, reference, source="Soper (2000), 298 K") \
    .savefig(figures / "10_partials_vs_experiment.png", dpi=150)

if not options.no_3d:
    from epsr.vis.viewer import ConfigurationViewer
    ConfigurationViewer(trajectory, gr="total", slab=10.0, reference=reference, off_screen=True) \
        .screenshot(figures / "11_configuration.png", step=-1)
```

---

## 3. Notes for the implementer

1. **Order of `system.add` calls = atom-type order = pair order.** Here OW comes before
   HW because water is the only molecule and OW is its first atom. Document this in
   `System.add`.
2. **`Molecule.from_geometry`** is a convenience for small molecules with a central atom
   (water, NH₃, CH₄). It is not essential. The explicit `atoms=` / `bonds=` / `angles=`
   form is the core API. If it's added, it must give the same coordinates as example 07
   today: O at the origin, H₁ on +x, H₂ in the xy plane at the H-O-H angle, rounded to
   float32. Otherwise the box, and with it byte-identity, changes.
3. **`coarse_grain={"rmax": 1.55}` on OW** is there only to match the validated `.ato`.
   That value comes from the EPSR-test `.mol` files. Consider making it the default for
   O (it is EPSR's convention), or keep it explicit. Decide while implementing Phase 2,
   by checking what makemole writes when it is omitted (0.5·σ).
4. **`settings={"rshmin": 0.7}`.** `rshmin` is a per-dataset keyword. Given globally, it
   applies to every sample unless a `NeutronSample(rshmin=...)` overrides it. `qmin` is
   left at its default (0.05), as in the validated input.
5. **`--quick` uses `EQUILIBRATION = 11`**, i.e. 1 scratch + 10, matching today's quick
   mode (scratch + 10 + 2 × 10).
6. **`equilibration` is a `History`** with two rows on the first call (scratch, then 99).
   `equilibration.r_factor[-1]` is the R-factor at the end of equilibration, as
   `equilibrated.r_factor` is today.
7. **`epsr.vis.available()`, `Summary.r_factor_for(name)` and `project.trajectory()`**
   are small conveniences proposed here. `project.trajectory()` returns one
   `Trajectory` concatenating every run recorded in the project, so equilibration and
   refinement appear together, as in today's single recorder. Add all three to Phase 5.
8. **The hydrogen-bond count** of today's script (`vis.metrics.neighbour_pairs`) is left
   out of the wireframe for brevity. It keeps working unchanged on `trajectory`.
9. **Isotopes (D7)** render one-to-one into the `.NWTS.dat`:
   - `{1: 0.5, 2: 0.5}` → `abundances 1 0.5 2 0.5`
   - an omitted type → `abundances 0 1`

   Keep the dict order as given: `nwts` doesn't sort, and identical text helps
   byte-identity. Validate against `_core.isotope_table()` (plan §4.5 item 2a); an
   unknown mass number would otherwise make `nwts` wait for typed input.
10. **Atom names (D4)** are Python-side only. The generated `.mol` numbers atoms
    `001…` in list order, so names never reach EPSR files and can't affect byte-identity.
    They must be unique within the molecule and are stored in `epsr_project.json`.
11. **Labels:** `OW` and `HW` meet the 3-character limit on atom-type labels (plan,
    validation rules).
