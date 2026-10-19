# User guide

This guide walks through the EPSR workflow in Python. The examples use the
DMSO/water test case in `EPSR-test/` (200 D6-DMSO + 400 water molecules, three
neutron datasets). Atom types `S`, `C`, `O`, `H` belong to DMSO and `O1`, `H1`
to water.

## 1. How EPSR works, in one paragraph

EPSR builds an atomistic model of a liquid or amorphous material that agrees with
measured diffraction data. A Monte Carlo simulation runs with a *reference
potential* (Lennard-Jones plus Coulomb, defined in the `.ato` file). After every
`ntimes` Monte Carlo cycles the simulated structure factors are compared with the
data, and the difference is used to refine an *empirical potential* (EP) that
pushes the simulation towards the data. The EP is limited in size by the
parameter `ereq`. Repeating this many times gives a model whose partial pair
distribution functions g(r) can be analysed.

A simulation folder holds:

| file | contents |
|---|---|
| `<name>.EPSR.inp` | run parameters (`ntimes`, `niter`, data sets, ...) |
| `<name>.ato` | the configuration, atom types, reference-potential parameters, random-number state |
| `<name>.pcof` | EP control parameters (`ereq`, ...) and the EP coefficients of every atom pair |
| `*.wts`, data files | neutron or X-ray weights and the measured data |
| `<name>.EPSR.*` | results written by each run (`.g01` partial g(r), `.u01` simulated F(Q), `.erg` history, ...) |

## 2. Setting up a system from Python

A new simulation needs a molecule template, a box, the neutron weights of every
sample and a run-control file. `epsr.Project.create` builds all of them from Python
objects, by writing the files EPSR's own setup programs read (`makemole`, `mixato`,
`randomise`, `fmole`, `nwts`) and running those programs, so the folder is exactly
what EPSR would build from the same numbers. You never write an EPSR input file.
The complete example is `examples/07_water_from_scratch/water_from_scratch.py`.

```python
import epsr
from epsr import AtomType, Molecule, NeutronSample, System

OW = AtomType("OW", element="O", epsilon=0.6503, sigma=3.165492, charge=-0.82,
              coarse_grain={"rmax": 1.55})
HW = AtomType("HW", element="H", epsilon=0.0, sigma=0.0, charge=0.41, mass=2.0,
              exchangeable=True)
water = Molecule.from_geometry("Water", centre=("O", OW),
                               arms=[("H1", HW, 0.976), ("H2", HW, 0.976)],
                               angle=107.134)

system = System("waterbox", temperature=300.0, density=0.1)   # atoms/Å³, cubic box
system.add(water, 1000)
print(system.pair_parameters())

samples = [
    NeutronSample("D2O", "data/SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}),
    NeutronSample("HDO", "data/SLS18500-HDO5050.mint01", isotopes={"HW": {1: 0.5, 2: 0.5}}),
    NeutronSample("H2O", "data/SLS18498-H2O.mint01"),          # natural isotopes
]

project = epsr.Project.create("work/water", system, samples,
                              settings={"rshmin": 0.7}, seed=None, overwrite=True)
project.equilibrate(iterations=100)                 # scratch start + 99, ereq = 0
history = project.refine(ereq=3.0, runs=2, iterations=100)
result = project.results()                          # read-only Simulation
result.fit("D2O"), result.gr_inter("OW", "HW"), result.potential.pair_params()
```

**Atom types** (`AtomType`) carry the reference potential: Lennard-Jones ε (kJ/mol)
and σ (Å), charge (e), mass (default: EPSR's atomic weight of the element) and
whether the atom is `exchangeable` (labile hydrogen). Labels have at most three
characters. EPSR always mixes the pair potentials from the two types
(Lorentz–Berthelot, ε_ij = √(ε_i ε_j), σ_ij = (σ_i + σ_j)/2, plus Coulomb
1390 q_i q_j / r); there is no per-pair Lennard-Jones input, and a type with ε = 0
has no Lennard-Jones site at all, so its pairs are Coulomb-only.
`system.pair_parameters()` prints the mixed values with a `terms` column.

**Molecules** (`Molecule`) are named atoms with positions, bonds, angles and
rotations; atoms are given by name or 0-based index:

```python
Molecule("Water",
         atoms=[("O", OW, (0, 0, 0)), ("H1", HW, (0.976, 0, 0)), ("H2", HW, (x2, y2, 0))],
         bonds=[("O", "H1"), ("O", "H2")],
         angles=None)                     # None: every angle makemole needs, from the geometry
Molecule.from_geometry(...)               # a central atom with 1-4 arms (water, NH3, CH4)
Molecule.from_xyz("water.xyz", types={"O": OW, "H": HW})   # bonds="auto" by covalent radii
Molecule.from_mol("DMSO.mol")             # an EPSR .mol file
Molecule.from_ato("DMSO.ato")             # a template .ato (restraints rounded to 3 digits)
Molecule.atom(ion_type)                   # a one-atom "molecule" (ions, silica)
```

Bonds and angles are EPSR's harmonic *restraints* to the target geometry, not a
force field: there are no per-bond force constants. Their stiffness is set per
molecule by `vibtemp`, `angtemp` and `dihtemp` and scaled by the masses, the bond
length and the temperature. Angles are restrained as 1-3 distances. `bond_lengths=`
and `angle_values=` set targets that differ from the coordinates. Dihedral
restraints are not supported yet.

**The box** (`System`) takes the temperature and the atomic number density in
atoms/Å³. **The order of the `add` calls, and of the atoms in each molecule, fixes
the atom-type order**, and with it the order of every pair-resolved result.

**Samples** (`NeutronSample`) are one data file each (a Gudrun `.mint01` file, checked
when the sample is created) with the isotopic composition as
`{atom type: {mass number: fraction}}`, mass number 0 being the natural element,
exactly EPSR's `abundances`. Types that are not given are natural. The mass numbers
are checked against the engine's scattering-length table before anything runs.
`exchangeable=[...]` on a sample replaces the types' own flags for that sample. The
sample name selects the dataset everywhere: `result.fit("D2O")`,
`result.summary.r_factor_for("D2O")`, `settings={"data": {"D2O": {"tweak": 1.05}}}`.

**Settings** use EPSR's own keyword names (`ntimes`, `qmin`, `rshmin`, `ereq`...),
checked against the engine's keyword table: `epsr.keywords()` lists them and
`epsr.Settings.describe("feedback")` explains one. A dataset keyword given globally
(`rshmin`) applies to every sample that does not set it. File names, the number of
datasets, the run-state flags and `niter` are set by the library.

**The schedule.** `Project.create` does not run EPSR. The first `equilibrate` makes
the *scratch start*: a one-iteration run without a `.pcof` file, which creates it,
followed by `iterations - 1` iterations on the reference potential (`ereq = 0`).
`refine(ereq=..., runs=..., iterations=...)` then switches the empirical potential
on and refines; `run(iterations=...)` is a single run. `run` and `refine` refuse a
project that has never run. `set(**settings)` changes keywords in the files between
runs, `settings` shows the current values and `status()` summarises the state.
`record=True` records every iteration with `epsr.vis` (`project.trajectory()`), and
`epsr.analysis.first_peak` / `coordination_number` cover the common g(r) analysis.

**Seeds.** EPSR's random-number state is stored in the box file, so every run is
reproducible. `seed=None` keeps EPSR's default (identical to `seed=1`); any other
seed (1 to 2³¹ − 2) gives a different starting box.

The folder also holds `epsr_project.json` (the definitions, for `project.system` and
`project.samples`) and the generated `.mol` files; EPSR ignores both. The setup
programs work in the folder as their working directory, so build several folders
in parallel with processes, not threads.

## 3. A first run

```python
import epsr

proj = epsr.Project.from_template("EPSR-test", "work/dmso")
sim = proj.simulation()
summary = sim.run()
print(summary.r_factor, summary.energy, summary.pressure)
```

`Project.from_template` copies the regular files of a folder into a new one
(it refuses to overwrite a non-empty folder unless `overwrite=True`). **Always run
on a copy**: EPSR rewrites `.ato`, `.pcof`, `.EPSR.inp` and `.EPSR.inpa` at the
end of each run and appends a line to `.erg`, and creating a `Simulation` already
rewrites `.qdr`.

`proj.simulation()` reads every input (the Rust `Epsr::setup`). `run()` performs
`niter` iterations and writes the results. It returns a `Summary`:

| field | meaning |
|---|---|
| `r_factor`, `r_factor_per_dataset` | overall and per-dataset R-factor (fit quality) |
| `energy`, `pressure` | intermolecular energy per molecule (kJ/mol) and pressure (kbar) |
| `energy_inter`, `energy_intra` | running intermolecular and intramolecular energies per molecule |
| `ep_amplitude` | absolute energy of the empirical potential, at most `ereq` |
| `ep_coefficient`, `gr_roughness`, `chisq`, ... | the other quantities EPSR writes to `.out` |

`summary.fortran` holds the same values under the Fortran variable names.

### Iterating step by step

```python
sim = proj.simulation(overrides={"niter": 3})
for s in sim.iterations():          # yields a Summary after each iteration
    print(s.iteration, s.r_factor)
```

`iterations()` writes exactly the same files as `run()`. Iterations must run in
order and only up to `niter`; a finished `Simulation` cannot be run again. Create
a new one for the next run.

## 4. Inspecting results

All arrays are float32 numpy copies.

```python
q, r = sim.q, sim.r                 # Q grid (1/Å), r grid (Å)
sim.pair_labels                     # [('S', 'S'), ('S', 'C'), ..., ('H1', 'H1')]
sim.gr_inter()                      # (n_pairs, nr) intermolecular partial g(r)
sim.gr_inter("O1", "H1")            # one pair; label order does not matter
sim.gr_intra(3), sim.sq_inter(), sim.sq_intra()
sim.current_gr()                    # (2, n_pairs, nr) histograms of the current configuration only
```

The partials are running averages over the last `inter` iterations, exactly
what EPSR writes to `.g01`, `.y01`, `.f01` and `.s01`.

Data and fits, per dataset (0-based index, in the order of the `data n` blocks):

```python
for i, ds in enumerate(sim.datasets):
    f = sim.fit(i)                  # restricted to the fitted Q range
    print(ds.name, f.q[0], f.q[-1], abs(f.residual).max())
sim.data_fq(0), sim.simulated_fq(0), sim.data_errors(0), sim.data_q_range(0)
```

Configuration:

```python
sim.positions                       # (n_atoms, 3) Å
sim.atom_labels, sim.molecule_atoms # per-atom labels; atom index range of each molecule
sim.box_size, sim.density           # Å, atoms/Å³
sim.get_parameter("ntimes")         # current value of any input keyword, as text
sim.get_parameter("tweak", 1)       # a per-dataset keyword of dataset 1 (0-based)
```

Result files can also be read directly, without a `Simulation`:

```python
g = proj.read_table("g01")          # epsr.io.Table
g.x, g["O1-H1"], g.errors["O1-H1"]
g.to_dataframe()                    # with the pandas extra
proj.history()                      # the whole .erg file as a History
```

## 5. Changing parameters: overrides

Overrides change input keywords after the files are read, as if the lines had
been appended to the `.EPSR.inp` file. Keywords of the `.pcof` header (`ereq`,
`rminpt`, `Gaussian`, ...) work the same way.

```python
sim = proj.simulation(overrides={
    "ntimes": 10,
    "ereq": 15.0,                       # EP amplitude limit
    "potfac": [1.0, 1.0],               # lists become space-separated values
    "data": {2: {"tweak": 1.05}},       # dataset blocks use EPSR's 1-based numbers
})
sim = proj.simulation(overrides=["ntimes 10", "data 2", "tweak 1.05"])   # raw lines also work
```

Note that `data` numbers in overrides are 1-based like the input file, while
every other dataset index in the API is 0-based. A dataset can also be given by
name (`{"data": {"HDO": {"tweak": 1.05}}}`): the sample names of a folder made by
`Project.create`, or else the data-file names without extension.

Overrides are written back to the input files at the end of the run, so they
persist. An unknown keyword raises `EpsrInputError`.

`settings=` is the checked form of the same thing: EPSR keyword names validated
against the engine's table before anything runs (with "did you mean" for typos),
datasets by name or 0-based index. `epsr.Settings` holds such values, and
`project.set(**settings)` writes them into the files without running.

```python
sim = proj.simulation(settings={"ntimes": 10, "data": {"D2O": {"tweak": 1.05}}})
```

## 6. Refining: repeated runs

A real refinement is hundreds of EPSR runs. `epsr.refine` runs them in a row,
each with a fresh `Simulation`:

```python
hist = epsr.refine(
    proj, n_runs=200,
    overrides={"ntimes": 5},            # applied to every run
    initial_overrides={"ireset": 1},    # first run only
    callback=lambda k, sim: print(k, sim.summary.r_factor),  # return True to stop
    until=lambda h: h.r_factor[-1] < 5e-5,
    snapshot_every=50,                  # copies to <folder>.snapshots/run0050, ...
)
hist.r_factor, hist.energy, hist.ep_amplitude      # numpy arrays, one value per run
hist.r_factor_per_dataset                          # (n_runs, n_datasets)
hist.to_csv("history.csv"); hist.to_dataframe()
epsr.plot.history(hist)
```

A progress bar appears when tqdm is installed and stdout is a terminal
(`progress=True/False` forces it).

## 7. Perturbing the potential

`sim.potential` gives access to the potentials of a live simulation.

```python
pot = sim.potential
pot.empirical("O1", "H1")           # EP of one pair (kJ/mol on sim.r)
pot.reference(), pot.total()        # LJ + Coulomb; reference + empirical
pot.params                          # ereq, potfac, refpotfac, feedback, ncoeffp, pdstep, ...
pot.type_params                     # per atom type: epsilon, sigma, charge, mass, fraction
```

The EP of each pair is a sum of basis functions (Poisson or Gaussian, spacing
`pdstep`) with `ncoeffp - 1` coefficients; the engine builds the r-space tables
from them. Perturbations therefore act on the coefficients:

```python
c = pot.ep_coefficients()           # (n_pairs, ncoeffp - 1)
c[pot.pair_index("O1", "H1")] *= 0.5
pot.set_ep_coefficients(c)          # rebuilds the EP immediately
pot.perturb_ep(("O1", "H1"), 0.002)  # add to one pair's coefficients
pot.scale_ep(0.9)                   # scale every pair
pot.zero_ep()                       # like ireset 2
sim.run()                           # the edited EP drives the run and is saved to .pcof
```

**The ereq limit.** Whenever the EP is rebuilt, the engine measures its absolute
energy and scales all coefficients down so it does not exceed `ereq`. Scaling
down, zeroing and perturbing within the limit act as given. Scaling the EP beyond
`ereq` is undone. To allow a larger EP, raise the limit when you create the
simulation: `proj.simulation(overrides={"ereq": 20.0})`.

To run Monte Carlo with a fixed EP, switch refinement off with
`overrides={"potfac": [0.0, 0.0]}`; the coefficients are then saved unchanged.

The reference-potential parameters (`type_params`) are read-only; edit the `.ato`
file, for example with the `changeato` tool, to change them.

## 8. The setup and analysis tools

All 35 programs of the EPSR suite are available through `epsr.tools.run(folder,
command, args, answers)`, with the same arguments as on the EPSR command line.

```python
from epsr import tools
tools.commands()                                       # the 35 names

# keyword-file programs take their .dat file as argument
tools.nwts(proj.folder, "DMSODWaterD.NWTS.dat")        # neutron weights
tools.write_keyword_file(proj.folder / "w.PARTIALS.dat", "partials",
                         {"fnameato": "200DMSO400Water.ato", "atom-c": "O1"})
tools.partials(proj.folder, "w.PARTIALS.dat")
epsr.read_table(proj.folder / "w.PARTIALS.g01")

# interactive programs read their answers from standard input
tools.run(proj.folder, "writexyz", ["DMSO.ato"], answers=["Y", "1", "10 10", "-10 10", "10 20 30", "3 4"])
```

* `tools.keywords("nwts")` shows a program's keyword table, defaults and
  descriptions. `write_keyword_file` takes global keywords as a dict and repeated
  blocks (`component`, `distribution`, ...) as a list of dicts; missing keywords
  take the defaults.
* The answers each interactive program expects are easiest to copy from the
  `*.stdin` files of the cases in `test_problems/200DMSO400Water/EPSR-test-ex*`
  (their `README.md` explains each answer).
* Like the Fortran, `tools.run` silently does nothing when a keyword file is
  missing; the typed helpers (`tools.nwts`, `tools.partials`, ...) raise
  `FileNotFoundError` instead.
* Several programs use the current working directory (they read
  `system_commands.txt` and write `dirlist.txt` there). `tools.run` therefore
  switches to the folder for the duration of the call; tool runs are serialised
  by a lock.
* The setup programs also have typed wrappers, which `Project.create` uses
  (section 2) and which check that their output was written:

  ```python
  tools.makemole(folder, "Water")                              # Water.mol -> Water.ato
  tools.mixato(folder, [("Water", 1000)], density=0.1, box="waterbox")
  tools.randomise(folder, "waterbox", seed=42)                 # seed=None: EPSR's default
  tools.fmole(folder, "waterbox", 100)                         # moves *.mol aside for the run
  tools.nwts(folder, "D2O.NWTS.dat", box="waterbox", components=[
      {"atom": "OW", "atsymbol": "O", "iexchange": 0, "abundances": {0: 1.0}},
      {"atom": "HW", "atsymbol": "H", "iexchange": 1, "abundances": {2: 1.0}},
  ])                                                           # -> D2O.NWTStot.wts
  ```

## 9. Plotting

`epsr.plot` replaces the PGPLOT/gnuplot front ends, which are not ported. It needs
the `plot` extra.

```python
epsr.plot.fits(sim)                                 # data and fit for every dataset
epsr.plot.partials(sim, [("O1", "H1"), ("O1", "O1")])
epsr.plot.potentials(sim, which="empirical")
epsr.plot.history(hist, ["r_factor", "energy"])
```

Each function accepts `ax=` to draw into existing axes and returns the axes.

### Following a refinement iteration by iteration: `epsr.vis`

EPSR writes its result files, and computes the R-factor, only at the end of each
run, so `epsr.plot` shows the final state. `epsr.vis` records **every iteration**
and analyses the whole refinement. Replace `Simulation.run()` and `epsr.refine`
by their recorded versions. They run the same iterations
(`Simulation.iterations()`), so the folder ends up exactly as without recording
(`tests/test_vis.py` checks this file by file).

```python
from epsr import vis

recorder = vis.TrajectoryRecorder(case="water")      # case: settings in vis.systems
vis.run_recorded(project.simulation(), recorder)     # instead of simulation.run()
history = vis.refine_recorded(project, 2, recorder,  # instead of epsr.refine(...)
                              overrides={"niter": 100},
                              initial_overrides={"ereq": 3.0})
trajectory = recorder.trajectory()                   # arrays [iteration, ...]
recorder.save("work/trajectory")                     # vis.Trajectory.load(...) later

vis.make_report(trajectory, "figures/")              # 9 figures + summary.csv/.md
fit = vis.metrics.fit_metrics(trajectory)            # R-factor, R², chi² per iteration
vis.plots.partials_vs_reference(trajectory, vis.load_reference_partials({...}))
```

The figures show convergence (energy, pressure, R-factor, R², EP amplitude
against `ereq`, density, MC acceptance), the F(Q) fits and residuals, the
partials, the empirical potential, coordination numbers and bond angles. With
the `viewer` extra, `epsr.vis.viewer.ConfigurationViewer(trajectory).show()`
opens a 3D view of the configuration next to the partial g(r), with a slider
over the iterations. `screenshot()` and `movie()` (GIF or MP4) work off-screen.
`examples/05_experimental_water.py` and `examples/06_experimental_SiO2.py`
end with these figures.

`epsr.vis` needs the `plot` extra (matplotlib, scipy); it can always be imported,
and `epsr.vis.available()` says whether the extra is installed (without it, any
other name of `epsr.vis` raises `ImportError`). With `Project`, pass `record=True`
to `equilibrate`, `refine` or `run` instead of using the recorded functions, and
get the trajectory with `project.trajectory()`. It is the visualisation
toolkit of `visualisation/python_visualisation_toolkit` (standalone package
`epsrvis`, with documented examples). The shared modules are identical in both
places.

## 10. Parallel runs

**Inside one simulation.** The engine parallelises its expensive loops (energies,
pressure, structure factors, the S(Q) to g(r) transform, the potential fit) over
all CPU cores with rayon. Results are byte-identical at any number of threads. Set
`RAYON_NUM_THREADS` before starting Python to limit the cores used, for example
`RAYON_NUM_THREADS=4 uv run python script.py`; the pool size is fixed once the
first simulation runs.

**Several folders at once.** The engine releases the GIL while it runs, so
independent folders can also be refined from threads:

```python
projects = [epsr.Project.from_template("EPSR-test", f"work/run{i}") for i in range(4)]
histories = epsr.refine_many(projects, n_runs=50, overrides={"ntimes": 5})
```

Never use one folder from two threads at once. Tool runs change the working
directory and are therefore serialised. Concurrent simulations share the engine's
one pool of worker threads, so running several at once gives the same results
but does not make each of them faster; it helps most when `RAYON_NUM_THREADS` is
set low or the machine has many cores.

## 11. Errors

| exception | when | also a |
|---|---|---|
| `epsr.EpsrIOError` | a file is missing or unreadable | `OSError` |
| `epsr.EpsrInputError` | a malformed input file or an unknown override keyword | `ValueError` |
| `epsr.EpsrUnsupportedError` | an input on which the original Fortran itself crashes or hangs | `NotImplementedError` |

All three derive from `epsr.EpsrError`. The list of unsupported inputs is in
`rust/docs/FORTRAN_MAP.md` ("Still unsupported").

## 12. Troubleshooting

* *The progress messages do not appear in my notebook.* With `verbose=True` the
  engine prints from Rust straight to the process's standard output, which
  Jupyter does not capture. Watch the terminal that runs the kernel.
* *My edits to the EP had no effect.* Check whether the EP was already at `ereq`
  (see section 7), and whether refinement (`potfac`) changed it again during the run.
* *Results differ slightly between machines.* EPSR computes in single precision
  and uses the platform's libm; the last digits can differ outside Linux/glibc.
