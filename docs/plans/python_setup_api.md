# Plan: a Python setup API for EPSR (investigation and design)

**Status:** design agreed (2026-10-19); implementation in progress on `feature/python-setup-api`.
See the *status* column of the checklist and *Implementation notes* below it.
**For:** the agent that implements it. Read this whole document first. The *capability
checklist* (after the decisions) is the definition of done and the scope boundary. §5 is
the work plan, and §8 specifies the next example, Cu-BTC.

**Goal (from the user):** a user should be able to define a system, its data and its
settings, and run a refinement, with clear Python syntax only. They should never write
or copy an EPSR input file. Example 07 (`python/examples/07_water_from_scratch`) is the
reference: today it works around the missing API by writing `.mol` and `.EPSR.inp` text
itself, which the user does not want.

**Hard constraints**

1. **Byte-identity with the Fortran must survive.** The engine is file-based: `Epsr::setup`
   reads files, and each run rewrites them. The new API writes the same files EPSR's own
   tools would write and then drives the existing engine. It must not re-implement EPSR in
   Python. Acceptance test: the new example 07 produces an `epsr_run/` folder that is
   byte-identical to the committed `python/examples/07_water_from_scratch/output/epsr_run/`
   (equivalently, to `validation/water`). The same two files are excepted as today: the
   date line of `.EPSR.out`, and `.qdr`.
2. **Don't break the existing API or its tests** (`python/tests`, 82 fast and more slow).
   The new layer goes on top of the current one. Existing calls stay as low-level escape
   hatches.
3. Files are an implementation detail the user can still inspect. They are never
   something the user has to author.

**Decisions confirmed by the user (2026-10-19).** These override anything else in this
document that disagrees.

| # | decision | where it is specified |
|---|---|---|
| D1 | **Keep EPSR's own keyword names** (`ntimes`, `ereq`, `qmin`, …) as the setting names, validated against the engine's keyword table, with the descriptions as documentation. No second vocabulary. | §4.1 principle 2, §4.3 *Settings* |
| D2 | **The scratch run happens in the first equilibration step**, not in `Project.create`. `create` builds the folder only. The first `project.equilibrate(iterations=N)` on a new project runs the 1-iteration scratch start and then N − 1 iterations, reproducing the validated 1 + 99 schedule. | §4.3 *Project*, §6 item 6 |
| D3 | **Add a `seed=` option.** Where the seeds come from was verified (§1.7). `seed=None` must keep today's byte-identical output. | §1.7, §4.3, §5 Phase 1 |
| D4 | **Atoms in a `Molecule` have unique names** (`"O"`, `"H1"`, `"H2"`). Bonds, angles and rotations accept those names *or* 0-based indices. | §4.3 *Molecule* |
| D5 | **`exchangeable` is declared on the `AtomType`** (`AtomType(..., exchangeable=True)`). It applies to every sample. A `NeutronSample(exchangeable=...)` may override it for that sample only. | §4.3 *AtomType*, *samples* |
| D6 | **`System` takes number density only:** `density=` in atoms/Å³, EPSR's native unit. No `mass_density=`. | §4.3 *System* |
| D7 | **Isotopes as mass numbers only:** `isotopes={"HW": {1: 0.5, 2: 0.5}}`, i.e. `{atom type: {mass_number: fraction}}`, exactly EPSR's `abundances` pairs. Mass number 0 = natural element; omitted types are natural. No symbol aliases (`"D"`). | §4.3 *samples* |

**Validation rules from the engine** (enforced when the objects are built, with clear
messages):
- **Atom-type labels:** at most **3 characters**, because `makemole` and the `.ato`
  truncate them (`fstr(..., 3)` in `tools/makemole.rs`), so longer labels would silently
  collide. They must be unique across the `System`.
- **Atom names within a molecule** (D4): Python-side only, and unique per molecule. In
  the generated `.mol`, atoms are numbered `001…`.
- **Engine limits** (`state.rs`): at most 40 atom types (`MCOMP`), 30 datasets (`MFILE`)
  and 50,000 molecules (`MMOL`). At most 5 isotopes per atom type per sample
  (`MISOTOPE`, `tools/nwts.rs`).
- **No blanks** in system, molecule or sample names. They become file names.

---

## Scope: the capability checklist

**This table is the definition of done.** Implement the capabilities in this list. Each
one names its phase (§5), where it shows up in the wireframe
([example_07_syntax_example.md](example_07_syntax_example.md), "WF") or in the Cu-BTC
example (§8), and the test that proves it. Anything in *Out of scope* below must not be
started without the user's agreement, even if it looks small.

### In scope: liquids and amorphous systems with neutron data (example 07)

| id | capability | phase | shown in | acceptance | status |
|---|---|---|---|---|---|
| C1 | `AtomType`: ε, σ, charge, mass (default by element), `exchangeable` (D5), `coarse_grain`; label ≤ 3 chars; docstring states the §1.8 mixing rules | 2 | WF §1 | unit tests for validation; water `.ato` type table identical |  |
| C2 | `Molecule`: named atoms + bonds/angles/rotations by name or index (D4); `vibtemp`/`angtemp`/`dihtemp`; `render_mol`; `from_geometry`, `from_xyz`, `from_mol`, `from_ato`, `Molecule.atom`; rejects dihedrals | 2 | WF §1 | `makemole` on the rendered water and silica `.mol` gives `.ato` files identical to `validation/` |  |
| C3 | `System`: `add` order = type order, `density` in atoms/Å³ (D6), `temperature`, engine limits, `pair_parameters()` with a `terms` column | 2 | WF §2 | the water table of §4.3; limit errors |  |
| C4 | `NeutronSample`: isotopes as `{mass_number: fraction}` (D7), checked against the isotope table; per-sample `exchangeable` override; per-dataset keywords; `.mint01` validation | 3 | WF §3 | rendered `.NWTS.dat` → `nwts` → `.NWTStot.wts` identical to `validation/water`; bad isotopes raise |  |
| C5 | `Settings`: the keyword table binding; EPSR names (D1); unknown keys rejected; routing to `.inp`/`.pcof`/data blocks; `.inp` rendering | 0, 4 | WF §4 | table round trip; `.inp` after the scratch run identical to the validated one | Phase 0 part done: `_core.epsr_keywords()` with file routing |
| C6 | Typed setup tools: `makemole`, `mixato`, `randomise`, `fmole` (`.mol` handled), `nwts` from structured input | 1 | inside `create` | the water box rebuilt with these calls is identical to `validation/water` |  |
| C7 | `seed=` (D3, §1.7): `set_rng_seed`; `None` = default | 0, 1 | WF §4 | the §1.7 table: default = seed 1, 42 reproducible and different, 0 rejected | Phase 0 part done: `_core.set_rng_seed`; §1.7 table reproduced |
| C8 | `Project.create`: steps 1–6 of §1.6, no EPSR run (D2), `epsr_project.json`, pending `.pcof` settings, rejects `ereq ≠ 0` | 5 | WF §4 | folder contents before the first run identical to the validated case at that step |  |
| C9 | `equilibrate`: scratch start + (N − 1) on the first call (D2); `History` with both rows | 5 | WF §5 | `create` + `equilibrate(100)` identical to steps 7–8 of `validation/water` |  |
| C10 | `refine(ereq=, runs=, iterations=, until=, callback=, snapshot_every=)` | 5 | WF §5 | full water schedule identical to `validation/water` |  |
| C11 | `run`, `set`, `settings`, `status`; `run`/`refine` on a project that has never run raise an error | 5 | WF §4 | unit tests |  |
| C12 | `record=True` and `project.trajectory()` (all recorded runs concatenated) | 5 | WF §5, §7 | the trajectory has scratch + 99 + 200 steps; recording does not change any file |  |
| C13 | `project.results()`: no file modified (`.qdr` included) | 0 (optional Rust), 5 | WF §6 | folder hashes are unchanged before and after the call | Phase 0: no Rust change (decision in *Implementation notes*) |
| C14 | Datasets by sample name everywhere (`fit("D2O")`, settings `data` blocks, `Summary.r_factor_for`) | 5 | WF §6 | unit tests; integers still work |  |
| C15 | `pair_params()` from the engine: ε_ij, σ_ij, q_i q_j, `coulfac`, active `.pcof` extras, `terms` | 0, 5 | WF §6 | water values; EPSR-test-ex1 shows the extras | Phase 0 part done: `_core.Simulation.pair_params()` |
| C16 | `epsr.analysis.first_peak`, `coordination_number` | 5 | WF §6 | O–O 2.79 Å, O···H ≈ 1.8 Å on the water result |  |
| C17 | `epsr.vis.available()` | 5 | WF §5 | importable with and without the `plot` extra |  |
| C18 | Example 07 rewritten to WF; docs (`user_guide`, `api`, READMEs) | 6 | WF | **new slow test:** `epsr_run/` byte-identical to the committed output |  |
| C19 | Backwards compatibility: every existing public call keeps working (`overrides` as an alias of `settings`) | all | — | the existing `python/tests` pass unchanged |  |

### In scope next: crystalline systems from a PDB file (example 08, Cu-BTC, §8)

| id | capability | phase | acceptance | status |
|---|---|---|---|---|
| C20 | `Framework.from_pdb`: a structure from a PDB file through `readpdb` (bond detection by vdW factor, optional `DIH` groups), with the `vanderWaalsRadii.txt` shipped by the package and a check of the label order | 7a | `.ato` after `readpdb` identical to the validated run |  |
| C21 | Force field and Monte Carlo parameters applied to an existing configuration, through `changeato` (per-type ε/σ/charge/mass; `temp`, `tol`, step sizes), using EPSR's `changeato` keyword names (D1) | 7a | `.ato` after `changeato` identical to the validated run |  |
| C22 | Bragg settings as `Settings` with tuple values (`qwidthqmax`, `hklqmin`, `mplicities`); results readers for `.bragg` and `.inten` (`epsr.io`) | 7a | 80 726 reflections at 940 unique Q, as in `validation/cu-btc/README.md` |  |
| C23 | Example 08 (Cu-BTC) in the new syntax | 7a | every file byte-identical to `validation/cu-btc/fortran_reference.tar.gz` after each EPSR run |  |

### Implementation notes (settled while implementing)

- **Phase 0 (2026-10-19).**
  - **Keyword routing**, read from the engine's writers (`save_inp_file`, `save_pcof_file`
    in `symbols.rs`), with 1-based table positions: `.inp` header = 1–2, 10–46 and 65
    (40 keywords); `.pcof` header = 3–9 and 47–64 (25 keywords, `ereq`…`sfreq`,
    `roverlap`…`refpotfac`); `data n` block = 66–73 (8 keywords). `_core.epsr_keywords()`
    returns this as the `file` field (`"inp"`, `"pcof"`, `"data"`). Tested against the
    validated water `.inp` and a `.pcof` the engine wrote.
  - **`set_rng_seed`** reads the box with `minitc` and writes it with `moloutc` (as
    `randomise` does). On the mixato water box this rewrite is byte-identical except for
    the RNG record, and the §1.7 table is reproduced exactly.
  - **Isotope table units:** the Sears values are in 10⁻¹² cm (H natural −0.3739, i.e.
    −3.739 fm), not fm as §4.3 says.
  - **C13, read-only results:** the `.qdr` is written by `lsparabola1` inside
    `Epsr::setup` whenever `iinit != 1`. Skipping it would need a new engine flag, so the
    Rust side stays unchanged and `project.results()` will save and restore the `.qdr`
    in Python (Phase 5).
  - `pair_params()` is bound on `_core.Simulation` with Fortran names (`ecor`, `dcor`,
    `qcor`, `coulfac`, `rminex`, `pcoefex`, `coordex`, `rmaxex`, `pcoefexmax`,
    `gaussians`, `pcof_entry`); the readable Python layer comes in Phase 5.

### Out of scope until agreed with the user

X-ray samples (`xwts`); typed inputs for the analysis tools (`partials`, `coord`,
`triangles`, `torangles`, `sdf`, `sharm`, …); the other box builders (`makeato`,
`makelattice`, `makelatticeato`, `addato`, `dockato`, `readxyz` as a box builder);
`PairTerms` (setting per-pair `.pcof` extras); dihedral restraints from `Molecule` and
non-cubic boxes (both need the Rust `makemole` `dihedral`/`unitcell` support ported);
two-sample neutron mixtures (`nsamples`, `xfrac`, `abund2`); force-field importers
(§6 item 5a); an in-memory, file-free engine setup (§4.5 item 4); consolidating the
plotting packages (§6 item 8).

A wireframe of example 07 in the new syntax is in
[example_07_syntax_example.md](example_07_syntax_example.md).

---

## 1. Anatomy of the EPSR-test input (`EPSR-test/`)

A simulation folder is a set of files. One EPSR invocation (`epsr <folder> epsr <name>.EPSR.inp`)
reads all of them, runs `niter` iterations and rewrites most of them in place. That
rewrite is how state carries over to the next run.

| file | role | written by | read by |
|---|---|---|---|
| `200DMSO400Water.EPSR.inp` | run control and dataset list (§1.1) | user or EPSRshell; rewritten by every run | `epsr` |
| `200DMSO400Water.pcof` | potential settings and the EP coefficients per pair (§1.2) | created by the first (scratch) run; rewritten by every run | `epsr` |
| `200DMSO400Water.ato` | the box: molecules, coordinates, intramolecular restraints, atom types with their reference potential, RNG state (§1.3) | `mixato`/`randomise`/`fmole`…; rewritten by every run | `epsr`, every tool |
| `DMSO.mol`, `WaterMolecule.mol` | molecule definitions, "gmol" format (§1.4) | user | `makemole` → `DMSO.ato`, `WaterMolecule.ato` (one-molecule templates) |
| `DMSODWater{D,H,Null}.NWTS.dat` | isotopic composition of each sample (§1.5) | user | `nwts` → `*.NWTStot.wts` (weights matrix) |
| `D6DMSO_{D2O,H2O,HDO}.mint01` | measured data, Gudrun output (Q, F(Q), error) | instrument reduction | `epsr` |
| `.EPSR.{inpa,erg,out,g01,f01,u01,…}` | outputs: copy of the inputs, energy history, log, result tables | `epsr` | user, `epsr.io` |
| `system_commands.txt`, `plot_defaults.txt`, `*.gnu`, `runepsr.txt` | EPSRshell/gnuplot configuration | shell | not read by the port |

### 1.1 `.EPSR.inp`: sections and what they control

The format is one `keyword value [description]` per line, then `data n` blocks, then `q`.
Keywords are positional in a symbol table (`rust/epsr-core/src/symtab.rs`, `EPSR_SYMBOLS`,
generated from `fortran_src/setup_epsr_rev1.f`). **The same table also holds the `.pcof`
header keywords**, so `.inp` and `.pcof` share one namespace, and the overrides mechanism
already accepts both. Missing keywords take the table default. That is why example 07's
minimal `.inp` works.

The keywords fall into these groups. Values are from EPSR-test, with the defaults in brackets.

**A. Identity and file names.** The title line `200DMSO400Water.EPSR` is the run name and
the prefix of every output file. `fnameato` and `fnamepcof` name the box and potential
files.

**B. Empirical-potential feedback.** How hard the data push on the model:
- `feedback 0.9`: confidence factor, the weight of the data relative to the simulation
  when the difference is inverted.
- `potfac 1 1`: 0 freezes EP refinement.
- `rspcrmin 1.0`, `rspcfrac 0.5`: the r-space quality coefficient that drives automatic
  `ereq` control.

(The EP *amplitude* `ereq` and its automatic control live in the `.pcof`, §1.2.)

**C. Monte Carlo schedule and run state.**
- `ntimes 5`: MC cycles per iteration.
- `niter 1`: iterations (EP updates) per invocation.
- `intra 100`: iterations between molecule "shakes".
- `rotfreq 5`: internal-rotation moves.
- `inter 5`: length of the running averages.
- `rejrate 0.75`: target rejection rate, used to adapt the step sizes.
- `nmolcell 10`: molecules per neighbour cell (performance).
- `num_threds 0`: threads.
- **State flags, rewritten by EPSR itself:**
  - `ireset`: 1 is a full reset; 2 zeroes the EP.
  - `iinit 0`: 1 zeroes the accumulators and rebuilds the grids.
  - `nsumt 102`: iterations already accumulated.

  EPSR sets `ireset` and `iinit` back to 0 after a run. Pitfall 6.4: never regenerate the
  `.inp` after the first run.

**D. Box and thermodynamics.**
- `rho 0.0881`: number density, re-derived from the `.ato`.
- `sizefactor 1 0.9 0`: box scaling, a compression schedule for starting from an expanded box.
- `control-p none 0 0`: pressure control along axes a, b and c.
- `ref_intra 0 0`: EP and Coulomb weighting on intramolecular structure.

**E. Q and r grids, resolution and data treatment.**
- `nq 600`, `qstep 0.05`: Q grid. The data must fit within `nq × qstep`.
- `qmin 0.05 0`: the low-Q limit of the Fourier transforms and of the EP fit.
- `qinscalefc 1 0`: Q scaling of the input data, and broadening of the self terms.
- `cellst 0.03`, `rmaxgr 0`: r grid. 0 means half the box.
- `ngrsamples 0`: origin molecules sampled for g(r).
- `fwhm 0`, `fwhmq 0.02`, `revlorch 0 0`: instrument resolution and broadening.
- `nsmoop 1`: background subtraction on or off.
- `porodpwr 0 0`: Porod low-Q scattering.

**F. Crystalline systems** (unused in EPSR-test; used by Cu-BTC):
- Bragg peaks: `qwidthqmax`, `mplicities`, `hklqmin`.
- Diffuse scattering: `diffuse`.
- Lattice refinement: `nlatrefine`, `latrefine`.

**G. Datasets.** `ndata 3`, then one `data n` block per measurement:
- `datafile`, `wtsfile`: the data file and the weights file from `nwts` or `xwts`.
- `nrtype 5`: data format, where 5 means a Gudrun `.mint01`. Q values are treated as bin
  edges (pitfall 6.5).
- `rshmin 0.7`: r below which g(r) structure counts as background.
- `szeros 0`: the Q = 0 limit.
- `tweak 1.0`: scale factor on the data.
- `efilereq 1.0`: this dataset's share of the EP amplitude.
- `porodadd 0`: extra Porod amplitude.

### 1.2 `.pcof`: the potential

**Header keywords** (the same symbol table, positions 2–9 and 47–64):

- **EP amplitude and its control.**
  - `ereq 10`: the EP's absolute-energy cap. 0 means no EP. A scratch start has `ereq 0`.
  - `ereqmin`, `ereqmax`, `ereqstep`, `thresh`, `bias`, `sfreq`: automatic `ereq` control.
- **Reference-potential shape.**
  - `roverlap`, `rminfac`: minimum separations.
  - `rminpt 9`, `rmaxpt 12`: truncation range.
  - `rbroad`: potential decay.
  - `expecf 0.3 0`, `reppottype exponential`, `power 12`: the repulsive core.
  - `addpottype Gaussian`: the type of the additional per-pair potential.
  - `rcharge`: polarisation.
  - `refpotfac 1`: a scale factor on the whole reference potential.
- **EP basis.**
  - `ncoeffp 116`, `pdmax 12`, `pdstep 0.1`: the basis grid. EPSR sets these.
  - `Gaussian` versus `Poisson`: the type of basis function.
  - `psigma2 0.01`: the basis function width.
  - `npitss 1000`, `paccept 0.0005`: the fit of the coefficients.

**Per-pair blocks** (21 pairs for 6 atom types). Each block has three lines:
1. The labels, the pair index and `coulfac`, the Coulomb factor for the pair.
2. Extra per-pair terms: a coordination-constraint or extra-potential amplitude, then
   optional additional Gaussian or Morse terms.
3. The refinement flag `pcoefp(1)`, followed by the EP coefficients.

**Pitfall 6.3:** a `.pcof` without pair blocks silently switches the charges off and
freezes the EP. A scratch start must have *no* `.pcof`, so that the engine creates it.

### 1.3 `.ato`: the box

The file is laid out as follows:

1. **Header.** `nmol`, the box side and the temperature for a cubic box. A non-cubic box
   writes `nmol` and the temperature, then a, b, c, then the cell angles.
2. **MC parameters.** `tol`, then the four step sizes:
   - `stepmi`: intramolecular translation
   - `stepri`: intramolecular rotation
   - `stepr`: molecular rotation
   - `stepm`: molecular translation

   The same line continues with the restraint "temperatures" (`vibtemp`, `angtemp`,
   `dihtemp`) and `ecore`/`dcore`.
3. **One block per molecule.** Each block lists its atoms, their coordinates and their
   bonded-distance restraints to other atoms. It ends with rotational groups (`ROT`) and
   dihedral groups.
4. **The atom-type table.** For each type: its label, its element and whether it is
   exchangeable, then ε, σ, mass and charge, followed by the coarse-graining sphere
   parameters.
5. **Trailer.** The random-number generator state (§1.7), a summary line per molecule type, the number density and
   the atom fractions.

**The atom-type order is set by the order in which molecule types are added** (in
`mixato`). It defines the pair order of every pair-resolved array. Pairs mix by
Lorentz–Berthelot: ε_ij = √(ε_i ε_j), σ_ij = (σ_i + σ_j)/2, plus Coulomb 1390 q_i q_j / r
(`rust/epsr-core/src/ato.rs:785-800`).

### 1.4 `.mol`: a molecule definition (`makemole` input)

These lines are read by `makemole` (`rust/epsr-core/src/tools/makemole.rs`):

- `atom i label x y z nneigh neigh...`: atoms, with their coordinates and bonded neighbours.
- `bond i j length`: bond length restraints.
- `angle i j k degrees`: angle restraints.
- `rot i j`: rotatable bonds. DMSO uses them for its methyl groups.
- `changelabel i TYPE`: assigns an atom to an atom type.
- `potential TYPE ε σ mass q element nsph rmin rmax cg`: the reference potential of a type.
  The last four values are coarse-graining parameters.
- `temperature`, `vibtemp`, `angtemp`, `dihtemp`, `ecoredcore`, `density`, `pseudofac`:
  the temperature, the restraint stiffnesses and the template density.

**The Rust port does not support `dihedral` or `unitcell` lines** (they stop with
"unsupported").

### 1.5 `.NWTS.dat`: the isotopic composition of a neutron sample (`nwts` input)

**Global keywords:**
- `fnameato`: the box file.
- `output`: weights per atom or per molecule.
- `nsamples`, `xfrac`: mixtures of samples.
- `normtot`, `normdif`: normalisation.
- `ncomponent`: the number of atom types.

**One `component` block per atom type:**
- `atom`, `atsymbol`: the type label and its element.
- `iexchange`: whether the hydrogen is labile.
- `abundances`: pairs of mass number and fraction, where 0 means natural. Natural is
  `0 1`, D is `2 1`, and a 50:50 H/D mix is `1 0.5 2 0.5`.
- `isubs`, `abund2`: the second sample.

EPSR-test's three samples differ only in their H and H1 abundances. `nwts` writes the
`.NWTStot.wts` weights matrix that `wtsfile` points to. X-ray data use `xwts` and
`.XWTS.dat` in the same way.

### 1.6 How the EPSR-test box was built (the setup pipeline)

1. `makemole DMSO`, then `makemole WaterMolecule`: write the one-molecule templates.
2. `mixato`: 200 DMSO + 400 water at the template density → `200DMSO400Water.ato`.
3. `randomise`: random positions and orientations.
4. `fmole <box> N`: relax the intramolecular geometry.
5. `nwts` for each sample: write the weights.
6. Write the `.inp`.
7. Run `epsr` from scratch with `ireset 1`, `iinit 1` and no `.pcof`: this creates the
   `.pcof`.
8. Equilibrate with `ereq 0`.
9. Set `ereq > 0`.
10. Refine repeatedly.

`validation/PYTHON_WRAPPER_WORKFLOW.md` documents this call by call. Its §6 lists the
pitfalls every builder must absorb.

### 1.7 Random numbers and seeds (verified 2026-10-19)

**The generator.** `rust/epsr-core/src/rng.rs` is the Numerical Recipes `ran1`
(Park–Miller with a Bays–Durham shuffle), translated from `auxiliary_routines.f`. Its
complete state is 35 integers: `mcall`, `idum`, `iy` and the 32-entry table `iv`. **No
program takes a seed argument.** The state is read from the `.ato` trailer when a file is
loaded (`ato.rs:628`, `initran1`) and written back on every save (`ato.rs:1170`,
`saveran1`). That is why a refinement can be continued and reproduced exactly.

**Where the first seed comes from.**
- `makemole` writes a template `.ato` with an all-zero state (`0 0 0 0 …`).
- `mixato` makes no random draws, and its box `.ato` keeps the zero state.
- The first draw happens in `randomise`. For a zero state, `ran1()` takes the
  `idum <= 0` branch and sets `idum = max(-idum, 1) = 1`. **The default seed is therefore
  effectively 1, for every system built this way.**
- Every later program continues from the saved state: `fmole`, every `epsr` run, and the
  analysis tools that use `randomno.dat` (`anal_common.rs`).

The only other setup tools that draw random numbers are `addato`, `changeato` and
`dockato`. `makelattice` reads a state from its input (`tools/makelattice.rs:264`).

**Experiment** (scratchpad only, no repository change): build the water box with
`makemole` and `mixato`, write `idum` into the box's RNG record, then run `randomise`.

| `idum` written before `randomise` | state after `randomise` (mcall, idum, iy) | box |
|---|---|---|
| unchanged (0), twice | 6000, 1715014484, 745054426 | identical both times |
| −1 | 6000, 1715014484, 745054426 | **identical to the default**, so the default is seed 1 |
| −42, twice | 6000, 1163647977, 1007254037 | different from the default, identical both times |

(`randomise` makes 6 draws per molecule, so 1000 molecules give `mcall = 6000`.)

**How `seed=` should work.**
- Before the first random draw (after `mixato`, before `randomise`), set `idum = -seed` in
  the box's RNG record, leaving the rest of the record at zero. A negative `idum` is the
  standard way to initialise `ran1`.
- `seed=None` leaves the record untouched, which is equivalent to `seed=1` and keeps
  byte-identity with `validation/water` and the committed example 07 output.
- Valid seeds are 1 ≤ seed ≤ 2³¹ − 2. Reject 0 and negative values: 0 would silently
  behave like 1.
- **Implementation:** add a small Rust helper, `set_rng_seed(ato_path, seed)` in
  `python/src/tools.rs`, that loads the `.ato` with the engine's own reader and writes it
  back with `saveran1`, so the file layout stays exactly the engine's. A Python text edit
  of the 35-integer record also works (the experiment did that, and `randomise` rewrites
  the whole file anyway), but the Rust route avoids format drift. Do **not** add a seed
  argument to `randomise` itself: that would change the tool's command line and its
  parity with the Fortran.
- **Systems that skip `randomise`** (crystals from `readpdb` or `makelattice`): apply the
  seed to the box `.ato` before the scratch run instead, because the first `epsr` run then
  makes the first draw. Store the seed in `epsr_project.json`.

### 1.8 What the model's interactions actually are (verified 2026-10-19)

The new API must make these rules visible. Several are not obvious from the inputs.

**Intermolecular pairs: always derived from the two atom types.**
(`ato.rs:785-800`, `generate_rp` in `potential.rs:571`)

- ε_ij = √(ε_i ε_j), σ_ij = (σ_i + σ_j)/2 (Lorentz–Berthelot), and a Coulomb prefactor of
  1390·q_i q_j kJ mol⁻¹ Å.
- **No per-pair Lennard-Jones input exists.** The LJ term (n–6, with n = `power`, default
  12) is computed only when ε_ij > 0. **A type with ε = 0 therefore has no LJ site:**
  every pair involving it is Coulomb-only. In the water model, OW–HW (−467.3/r) and
  HW–HW (+233.7/r) are pure Coulomb. That is the SPC convention: H sits inside O's LJ
  sphere, and the hydrogen bond comes from the O–H Coulomb attraction.
- **Coulomb** is between point charges, or uniformly charged spheres if `qradius` is set
  (`chsph`). There is no Ewald sum. A polynomial taper takes it to zero at `rmaxpt`
  (`truncatec`). The LJ term has a cosine taper between `rminpt` and `rmaxpt`
  (`truncate`).
- **`refpotfac`** scales the whole reference potential.
- **The EP** is added per pair once `ereq > 0`. That's where the data reshape each pair,
  including O···H.
- **Optional per-pair extras in the `.pcof`**, the only per-pair knobs that exist:
  - `coulfac`, a Coulomb scale factor
  - an extra repulsive core below `rminex` (exponential or harmonic, `reppottype`, with
    `expecf`)
  - additional Gaussian or modified-Morse terms (`addpottype`)
  - minimum separations (`roverlap`, `rminfac`)

  Today they can only be set by editing the `.pcof` by hand.

**Intramolecular interactions: harmonic restraints, EPSR style.**
(`en_atom_intra` and `en_mole_dihedral` in `energy.rs:138-290`)

- **Bonds:** each bonded pair is restrained harmonically about its target distance r₀:
  E ∝ `vibtemp` · RT · √μ_ab / r₀ · (r − r₀)², where μ_ab = m_a m_b/(m_a + m_b) is the
  reduced mass of the two atom types and RT = 0.008314 · T kJ/mol.
- **Angles** are *not* harmonic in θ. Each angle is converted by `makemole` into a 1–3
  distance restraint (stored with a negative sign in the `.ato`) and restrained the same
  way, scaled additionally by `angtemp`.
- **Dihedrals:** for `DIH` rotational groups,
  E = `dihtemp` · `vibtemp` · RT · [(cos φ − cos φ₀)² + (sin φ − sin φ₀)²]. The engine
  implements this, but the Rust `makemole` cannot generate `dihedral` lines, so today
  dihedrals can only enter through a ready-made `.ato`.
- **Stiffness is not a per-bond force constant.** It is set once per molecule type
  (`vibtemp`, `angtemp`, `dihtemp`) and scaled by √μ / r₀ and the temperature. This is
  Soper's design: it gives bond-length fluctuations of roughly the zero-point
  vibrational width. **Consequence:** Dissolve-style per-bond constants such as
  `Harmonic k=4431.53` in `experimental_data/water/dissolve/input.txt` have no EPSR
  equivalent and cannot be transferred.
- **Units caveat:** the restraint energy is evaluated in the engine's reduced
  coordinates (`stepcl = cellst/dmaxmol`). Before the API reports an *effective* force
  constant in kJ mol⁻¹ Å⁻², verify the conversion against `en_atom_intra` and an `.erg`
  intramolecular energy.
- **Non-bonded atom pairs within a molecule** (further apart than 1–3) feel:
  - a soft WCA repulsive core with ε = `ecore`, σ = `dcore` (`ecoredcore`; default 0 0,
    i.e. off; `generate_mp`)
  - optionally the reference potential and the EP, weighted by `ref_intra` (default
    0 0, i.e. off; `epsr.rs:300-318`)

  So by default nothing acts between them.
- **Moves:** intramolecular moves are atom displacements every `intra` iterations (the
  "shake") and rotations of `rot` groups every `rotfreq` iterations (`molfit` in
  `montecarlo.rs:232`), accepted by Metropolis on the restraint energy.

---

## 2. What the Python wrapper offers today, per functionality

`python/epsr/`: `simulation.py`, `project.py`, `refine.py`, `potential.py`, `tools.py`,
`io.py`, `plot.py` and `vis/`, with the bindings in `python/src/*.rs`.

| functionality | current API | gap |
|---|---|---|
| **Open a folder and copy a template** | `Project(folder)`, `Project.from_template(src, work)`, `copy_folder`, `project.snapshot(label)` | none for existing folders; **no way to create a project from definitions** |
| **A run (setup + `niter` iterations)** | `Simulation(folder, inp, overrides=…)`, `.run()`, `.iterate(i)`, `.iterations()`, `.completed` | creating a `Simulation` has a side effect: it rewrites `.qdr` |
| **Change `.inp`/`.pcof` keywords** (groups A–F and the `.pcof` header) | `overrides={"niter": 99, "ereq": 3.0, "data": {2: {"tweak": 1.05}}}`, raw lines, `render_overrides` | untyped and unvalidated, so a typo surfaces as an engine error or silently not at all; dataset numbers are 1-based here but 0-based everywhere else; overrides **persist** (EPSR writes them back), which the name hides |
| **Read the current settings** | `sim.control` (`Control`, Fortran names), `sim.get_parameter(name, dataset)`, `sim.potential.params` (`PotParams`), `project.input_file()` (`KeywordFile`) | no single typed "settings" view; `Control` covers only part of the table |
| **Datasets (G)** | `sim.datasets` (`Dataset`), `sim.fit(i)`, `data_fq`, `simulated_fq`, `data_errors`, `data_q_range` | selection by integer only, not by sample name |
| **Repeated runs and the EP schedule** | `epsr.refine(project, n_runs, overrides, initial_overrides, callback, until, snapshot_every)`, `refine_many`, `History` | "switch on `ereq` once" needs `initial_overrides`; there is no equilibrate/refine vocabulary; recording needs separate functions (`vis.run_recorded`, `vis.refine_recorded`) |
| **Potentials and the EP** | `sim.potential`: `.empirical()`, `.reference()`, `.total()`, `.ep_coefficients()`, `.set_ep_coefficients()`, `.perturb_ep()`, `.scale_ep()`, `.zero_ep()`, `.ep_flags()`, `.shared_pairs()`, `.type_params` | no mixed pair parameters (example 07 recomputes Lorentz–Berthelot by hand); no access to `coulfac` or the extra per-pair terms |
| **Results** | `gr_inter/gr_intra/sq_inter/sq_intra(pair)`, `current_gr/sq`, `positions`, `atom_labels`, `molecule_atoms`, `box_size`, `density`, `summary` (`Summary`); `project.read_table("g01")`, `project.history()` | fine; peak positions and coordination numbers are done by hand (helpers exist in `epsr.vis.metrics`) |
| **Molecule definition (`.mol`, §1.4)** | **none.** Example 07 writes the text by hand (`mol_file_text()`), then calls `tools.run(folder, "makemole", ["Water"])`. | missing entirely |
| **Building the box (`mixato`, `randomise`, `fmole`; `makeato`, `makelattice`, `readpdb`, `readxyz`, `addato`, `dockato`)** | `tools.run(folder, cmd, args, answers=[…])` with positional stdin answers | no typed helpers; the user must know the prompt order; the `fmole`/`.mol` workaround (pitfall 6.2) is left to the user |
| **Isotopic composition and weights (`.NWTS.dat`, §1.5)** | `tools.write_keyword_file(path, "nwts", globals, blocks)` + `tools.nwts(folder, dat)`; `tools.keywords("nwts")` | Fortran keyword names and abundance strings (`"1 0.5 2 0.5"`); one atom type per block, by hand |
| **Writing the `.inp`** | **none.** Example 07 writes the text (`epsr_input_text()`). | missing; and **the EPSR keyword table is not exposed to Python**: `tools.keywords("epsr")` fails, only the 23 tool tables are exposed (`python/src/tools.rs`) |
| **Analysis tools** (`partials`, `coord`, `triangles`, `torangles`, `sharm`, `sdf`…) | `tools.run`; typed thin wrappers for 6 programs that take a pre-written `.dat` | the `.dat` still has to be written via `write_keyword_file` |
| **Plots and iteration recording** | `epsr.plot.fits/partials/potentials/history`; `epsr.vis` (`TrajectoryRecorder`, `run_recorded`, `refine_recorded`, `make_report`, `metrics`, `plots`, viewer) | recording is a parallel API, not an option of run/refine |
| **Errors** | `EpsrError`, `EpsrIOError`, `EpsrInputError`, `EpsrUnsupportedError` | fine |

**Summary:** everything from "the folder exists" onwards (running, reading, perturbing,
refining, plotting) is well covered. Everything *before* that (molecules, box, weights,
input file) is either raw file text or positional stdin.

---

## 3. Example 07 today: where the friction is

`python/examples/07_water_from_scratch/water_from_scratch.py` (377 lines) does the
following:

1. **Constants as dicts** (`ATOM_TYPES`, `ATOMS`, `BONDS`, `ANGLES`, `SETTINGS`, `DATASETS`).
   Good content, but no types or validation.
2. **`mol_file_text()`** renders the gmol format by hand. This includes the f32 rounding
   trick (`repr(float(np.float32(v)))`), the neighbour lists and the four magic trailing
   values on `potential` lines (`1 0 RADMAXSPH 1`, coarse-graining parameters copied to
   match the validated `.ato`).
3. **`tools.run(..., "mixato", answers=["1", "Water", "1000", "0.1", "waterbox"])`** uses
   positional prompts.
4. **`(RUN_FOLDER / "Water.mol").unlink()`** works around pitfall 6.2 by hand.
5. **NWTS via `write_keyword_file`** with Fortran keywords (`atsymbol`, `iexchange`,
   `abundances "1 0.5 2 0.5"`). `ncomponent` must equal the number of types.
6. **`epsr_input_text()`** hand-writes a minimal `.inp`, the part the user objects to. The
   user has to know that `fnamepcof` must name a file that does not exist yet (pitfall 6.3).
7. **The schedule** is a scratch run (`project.simulation()` + `run`), then
   `overrides={"niter": 99}`, then `refine(..., initial_overrides={"ereq": 3.0})`, with
   `if recorder is None` branching between `epsr.refine` and `vis.refine_recorded`.
8. **Results** recompute Lorentz–Berthelot pair parameters, and the peak positions use
   manual windows.
9. **Data paths** read `.mint01` from `validation/water`. They should come from
   `experimental_data/water/data/`, which holds identical files.

Items 2–7 are what the new API removes.

---

## 4. Proposed API

### 4.1 Design principles

1. **Declarative objects → generated files → existing engine.** New Python value objects
   describe the system. One call renders them to exactly the files EPSR's own tools expect
   and runs those tools. Running, reading and refining stay on the current `Project` and
   `Simulation`.
2. **Fortran keyword names stay canonical for settings.** EPSR users know them, the
   engine and docs use them, and the keyword table gives free validation and descriptions.
   Put readable names in docstrings, not in a second vocabulary. Physical-object
   parameters (atom types, molecules, samples) get readable Python names (`epsilon`,
   `sigma`, `charge`, `isotopes`, `exchangeable`), because there is no table behind them.
3. **Every pitfall in `PYTHON_WRAPPER_WORKFLOW.md` §6 is handled inside the library**, not
   by the user: chdir and lock, `.mol` removal before `fmole`, no `.pcof` on a scratch
   start, and never regenerating the `.inp` after the first run.
4. **Datasets are addressed by sample name** (`"D2O"`) everywhere. Integers still work.
5. **Low-level calls remain** (`tools.run`, overrides as lines, `write_keyword_file`) as
   escape hatches.

### 4.2 Target user code: example 07 rewritten

> **Superseded for the details:** the sketch below predates decisions D2 and D4–D6. The
> canonical version is [example_07_syntax_example.md](example_07_syntax_example.md), which
> has named atoms, `exchangeable` on `AtomType` and the scratch start inside `equilibrate`.

```python
import epsr
from epsr.system import AtomType, Molecule, System
from epsr.samples import NeutronSample

OW = AtomType("OW", element="O", epsilon=0.6503, sigma=3.165492, charge=-0.82)
HW = AtomType("HW", element="H", epsilon=0.0, sigma=0.0, charge=0.41, mass=2.0)

water = Molecule(
    "Water",
    atoms=[(OW, (0.0, 0.0, 0.0)), (HW, (0.976, 0.0, 0.0)), (HW, (-0.28754, 0.93268, 0.0))],
    bonds=[(0, 1), (0, 2)],          # 0-based; restraint lengths default to the geometry
    angles=[(1, 0, 2)],              # restraint angle defaults to the geometry
    vibtemp=65, angtemp=1, dihtemp=10,
)
# alternatives: Molecule.from_xyz("water.xyz", types={"O": OW, "H": HW}, bonds="auto"),
#               Molecule.from_mol("Water.mol"), Molecule.from_ato("EPSR-test/WaterMolecule.ato")

system = System("waterbox", temperature=300.0, density=0.1)   # atoms/Å³, cubic box
system.add(water, 1000)

data = Path("experimental_data/water/data")   # any folder of .mint01 files
samples = [
    NeutronSample("D2O", data / "SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}, exchangeable=["HW"]),
    NeutronSample("HDO", data / "SLS18500-HDO5050.mint01", isotopes={"HW": {1: 0.5, 2: 0.5}}, exchangeable=["HW"]),
    NeutronSample("H2O", data / "SLS18498-H2O.mint01", exchangeable=["HW"]),   # natural abundance
]

project = epsr.Project.create(
    "work/water", system, samples,
    settings={"qmin": 0.05, "rshmin": 0.7},   # EPSR keyword names, validated against the keyword table
    relax_iterations=100,                     # fmole
    seed=None,                                # default: effectively seed 1, byte-identical to validation/water
    overwrite=True,
)   # makemole → mixato → (seed) → randomise → fmole → nwts → .EPSR.inp; EPSR has not run yet

project.equilibrate(iterations=100)   # first call: 1-iteration scratch start + 99 on the reference potential (ereq = 0)
history = project.refine(ereq=3.0, runs=2, iterations=100, record=True)

result = project.results()                # a Simulation over the finished folder
result.fit("D2O").residual
result.gr_inter("OW", "OW")
result.potential.pair_params()            # mixed ε_ij, σ_ij, q_i q_j per pair
```

This replaces items 2–7 of §3 with about 25 lines. The figures part of example 07 stays
as it is, apart from using `history.trajectory` and `result`.

### 4.3 New objects and functions

**`epsr.system`** (new module; pure Python, rendering only):

- `AtomType(label, element, epsilon, sigma, charge, mass=None, *, exchangeable=False, coarse_grain=None)`.
  `exchangeable` (D5) becomes `iexchange 1` for this type in every sample's `.NWTS.dat`.
  The **docstring must state the §1.8 rules**: pairs are always Lorentz–Berthelot mixed;
  ε = 0 means the type has no LJ site, so all its pairs are Coulomb-only; and per-pair
  changes are only possible through `PairTerms` (Phase 7).
  `mass` defaults to the natural isotopic mass of `element`, which needs a small element
  table; check `rust/epsr-core/src/tools/epsrwts.rs` and `fortran_src/chemical_symbols.inc`
  for one to reuse. `coarse_grain` holds `nsph`, `rmin`, `rmax` and `cg`. Its default must
  reproduce what makemole writes when the values are omitted: (1, 0, 0.5·σ, 1). Example 07
  passes `rmax = 1.55` for OW explicitly to match the validated `.ato`, so expose `rmax`.
- `Molecule(name, atoms, bonds=(), angles=(), rotations=(), *, bond_lengths=None, angle_values=None, vibtemp=65, angtemp=1, dihtemp=10, ecore=0, dcore=0)`.
  `atoms` is a list of `(atom_name, atom_type, (x, y, z))`. Bonds, angles and rotations
  refer to atoms by `atom_name` or 0-based index (D4). Example:
  `bonds=[("O", "H1"), ("O", "H2")]`, `angles=[("H1", "O", "H2")]`.
  `molecule.atom("H1")` returns the atom, and the names are kept in `epsr_project.json`
  for geometry-table round trips. Methods:
  - `.from_xyz(path, types, bonds="auto"|list)`
  - `.from_mol(path)`
  - `.from_ato(path)`, for the EPSR-test templates
  - `.atom(atom_type)`, a one-atom "molecule" for silica-type systems
  - `.render_mol() -> str`, which also handles the f32 formatting

  Dihedrals: raise `EpsrUnsupportedError` at construction, because the Rust makemole
  does not support them.

  The **docstring must explain the restraint model** (§1.8):
  - bonds and angles take *target geometry only*, not force constants
  - stiffness comes from `vibtemp`, `angtemp` and `dihtemp` per molecule
  - angles are restrained as 1–3 distances

  Accepting a `k=` argument would be misleading, so don't. `Molecule.from_dissolve(...)`
  or a geometry table importer should drop force constants with a warning.
- `System(name, temperature, density, *, box=None)`. `density` is the atomic number
  density in atoms/Å³ only (D6), as passed to `mixato`.
  - `.add(molecule, count)`. **The order of the `add` calls fixes the atom-type order.**
    Document this, because it sets the pair order.
  - `.atom_types`
  - `.pair_parameters()`: Lorentz–Berthelot, as in `ato.rs`, for printing before a run.
    It returns one row per pair, in engine pair order: ε_ij, σ_ij, 1390·q_i q_j, and a
    **`terms` column** (`"LJ + Coulomb"`, `"Coulomb only"`, `"LJ only"` or `"none"`),
    so that a zero ε is visible. For water:

    | pair | ε_ij | σ_ij | 1390 q_i q_j | terms |
    |---|---|---|---|---|
    | OW–OW | 0.6503 | 3.165 | +934.6 | LJ + Coulomb |
    | OW–HW | 0 | 1.583 | −467.3 | Coulomb only |
    | HW–HW | 0 | 0 | +233.7 | Coulomb only |

    `result.potential.pair_params()` (§4.4 #7) returns the same table from the running
    engine, plus `coulfac` and any `.pcof` per-pair extras that are active. Both should
    also be available as a pandas `DataFrame` when pandas is installed.

  Later: `System.from_pdb(...)` and `System.lattice(...)` for crystalline cases (Cu-BTC).

**`epsr.samples`** (new module):

- `NeutronSample(name, data_file, isotopes=None, *, exchangeable=None, nrtype=5, rshmin=None, tweak=1.0, efilereq=1.0, szeros=0.0, porodadd=0.0)`.
  - `exchangeable=None` takes each `AtomType.exchangeable` (D5). A list overrides it for
    this sample.
  - `isotopes` maps a type label to `{mass_number: fraction}` (D7). Types not given are
    natural. Validate: fractions sum to 1 (within 1e-6), at most 5 isotopes, each mass
    number present for the type's element in the isotope table. It renders one-to-one,
    in the given order, to EPSR's format. EPSR's format, for reference: in the
    `.NWTS.dat`, each component has `abundances` as pairs of mass number and fraction,
    where mass number 0 means the natural element. So natural is `0 1`, pure D is
    `2 1`, and 50:50 H/D is `1 0.5 2 0.5`.
    - Fractions are normalised by `nwts`.
    - Scattering lengths come from the Sears (1991) table compiled into the engine
      (`tools/sears91.rs`; e.g. H: 0 → −0.3739, 1 → −0.37406, 2 → 0.6671, 3 → 0.4792 fm).
    - **An unknown mass number makes `nwts` prompt on stdin** ("Did not find massno …,
      try typing again"). The wrapper must validate isotopes against that table first.
      Expose `sears91::SEARS91` through a small binding, `_core.isotope_table()`.
    - There are also two-sample mixture keywords (`nsamples`, `xfrac`, `isubs`,
      `abund2`). They are out of scope for now; keep `nsamples 1`.
  - `rshmin=None` falls back to the global setting or the default, 0.7.
  - Validate the data file on construction. Port `archive/mellor/src/mellor/data.py`
    (`read_mint01`) and the checks from the deleted web app at commit
    `f1fe713:webapp/src/epsrweb/storage.py`: at most 5000 points, data beyond
    `nq × qstep`, no blanks in file names.
- `XraySample(...)`, the same with `xwts`. This can come later.
- Each sample renders its `.NWTS.dat` from the `System` (one component block per atom
  type, in system order).

**`epsr.tools` typed setup wrappers** (thin, testable on their own; `Project.create` uses
them):

- `tools.makemole(folder, molecule_or_name)`, which writes the `.mol` when given a `Molecule`
- `tools.mixato(folder, components=[("Water", 1000)], density=0.1, name="waterbox")`, which
  builds the answers
- `tools.randomise(folder, box)`
- `tools.fmole(folder, box, iterations)`, which moves any template `.mol` out of the way
  (pitfall 6.2), runs, and restores it
- `tools.nwts(folder, sample=..., system=...)`, overloaded or as a new name, keeping the
  existing `nwts(folder, dat_file)` signature
- Later, as needed: `makeato`, `makelattice`, `readpdb`, `changeato`, `readxyz`, `addato`.
  The `*.stdin` files of `test_problems/*/EPSR-test-ex*` give the answer order for each.

**Settings** (`epsr.settings` or inside `simulation.py`):

- `epsr.keywords()` returns the EPSR keyword table: name, type, default, description, and
  whether the keyword goes in the `.inp` header, a data block or the `.pcof`. This **needs
  a new binding** (§4.5).
- `Settings`, a validated mapping:
  - It rejects unknown keywords, with "did you mean".
  - It knows which file and which dataset block each keyword belongs to.
  - It renders a complete `.inp` with the shared layout from `tools.write_keyword_file`
    (`_keyword_line`).
  - `Settings.describe("feedback")` prints the description.

  Keep Fortran names.

**`Project` additions** (`project.py`):

- `Project.create(folder, system, samples, settings=None, *, relax_iterations=100, randomise=True, seed=None, overwrite=False, verbose=False)`:
  - It does steps 1–6 of the §1.6 pipeline: templates, box, seed, randomise, fmole,
    weights and `.EPSR.inp`. **It does not run EPSR (D2).** There is no `.pcof` yet,
    which is the required scratch-start state (pitfall 6.3).
  - `seed`: see §1.7. It is applied between `mixato` and `randomise`.
  - `settings` that belong to the `.pcof` (e.g. `rminpt`, `psigma2`) are kept in
    `epsr_project.json` and passed as overrides to the scratch run. **Verified:** an
    override of `rminpt 8.0` on a scratch run is written into the new `.pcof`.
  - `ereq` is the exception: the scratch start forces it to 0. A `.pcof` override of
    `ereq 2.5` came out as 0.0, both in memory and in the file. `create` should reject a
    non-zero `ereq` and point the user to `refine(ereq=...)`.
  - It writes a small `epsr_project.json` with the definitions, the seed and any pending
    settings, for provenance and `Project.load`. That file is ignored by EPSR.
- `project.run(iterations=None, **settings) -> Summary`: one invocation. On a project
  that has never run, it raises an error pointing to `equilibrate`.
- `project.equilibrate(iterations, *, scratch_iterations=1, record=False, **settings) -> History`:
  runs with `ereq 0`.
  - **First call on a new project (D2):** it runs the scratch start (`ireset 1`,
    `iinit 1`, `nsumt -1`, `niter = scratch_iterations`, plus the pending `.pcof`
    settings), then one run of `iterations − scratch_iterations`.
    `equilibrate(iterations=100)` therefore reproduces steps 7–8 of
    `validation/water/steps.txt` (1 + 99) byte for byte.
  - Later calls: one run of `iterations`.
  - "New project" means no `.pcof` in the folder, the same test the engine uses.
  - The returned `History` has one row per run (two on the first call), so the scratch
    iteration is visible.
  - `iterations` must be greater than `scratch_iterations` on the first call.
- `project.refine(ereq=None, runs=1, iterations=None, *, record=False, until=None, callback=None, snapshot_every=None, **settings) -> History`.
  Under the hood it calls the existing `epsr.refine`, with `ereq` applied on the first run
  only. That replaces `initial_overrides`.
- `project.set(**settings)`: a persistent edit of `.inp`/`.pcof` without running, which
  `validation/set_keyword.py` does today. Implement it as a no-iteration path, or as a
  file edit through the existing `KeywordFile` parser. **Decide which**; a file edit is
  simpler and matches `set_keyword.py`.
- `project.settings` returns the current values as `Settings`.
- `project.results() -> Simulation`, for reading only. Document the `.qdr` side effect, or
  add a read-only mode (§4.5).
- `project.samples` and `project.system`, loaded from `epsr_project.json` when present.
- `project.status()`: whether EPSR has run yet (a `.pcof` exists), iterations
  accumulated (`nsumt`), the current `ereq`, the last R-factor and the seed.
- `project.refine(...)` on a project that has never run raises an error pointing to
  `equilibrate`, because the reference-potential stage should be explicit.
- With `record=True`, `refine`, `run` and `equilibrate` use `epsr.vis` and return or
  attach the trajectory (`history.trajectory`, or `project.trajectory()`).

### 4.4 Changes to the existing syntax

Each item is backwards compatible: keep the old form and add or deprecate.

| # | current | proposed | why |
|---|---|---|---|
| 1 | `sim.fit(0)`, `data_fq(0)` (0-based) and `overrides={"data": {1: …}}` (1-based) | also accept a sample name: `sim.fit("D2O")`, `settings={"data": {"D2O": {"tweak": 1.05}}}`; names come from `epsr_project.json`, or else the data-file stem | one consistent, readable way to address datasets |
| 2 | `overrides=` (they persist, silently) | `settings=` as the preferred name, documented as persistent; `overrides` kept as an alias | the current name suggests a one-off change |
| 3 | `refine(..., initial_overrides={"ereq": 3})` | `project.refine(ereq=3, ...)` | a schedule vocabulary instead of a mechanism |
| 4 | `vis.run_recorded(sim, rec)`, `vis.refine_recorded(project, recorder=rec)` | `record=True` or `recorder=` on run, refine and equilibrate | one entry point; removes `if recorder is None` branching |
| 5 | `tools.run(folder, "mixato", answers=[…])` | typed wrappers (§4.3); `tools.run` stays | no positional prompts |
| 6 | `Simulation(folder)` used only to read results, which rewrites `.qdr` | `project.results()`, with the side effect documented or removed | reading should not modify files |
| 7 | `sim.potential.type_params` only | add `pair_params()`: ε_ij, σ_ij, q_i q_j, `coulfac`, the active `.pcof` extras and a `terms` column (`LJ + Coulomb` / `Coulomb only` / …) per pair, from new accessors | example 07 recomputes this by hand, and a zero ε silently turning a pair Coulomb-only is not visible anywhere today (§1.8) |
| 8 | peak and coordination analysis by hand | `epsr.analysis.first_peak(r, g, window)`, `coordination_number(...)`, re-exported from `epsr.vis.metrics` where possible | common post-processing |
| 9 | `PYTHON_WRAPPER_WORKFLOW.md` §4 says "`ereq` lives in the `.pcof`, so the `.pcof` file must be edited" | the overrides already accept `.pcof` keywords (see the `render_overrides` docstring); correct the doc | stale guidance |

### 4.5 Rust and binding additions needed

Keep them small. Accessors only; no change to engine behaviour.

1. **`_core.epsr_keywords()`** exposes `EPSR_SYMBOLS` plus its marks (`MARK_NSYMBOLUPTODATA`
   etc. in `symtab.rs`) and the `.pcof` split (`MARK_IEREQFIRST..LAST`, positions
   `MARK_IROVERLAP..` up to `MARK_INCOEFFP`). Mirror `tool_keywords` in `python/src/tools.rs`.
   **Confirm the exact `.pcof` routing** in `symtab.rs` (`setup_input_file` and the save
   path) before relying on positions.
2. **Pair-parameter accessors** in `api.rs`: `ecor`, `dcor`, `qcor` and `coulfac` per pair
   (`st.mol.*`, `st.pot.coulfac`). Also expose the per-pair `.pcof` extras: `rminex` and
   `pcoefex` (the repulsive core), and `npgauss`, `hpgauss`, `rpgauss` and `wpgauss` (the
   additional Gaussian or Morse terms), mapped from `.pcof` entries to pairs through
   `irefpcof`. Bind them in `python/src/sim.rs`. Add them to `_core.pyi`.
2a. **`_core.isotope_table()`** exposes `tools/sears91.rs` (`SEARS91`: symbol, mass
   number, scattering length) so that `NeutronSample` can validate isotopes before `nwts`
   runs. An unknown mass number makes `nwts` prompt on stdin instead of failing (§4.3).
3. (Optional) **A read-only setup**, so that reading results doesn't rewrite `.qdr`. Check
   where `.qdr` is written in `Epsr::setup_with_lines` (`epsr.rs:165` onwards) and whether
   a flag can skip it. Without one, `project.results()` can copy the `.qdr` and restore it.
4. Not needed: an in-memory, file-free `Epsr::from_config`. It would duplicate the parsers
   and risk byte-identity for no user-visible gain, since the files are hidden anyway.
   Revisit only if file I/O becomes a bottleneck.

---

## 5. Work plan (phased, each phase independently mergeable)

**Phase 0: prerequisites.**
- §4.5 items 1, 2 and 2a, plus the `set_rng_seed(ato_path, seed)` helper from §1.7.
- Seed tests:
  - `seed=None` and `seed=1` give byte-identical boxes.
  - `seed=42` is reproducible and differs from the default.
  - `seed=0` and negative seeds raise.
  - The state after `randomise` matches the §1.7 table (`mcall 6000`; `idum 1715014484`
    for the default, `1163647977` for seed 42).
- Tests: the keyword-table round trip (every default parses), and pair params against
  hand-computed Lorentz–Berthelot for EPSR-test.

**Phase 1: typed setup wrappers in `epsr.tools`.**
- `makemole` (from a `.mol` name), `mixato`, `randomise(folder, box, seed=None)` (the
  seed applied through the helper before the tool runs), `fmole` with the `.mol`
  handling, and `nwts` from structured components.
- Acceptance: rebuilding `validation/water` with only these calls gives a box and weights
  byte-identical to the validated run. Reuse `python/tests/epsr_workflow.py` to compare.

**Phase 2: `epsr.system`.**
- `AtomType`, `Molecule` (`render_mol`, `from_mol`, `from_xyz`, `from_ato`, `atom`) and
  `System`.
- Port the rendering logic from `archive/mellor/src/mellor/inputs.py` (`render_mol`,
  `mixato_answers`, the f32 and number formatting) and the data model from
  `archive/mellor/src/mellor/spec.py`. Both were working and tested; adapt the names.
- Acceptance: `render_mol` of the water molecule passes through `makemole` to an `.ato`
  identical to the validated `Water.ato`. Silica atoms do the same against
  `validation/silica_bulk`.

**Phase 3: `epsr.samples`.**
- `NeutronSample` and `.mint01` validation. Port `mellor.data.read_mint01`, `_abundances`
  and `render_nwts`.
- Acceptance: the generated `.NWTS.dat` passes through `nwts` to the same `.NWTStot.wts`
  as `validation/water`.

**Phase 4: `Settings` and `.inp` rendering.**
- Build `Settings` on `epsr_keywords()`. Port `mellor.inputs.render_inp`,
  `dataset_keywords` and `auto_rshmin` (the last as an opt-in default only).
- Acceptance: after the scratch run, the generated `.inp` equals the one EPSR rewrites for
  the validated water case.

**Phase 5: `Project` workflow.**
- `create`, `run`, `equilibrate` (with the scratch start, D2), `refine`, `set`,
  `settings`, `results`, `status` and `epsr_project.json`.
- Acceptance: `create` + `equilibrate(100)` + `refine(ereq=3.0, runs=2, iterations=100)`
  reproduce `validation/water` byte for byte.
- Dataset-by-name (§4.4 #1), `settings=` (#2), `record=` (#4) and `pair_params` (#7).

**Phase 6: example and docs.**
- Rewrite example 07 as in §4.2, reading its data from `experimental_data/water/data/`.
- Its `output/epsr_run/` must stay byte-identical to the committed one (the full schedule,
  not `--quick`). `python/tests/test_examples.py` only checks that the examples exit
  cleanly, so add a slow test that compares `epsr_run/` with the committed output.
- Update `python/docs/user_guide.md` (a new "Setting up a system" section), `api.md`,
  `python/README.md` and the top-level README snippet.
- Optionally port example 06 (silica) to the new API as a second acceptance case.

**Phase 7a: Cu-BTC, crystalline systems from a PDB file (C20–C23, §8).** Start only
after Phase 6 is complete and its byte-identity test passes.

**Phase 7b: later extensions, not now** (out of scope until agreed; see the checklist).
- **`PairTerms`: per-pair extras without editing files.** For example:

  ```python
  PairTerms(("OW", "HW"), coulfac=1.0,
            repulsive_core=RepulsiveCore(rmin=1.4, amplitude=...),
            extra=[Gaussian(r=1.8, width=0.1, height=-2.0)])
  ```

  It is passed to `Project.create(..., pair_terms=[...])` or `project.set_pair_terms(...)`.
  It is written into the per-pair blocks of the `.pcof`, so it is only possible after the
  scratch run. `create` stores it and the scratch step of `equilibrate` applies it, the
  same way as the other pending `.pcof` settings.
  - Covers only what EPSR supports (§1.8): `coulfac`, the repulsive core, and additional
    Gaussian or modified-Morse terms. **Not** per-pair LJ.
  - Fits a future table-driven force-field layer: a table row can map to `AtomType`
    values plus optional `PairTerms`.
  - Acceptance: EPSR-test-ex1 already exercises `rminex > 0` on three pairs, so use it as
    the oracle. Add the `.pcof` write path through the engine's own writer, not as Python
    text.
- X-ray samples
- Typed configs for the analysis tools (`partials`, `coord`, `triangles`, `sdf`) from
  `System` labels
- Non-cubic boxes
- `Project.load` from `epsr_project.json`

---

## 6. Risks and open questions

1. **Byte-identity of generated text.** `makemole` rewrites the `.mol` and writes the
   `.ato` from parsed values, so identity depends on the *values* reaching makemole (f32
   rounding), not on whitespace. The `.inp` is rewritten by EPSR after the first run. The
   main risk is the `potential`-line coarse-graining values and the geometry rounding.
   Mellor's `inputs.py` is tested byte-identical against
   `validation/` for `.mol` (water, silica) and `.NWTS.dat`
   (`archive/mellor/tests/test_inputs.py`), so start from it and port those tests too.
2. **The working directory and threads.** The setup tools chdir under a process-wide lock
   (`tools._WORKING_DIRECTORY_LOCK`). `Project.create` is therefore serial within a
   process. Document that parallel setups need processes. `refine_many` on prepared
   folders is unaffected.
3. **Determinism and seeds.** Resolved (D3): see §1.7. One thing stays open: whether
   `seed=` should also be offered on `Project.from_template` (re-seeding an existing
   box). It is not needed now.
4. **Box shapes.** `mixato` makes cubic boxes only. Non-cubic boxes come from
   `makelattice`/`readpdb`. Keep `System(box=None)` cubic and fail clearly otherwise
   until Phase 7.
5. **Unsupported Rust makemole features** (`dihedral`, `unitcell`). Raise
   `EpsrUnsupportedError` while building the `Molecule`, not inside the tool. Porting
   `makemole`'s `dihedral` line handling would unlock dihedral restraints, which the
   engine already evaluates (§1.8). It is worth a separate task if flexible molecules
   are needed.
5a. **Force-field import.** Per-bond force constants (Dissolve, OPLS/GAFF tables) have no
   EPSR equivalent (§1.8). Any importer must keep the geometry and the per-type LJ and
   charge values, drop the bond/angle constants with a warning, and set `vibtemp`,
   `angtemp` and `dihtemp` from EPSR defaults or user input. Per-pair LJ overrides
   (NBFIX-style) cannot be represented either.
6. **Where the scratch run goes.** Resolved (D2): in the first `equilibrate` call. What
   remains is to keep `equilibrate`'s docstring and the returned `History` explicit about
   the 1 + (N − 1) split.
7. **Naming.** Should the module names be `epsr.system` / `epsr.samples`, or should
   everything be exported at the top level (`epsr.AtomType`, …)? The recommendation is
   both: modules for organisation, re-exports in `epsr/__init__.py`.
8. **`epsr.vis` versus the standalone `epsrvis`.** `record=True` should use `epsr.vis`
   only. The three plotting generations (`visualisation/src/epsr_vis`, `epsrvis` and
   `epsr.vis`) are a separate cleanup and out of scope here.

## 7. Files to read before starting

- `python/epsr/{simulation,project,refine,potential,tools,io}.py`, `python/epsr/_core.pyi`
  and `python/src/{sim,tools}.rs`: the current API
- `python/examples/07_water_from_scratch/water_from_scratch.py`: the target to replace
- `validation/PYTHON_WRAPPER_WORKFLOW.md` (§6 pitfalls) and `validation/water/`: the
  pipeline and its oracle
- `rust/epsr-core/src/symtab.rs` (the keyword table), `symtab_tools.rs` (`NWTS`),
  `tools/makemole.rs` (the `.mol` grammar), `ato.rs:785-800` (the mixing rule) and
  `api.rs` (where to add accessors)
- `archive/mellor/src/mellor/{spec,inputs,data}.py`: working prior art for exactly this
  rendering layer
- `rust/docs/PYTHON_WRAPPER_NOTES.md`: the engine API contract
- For Phase 7a (Cu-BTC): `validation/cu-btc/` (`README.md`, `steps.txt`, `*.stdin`,
  `CuBTC.EPSR.inp`), `python/tests/test_cu_btc.py`, `python/tests/epsr_workflow.py`,
  `rust/epsr-core/src/tools/{readpdb,changeato}.rs` (supported options) and
  `archive/mellor/src/mellor/inputs.py` (`readpdb_answers`, `changeato_commands`)

---

## 8. Next example: Cu-BTC (example 08), specification

**Purpose.** This is the second development target after example 07, and it starts in
Phase 7a. Cu-BTC (HKUST-1) is the third validated real-data case. It is the first that
**isn't built from molecules by `mixato`**: the structure comes from a crystal (PDB
file), and EPSR computes Bragg peaks. Expressing it in the new syntax fills the gaps
for crystalline and pre-built configurations (C20–C23). The oracle already exists, so
byte-identity can be proved after every EPSR run.

### 8.1 The system (from `validation/cu-btc/`)

- **Structure:** one cubic unit cell, a = 26.3336 Å, 624 atoms (48 Cu, 192 O, 288 C, 96 H)
  at 300 K. `validation/cu-btc/make_pdb.py` converted the Dissolve species
  (`experimental_data/cu-btc/dissolve/input.txt`) into `validation/cu-btc/CuBTC.pdb`,
  splitting the framework at its Cu–O bonds into **24 Cu₂ paddlewheels** (residue `CU2`)
  and **32 benzene-tricarboxylate linkers** (residue `BTC`, 18 atoms). The crystal
  positions are kept: there is no `randomise` and no `fmole`.
- **Atom types**, numbered by `readpdb` in order of first appearance in the PDB; the
  labels are the PDB atom names:

  | # | label | element | ε (kJ/mol) | σ (Å) | charge (e) | mass |
  |---|---|---|---|---|---|---|
  | 1 | `CU` | Cu | 0.02092 | 3.144 | +1.098 | element default |
  | 2 | `O` | O | 0.401664 | 3.033 | −0.665 | default |
  | 3 | `C1` | C | 0.39748 | 3.473 | +0.778 | default |
  | 4 | `C2` | C | 0.39748 | 3.473 | −0.092 | default |
  | 5 | `C3` | C | 0.39748 | 3.473 | −0.014 | default |
  | 6 | `H` | H | 0.06276 | 2.846 | +0.109 | **1.008** (`atmass`) |

- **Bond detection:** `readpdb` with van der Waals factor **0.68** and no dihedrals. This
  finds exactly the Dissolve bonds plus the 24 Cu–Cu pairs. The 192 Cu–O bonds become
  intermolecular (LJ + Coulomb only), which is the main approximation of the model.
  `vanderWaalsRadii.txt` must be in the folder. It is byte-identical in `EPSR-test/` and
  `validation/cu-btc/`.
- **Monte Carlo parameters** (`changeato.stdin`): `temp 300`, `tol 0.0`, `stepmi 0.1`,
  `stepri 0.3`, `stepr 0.3`, `stepm 1.0`.
- **Data:** `experimental_data/cu-btc/data/Empty-CuBTC.mint01` (NIMROD; identical to the
  copy in `validation/cu-btc/`), natural isotopes, no exchangeable hydrogen, `nrtype 5`,
  `rshmin 0.7`.
- **Settings** that differ from the defaults: `qwidthqmax 0.2 5.0` (Bragg peaks up to
  5 Å⁻¹, Gaussian width 0.2 Å⁻¹), `hklqmin 0.0 10.0 0.0`, `qmin 0.05`.
- **Schedule** (`steps.txt`): scratch run (1) + 99 on the reference potential, then
  `ereq 10.0` and 2 × 100. That is `equilibrate(100)` then
  `refine(ereq=10.0, runs=2, iterations=100)` in the new API.
- **Expected results** (`validation/cu-btc/README.md`): R-factor 3.08×10⁻² → 2.29×10⁻²
  → 1.69×10⁻² → 1.64×10⁻²; 80 726 Bragg reflections at 940 unique Q; Cu–Cu 2.64 Å; Cu–O
  relaxes from 1.95 to 2.25 Å. The full run takes about 4 minutes.

### 8.2 Target syntax (wireframe)

```python
import epsr
from epsr.samples import NeutronSample
from epsr.system import AtomType, Framework, System          # Framework: 🆕 C20

CASE = REPOSITORY / "validation" / "cu-btc"                  # CuBTC.pdb lives here (made by make_pdb.py)
DATA = REPOSITORY / "experimental_data" / "cu-btc" / "data"

# Atom types in readpdb's order (first appearance in the PDB); labels = PDB atom names.
CU = AtomType("CU", element="Cu", epsilon=0.02092,  sigma=3.144, charge=1.098)
O  = AtomType("O",  element="O",  epsilon=0.401664, sigma=3.033, charge=-0.665)
C1 = AtomType("C1", element="C",  epsilon=0.39748,  sigma=3.473, charge=0.778)
C2 = AtomType("C2", element="C",  epsilon=0.39748,  sigma=3.473, charge=-0.092)
C3 = AtomType("C3", element="C",  epsilon=0.39748,  sigma=3.473, charge=-0.014)
H  = AtomType("H",  element="H",  epsilon=0.06276,  sigma=2.846, charge=0.109, mass=1.008)

framework = Framework.from_pdb(                              # 🆕 C20: readpdb
    CASE / "CuBTC.pdb",
    types=[CU, O, C1, C2, C3, H],      # checked against the PDB label order
    vdw_factor=0.68,                   # bond detection; finds Dissolve bonds + Cu-Cu
    dihedrals=False,
)
print(framework)        # 🆕 cubic cell 26.3336 Å, 624 atoms, 56 molecules (24 × 2 atoms, 32 × 18 atoms)

system = System("CuBTC", temperature=300.0, structure=framework)   # 🆕 density taken from the cell
system.monte_carlo(tol=0.0, stepmi=0.1, stepri=0.3, stepr=0.3, stepm=1.0)   # 🆕 C21, changeato names

samples = [NeutronSample("CuBTC", DATA / "Empty-CuBTC.mint01")]    # natural, no exchange

project = epsr.Project.create(
    OUTPUT / "epsr_run", system, samples,
    settings={"qwidthqmax": (0.2, 5.0), "hklqmin": (0.0, 10.0, 0.0),   # 🆕 C22 tuple values
              "qmin": 0.05, "rshmin": 0.7},
    overwrite=True,
)   # readpdb → changeato → nwts → .EPSR.inp; no randomise/fmole for a Framework

project.equilibrate(iterations=100)                          # scratch + 99, Bragg peaks from the start
history = project.refine(ereq=10.0, runs=2, iterations=100, record=record)

result = project.results()
bragg = epsr.io.read_bragg(project.output("bragg"))          # 🆕 C22: reflections, unique Q
print(result.potential.pair_params())                         # CU-O row: LJ + Coulomb (the Cu-O "bond")
```

Open design points for the implementer. Decide them in Phase 7a and record them in this
section:
- **Is `Framework` a separate class, or `System.from_pdb(...)`?** The recommendation is
  a separate `Framework`, because its molecules are *derived* by bond detection rather
  than declared. A `System` can then hold either `add(molecule, n)` entries or one
  `structure=`. Mixing them (molecules inside a framework) is out of scope.
- **`system.monte_carlo(...)` versus a `settings`-style dict.** These values live in the
  `.ato`, not the `.inp`, and are set through `changeato`. Keep EPSR's `changeato`
  keyword names (D1). The `changeato` keyword table is already exposed
  (`epsr.tools.keywords("changeato")`), so validate against it. The same call should
  also work for `mixato`-built systems, where today the defaults are used.
- **The `seed`** for a `Framework` is applied before the scratch run (§1.7), since
  nothing random happens before it.
- **The PDB file** is an input, like a `.mint01`. Converting Dissolve to PDB
  (`make_pdb.py`) stays a data-preparation script and is out of scope.

### 8.3 Gaps this example fills

| gap in today's wrapper (§2) | how Cu-BTC fills it | capability |
|---|---|---|
| Configurations not built by `mixato` (crystals, PDB structures) need `tools.run(…, "readpdb", answers=["0.68", "n"])` with positional prompts and a hand-copied `vanderWaalsRadii.txt` | `Framework.from_pdb(path, types, vdw_factor, dihedrals)`; the vdW radii file is shipped as package data and copied in; the label order is checked against `types` | C20 |
| Changing the force field or the MC step sizes of an existing `.ato` needs `changeato` with an interactive command script (`atom 1`, `epsilon …`, `save`, `y`, `e`) | `AtomType` values are applied by a typed `changeato` wrapper; `system.monte_carlo(...)` sets `temp`/`tol`/step sizes; `atmass` only when it differs from the element default | C21 |
| `.ato` Monte Carlo parameters (`tol`, `stepmi`, `stepri`, `stepr`, `stepm`) can't be set at all in the example 07 flow | `system.monte_carlo(...)`, also for molecular systems | C21 |
| Crystalline settings (`qwidthqmax`, `hklqmin`, `mplicities`) are raw strings in overrides | `Settings` with tuple values, validated against the keyword table | C22 |
| No readers for Bragg outputs (`.bragg`, `.inten`); the test helper `read_bragg_reflections` lives in `python/tests/epsr_workflow.py` | promote it to `epsr.io.read_bragg` (and `read_inten`) | C22 |
| Atom types whose label differs from the element (`CU`/Cu, `C1`–`C3`/C) and several types per element | exercised end to end (labels ≤ 3 chars) | C1, C20 |
| Density given implicitly by a cell rather than `density=` | `System(structure=framework)` takes the cell; `density=` is rejected when given together with a structure | C3, C20 |
| The setup path without `randomise`/`fmole` | `Project.create` chooses the pipeline from the system type; the seed is applied before the scratch run | C8, C7 |

**Still not covered after Cu-BTC:** non-cubic cells (`makelattice`, `unitcell`), X-ray
data, analysis-tool inputs, `PairTerms` (e.g. to model the Cu–O bond), and molecules
mixed into a framework. See the out-of-scope list.

### 8.4 Implementation notes and acceptance

- **Prior art:** `archive/mellor/src/mellor/inputs.py` (`readpdb_answers`,
  `changeato_commands` (`atmass` only when it differs from `EPSR_ELEMENT_MASSES`),
  `pdb_atom_labels`, `_check_pdb_atom_types`, `render_inp`) and
  `archive/mellor/tests/test_inputs.py`. Those tests already pin the answers and the
  `.inp` to the validated Cu-BTC files, so port them with the code. **Do not** port
  `render_vdw_radii`: its text differs from the validated file (`0` versus `0.0`, and the
  trailing newline). Ship `EPSR-test/vanderWaalsRadii.txt` as package data and copy it
  verbatim.
- **The existing oracle:** `python/tests/test_cu_btc.py` already runs `steps.txt` through
  the wrapper and compares every file with `fortran_reference.tar.gz` after each EPSR
  run, using `compare_dirs` and `extract_archive` from `epsr_cases.py`. Example 08's test
  must make the same comparison for a folder built through the new API.
- **Acceptance (C23):**
  - after `readpdb` and after `changeato`, the `.ato` is identical to the validated
    intermediate (generate it once with `run_steps_with_wrapper`)
  - after each of the four EPSR runs, every file matches the Fortran archive (same
    masking as `test_cu_btc.py`)
  - 80 726 reflections at 940 unique Q
  - the R-factors of §8.1
- **Location:** `python/examples/08_cu_btc/cu_btc.py` + `README.md`, and a slow test in
  `python/tests/` marked `experimental`/`slow` with `xdist_group("cu-btc")`, like the
  existing one.
