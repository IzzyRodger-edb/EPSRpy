# Benzene from scratch

EPSR refinement of liquid benzene against three NIMROD measurements: C₆H₆, C₆D₆ and a
50:50 mixture of the two (`experimental_data/benzene/data/*.mint01`, ISIS runs 27749,
14901 and 27731). The only input files are the measured data and the molecule's
coordinates, [`benzene.xyz`](benzene.xyz). [`benzene_from_scratch.py`](benzene_from_scratch.py)
defines everything else in Python with the setup API. The force field, density and
EReq are those of the Dissolve benzene example. Unlike examples 07–09, this system has
**no validated EPSR reference**: it shows the setup API on a new system.

| input | defined in the script as | written by |
|---|---|---|
| reference potential (OPLS-AA ε, σ, charge per atom type) | `AtomType("CA", element="C", ...)`, `AtomType("HA", element="H", ...)` | the `potential` lines of `Benzene.mol` |
| molecule geometry and restraints | `Molecule.from_xyz("benzene.xyz", types={"C": CA, "H": HA})`: bonds detected from covalent radii, the 18 angles generated | `Benzene.mol` → `makemole` |
| box (200 molecules, 0.876 g/cm³ = 0.0810 atoms/Å³, 300 K) | `System("benzene", ...).add(benzene, 200)`, density converted from g/cm³ in the script | `mixato`, `randomise`, `fmole` |
| neutron weights | `NeutronSample(..., isotopes={"HA": {2: 1.0}})` (C₆D₆), `{"HA": {1: 0.5, 2: 0.5}}` (50:50) | `<sample>.NWTS.dat` → `nwts` |
| EPSR control file | `settings={"rshmin": 0.7}` | `benzene.EPSR.inp`, every other keyword at EPSR's default |

`benzene.xyz` is an ideal D6h ring with the gas-phase bond lengths C–C 1.397 Å and
C–H 1.084 Å. The 120° angle restraints keep the ring planar, because a hexagon whose
six angles are all 120° cannot pucker.

Run from the `python/` folder (the `plot` extra gives the figures, `viewer` the 3D image):

```bash
uv run python examples/10_benzene_from_scratch/benzene_from_scratch.py            # ~10 minutes
uv run python examples/10_benzene_from_scratch/benzene_from_scratch.py --quick    # 10 iterations per run
uv run python examples/10_benzene_from_scratch/benzene_from_scratch.py --seed 42  # another starting box
```

With `epsr-py` installed by pip instead (`pip install "./python[plot]"` from the
repository root, see the [top-level README](../../../README.md#1-installation)), run
`python python/examples/10_benzene_from_scratch/benzene_from_scratch.py` from the root.

The schedule follows the other liquids: `project.equilibrate(iterations=100)` (scratch
start + 99 iterations on the reference potential), then
`project.refine(ereq=20.0, runs=2, iterations=100)`.

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the EPSR folder (`.ato`, `.pcof`, `.EPSR.inp`, `.erg`, every `.EPSR.*`
  result), plus `Benzene.mol` and `epsr_project.json` (the Python definitions), which
  EPSR does not read.
- `trajectory/`: every iteration, recorded by `epsr.vis` (reload with `vis.Trajectory.load`).
- `figures/`: `01`–`09` and `summary.csv`/`summary.md` from `vis.make_report`; `10` the
  ring-ring structure of the final configuration (the g(r) of the ring centres, and the
  angle between ring planes against the centre-centre distance: 0° stacked, 90°
  T-shaped); `11` the 3D configuration; `12` the reference and empirical potentials;
  `fits.png` the final fit to each dataset.
- `summary.txt`: the console summary (system, pair parameters, R-factors, peaks, ring orientation).

The R-factor falls from 2.0×10⁻³ (reference potential) to 1.4×10⁻³ (empirical
potential). By sample: C₆H₆ 2.2×10⁻⁴, C₆D₆ 5.3×10⁻⁴, 50:50 3.4×10⁻³. The ring centres
have their first peak at 5.7 Å. Of the rings 5–7.5 Å apart, about half are near
perpendicular (> 60°) and 15 % near parallel (< 30°).

## Limitations

- **The 50:50 sample is fitted worst.** Its data rise steeply below Q ≈ 0.3 Å⁻¹, and the
  model does not reproduce this (see `fits.png`). The weights also treat H and D as mixed
  at random over every H site, while the sample is a mixture of whole C₆H₆ and C₆D₆
  molecules. The two differ only in the intramolecular H–H weights. Modelling the mixture
  exactly, with two species whose H atoms are different types (`HA` and `DA`, with
  identical potentials), was tried and fitted worse. The empirical potential then refines
  the `HA` and `DA` pairs separately, and the two physically identical species partly
  demix in the box.
- **Interlocked rings.** `randomise` places molecules at random, so in the starting box
  some rings pass through each other (78 pairs at the start). The Monte Carlo separates
  most of them within 50 iterations, but 5 pairs stay interlocked, with an H atom through
  the neighbouring ring. Their repulsion makes the reported energy meaningless (+13 000 kJ/mol
  per molecule at the end). The script counts these pairs, reports them and
  leaves them out of the orientation statistics; they are the spike near 1.3 Å in
  figure `10`. A softer Lennard-Jones repulsion during equilibration (`power` 8 or 10
  instead of 12) does not free them, and keeping `mixato`'s box (`randomise=False`)
  stacks every molecule at the same place. Removing them needs an overlap-free insertion
  step, which the setup API does not have yet.
