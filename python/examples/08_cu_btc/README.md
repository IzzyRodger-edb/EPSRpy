# Cu-BTC: a crystalline framework from a PDB file

EPSR refinement of the metal-organic framework Cu-BTC (HKUST-1) against the NIMROD data
of the empty framework (`experimental_data/cu-btc/data/Empty-CuBTC.mint01`), with Bragg
peaks calculated up to Q = 5 Å⁻¹. [`cu_btc.py`](cu_btc.py) defines the system with the
setup API; unlike example 07 the configuration comes from a crystal structure, not from
molecules:

| input | defined in the script as | written by |
|---|---|---|
| the unit cell (624 atoms, a = 26.3336 Å) | `Framework.from_pdb("CuBTC.pdb", types=[...], vdw_factor=0.68)` | `readpdb`: one molecule per PDB residue, bonds by van der Waals radii |
| force field (ε, σ, charge, mass per type) | `AtomType("CU", element="Cu", ...)` ... in the PDB's label order | `changeato` |
| temperature and Monte Carlo step sizes | `System("CuBTC", 300.0, structure=framework).monte_carlo(tol=..., stepm=...)` | `changeato` (EPSR's keyword names) |
| neutron weights | `NeutronSample("CuBTC", ...)` (natural isotopes) | `nwts` |
| Bragg peaks and data treatment | `settings={"qwidthqmax": (0.2, 5.0), "hklqmin": (0.0, 10.0, 0.0), ...}` | `CuBTC.EPSR.inp` |

The PDB file is `validation/cu-btc/CuBTC.pdb`, made from the Dissolve example by
`validation/cu-btc/make_pdb.py`: the framework is split at its Cu–O bonds into 24 Cu₂
paddlewheels and 32 benzene-tricarboxylate linkers. The 192 Cu–O bonds are therefore
intermolecular in this model (Lennard-Jones + Coulomb only), its main approximation;
see `validation/cu-btc/README.md`.

Run from the `python/` folder:

```bash
uv run python examples/08_cu_btc/cu_btc.py            # ~4 minutes
uv run python examples/08_cu_btc/cu_btc.py --quick    # 10 iterations per run
```

The schedule is that of `validation/cu-btc/steps.txt`: `equilibrate(iterations=100)`
(scratch start + 99 on the reference potential, Bragg peaks from the start), then
`refine(ereq=10.0, runs=2, iterations=100)`. The run folder is `output/epsr_run/`.

The R-factor falls from 3.1×10⁻² (scratch start) to 2.3×10⁻² (reference potential) and
1.6×10⁻² (empirical potential); EPSR enumerates 80 726 Bragg reflections at 940 unique
Q values (`epsr.io.read_bragg`). The test `tests/test_framework.py` (slow) runs the same
schedule through the API and requires every file to be byte-identical to the validated
steps (`validation/cu-btc/steps.txt`) after each of the four EPSR runs.
