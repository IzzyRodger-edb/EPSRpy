# Perturbing the potential

Three ways to change the potential of a copy of the DMSO/water test case
(`EPSR-test/`) from Python. [`perturb_potential.py`](perturb_potential.py) uses the
water hydrogen bond (atom types `O1`–`H1`) as its example pair:

| way | in the script | what it changes |
|---|---|---|
| (a) EPSR keywords | `project.set(ereq=15.0, feedback=0.8, potfac=(1.0, 1.0))` | edits the `.pcof` / `.EPSR.inp` files, as a user would between runs; read back with `project.settings` |
| (b) the EP coefficients, in memory | `simulation = project.simulation()`, then `simulation.potential.set_ep_coefficients(...)`, `perturb_ep(pair, value)`, `scale_ep(factor)`, then `simulation.run()` | damps the O1–H1 empirical potential by half, adds a constant to its coefficients and scales the whole EP; the run refines the modified EP further and saves it |
| (c) a fixed EP | `project.run(potfac=(0.0, 0.0))` | Monte Carlo on the current potential, with refinement switched off; `potfac` persists, so set it back to refine again |

The in-memory route (b) is the one place where a `Simulation` is used directly. For
everything else, `Project` reads and writes the folder. The script also checks that
total = reference + empirical, and lists the reference-potential parameters of each atom
type (`potential.type_params`). These come from the `.ato` file and are read-only here.

`template/` holds the files the run needs, copied from `EPSR-test/`: the box
`200DMSO400Water.ato`, the potential `200DMSO400Water.pcof`, the control files
`.EPSR.inp` and `.EPSR.inpa`, the accumulated results that EPSR keeps averaging
(`.EPSR.g01`, `.f01`, `.s01`, `.y01`), the energy record `.EPSR.erg`, the three data
files `D6DMSO_*.mint01` and their weights `DMSODWater*.NWTStot.wts`.

Run from the `python/` folder:

```bash
uv run python examples/03_perturb_potential/perturb_potential.py    # ~10 s
```

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the copy of `EPSR-test/` after the two runs.

The console shows the O1–H1 EP minimum before (−10.3 kJ/mol) and after (−5.1 kJ/mol)
damping, the R-factor and EP amplitude after the run, and confirms that the run with
`potfac 0` leaves the EP unchanged.
