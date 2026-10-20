# A refinement loop

Repeated EPSR runs on a copy of the DMSO/water test case (`EPSR-test/`), with
convergence monitoring, as the EPSR shell loop does. [`refinement_loop.py`](refinement_loop.py)
uses `project.refine`, the same call that examples 05–10 use after equilibration:

| option of `project.refine` | in the script | what it does |
|---|---|---|
| `runs` | `--runs` (default 10) | number of EPSR runs |
| EPSR keywords | `ntimes=5` | Monte Carlo cycles per iteration; persists in the files |
| `callback` | `report_progress(run_number, simulation)` | called after every run (R-factor, energy, pressure, EP amplitude); returning true stops the loop |
| `until` | `target_reached(history)` | stops early once the R-factor is below 10⁻⁴ |
| `snapshot_every` | `5` | copies the folder to `output/epsr_run.snapshots/run0005/`... |

`ereq` is not given, so the files' value (10 kJ/mol) is kept. The returned
`epsr.History` has one row per run, and `project.history()` reads the folder's whole
`.erg` record (every iteration since the test case was first made).

`template/` holds the files the run needs, copied from `EPSR-test/`: the box
`200DMSO400Water.ato`, the potential `200DMSO400Water.pcof`, the control files
`.EPSR.inp` and `.EPSR.inpa`, the accumulated results that EPSR keeps averaging
(`.EPSR.g01`, `.f01`, `.s01`, `.y01`), the energy record `.EPSR.erg`, the three data
files `D6DMSO_*.mint01` and their weights `DMSODWater*.NWTStot.wts`.

Run from the `python/` folder:

```bash
uv run python examples/02_refinement_loop/refinement_loop.py            # 10 runs, ~30 s
uv run python examples/02_refinement_loop/refinement_loop.py --runs 3
```

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the copy of `EPSR-test/` after the runs.
- `epsr_run.snapshots/`: a copy of the folder every five runs.
- `history.csv`: the `History` table (R-factor, energies, pressure, EP amplitude, per-dataset R-factors).
- `history.png`: the same, plotted by `epsr.plot.history` (needs the `plot` extra).

The test case is already well refined, so the R-factor stays near 1.15 × 10⁻⁴ and the
target of 10⁻⁴ is not reached: all the runs are made.
