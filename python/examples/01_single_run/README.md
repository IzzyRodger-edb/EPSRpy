# A single EPSR run

One EPSR run on a copy of the DMSO/water test case (`EPSR-test/`: 200 DMSO + 400 water,
three SANDALS datasets), followed by a look at the results. `EPSR-test/` is a finished
EPSR folder with a box, a potential and weights, so EPSR can continue from it straight
away. [`single_run.py`](single_run.py) shows the calls that every later example builds on:

| step | in the script | what it does |
|---|---|---|
| open a folder | `epsr.Project.from_template("template", "output/epsr_run")` | copies the folder (EPSR rewrites its files on every run) and opens the copy |
| inspect it | `project.status()` | iterations so far, `ereq`, last R-factor, the datasets |
| run EPSR | `project.run(iterations=1, ntimes=5)` | one run, with EPSR keywords by their own names (checked, and saved to the files) |
| read the results | `project.results()` | a read-only `Simulation` of the folder: grids, `fit(name)`, `gr_inter(a, b)`, `potential`, `positions` |
| analyse | `epsr.analysis.first_peak(result.r, result.gr_inter("O1", "H1"), (1.5, 2.3))` | the water hydrogen-bond distance |

`template/` holds the files the run needs, copied from `EPSR-test/`: the box
`200DMSO400Water.ato`, the potential `200DMSO400Water.pcof`, the control files
`.EPSR.inp` and `.EPSR.inpa`, the accumulated results that EPSR keeps averaging
(`.EPSR.g01`, `.f01`, `.s01`, `.y01`), the energy record `.EPSR.erg`, the three data
files `D6DMSO_*.mint01` and their weights `DMSODWater*.NWTStot.wts`.

Run from the `python/` folder:

```bash
uv run python examples/01_single_run/single_run.py    # a few seconds
```

With `epsr-py` installed by pip instead, run
`python python/examples/01_single_run/single_run.py` from the repository root.

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the copy of `EPSR-test/` after the run.

The console shows the R-factor, energy and pressure of the run, the fit of each dataset
(`D6DMSO_D2O`, `D6DMSO_HDO`, `D6DMSO_H2O`, named after their data files), the
intermolecular O1–H1 peak at 1.80 Å (water–water hydrogen bond) and the O–H1 peak at
1.65 Å (DMSO oxygen to water hydrogen).
