# The auxiliary EPSR programs

EPSR comes with 35 auxiliary programs besides the refinement: setup programs
(`makemole`, `mixato`, `randomise`, `fmole`, `nwts`, `readpdb`, `changeato`...) and
analysis programs (`partials`, `coord`, `triangles`, `sdf`, `sharm`...). `epsr.tools`
runs all of them through the Rust port, on a folder, as EPSR does.
[`auxiliary_tools.py`](auxiliary_tools.py) runs them on a copy of the DMSO/water test
case (`EPSR-test/`), in each of the three ways they take their input:

| kind of program | in the script | example |
|---|---|---|
| keyword-file setup program | `tools.nwts(folder, "DMSODWaterD.NWTS.dat")` | neutron weights from an existing `.NWTS.dat` |
| keyword-file analysis program | `tools.keywords("partials")`, `tools.write_keyword_file(path, "partials", {...})`, `tools.partials(folder, file)` | the input file written from a dict (keywords checked), then the partial g(r) around water O, read back with `epsr.read_table` |
| command-line arguments | `tools.run(folder, "writexyz", [...])` | the whole box as `.xyz` |
| questions on standard input | `tools.run(folder, "writexyz", ["DMSO.ato"], answers=[...])` | one answer per line, as typed at the prompt |

Errors are Python exceptions: a missing input file raises `FileNotFoundError`.

The setup programs also have typed wrappers that `Project.create` uses: `tools.makemole`,
`tools.mixato`, `tools.randomise`, `tools.fmole`, `tools.readpdb`, `tools.changeato`.
Examples 05 and 06 call them directly on existing EPSR input files.

`template/` holds the files the run needs, copied from `EPSR-test/`: the box
`200DMSO400Water.ato`, the potential `200DMSO400Water.pcof`, the control files
`.EPSR.inp` and `.EPSR.inpa`, the accumulated results that EPSR keeps averaging
(`.EPSR.g01`, `.f01`, `.s01`, `.y01`), the energy record `.EPSR.erg`, the three data
files `D6DMSO_*.mint01` and their weights `DMSODWater*.NWTStot.wts`. The tools also
use the weights input `DMSODWaterD.NWTS.dat` and the DMSO template `DMSO.ato`.

Run from the `python/` folder:

```bash
uv run python examples/04_tools/auxiliary_tools.py    # a few seconds
```

## Output

Everything is written to `output/` (or `$EPSR_EXAMPLES_WORK/output/`):

- `epsr_run/`: the copy of `EPSR-test/`, with the files the programs wrote:
  `DMSODWaterD.NWTStot.wts`, `water.PARTIALS.dat` and its results `water.PARTIALS.g01`...,
  `200DMSO400Water00001.xyz` and `DMSOappend.xyz`.
