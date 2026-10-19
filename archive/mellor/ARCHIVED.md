# Mellor: archived

**Status:** archived on 2026-10-19. Mellor is no longer part of the active codebase;
development now centres on the Python wrapper (`python/`) around the Rust engine
(`rust/`). The code here is kept complete and working, so that an agent can be
reinitialised from it later.

- **Last commit with Mellor at the top level:** `f1fe713` ("Put EPSR Mellor logo onto
  webapp"). `git show f1fe713:mellor/<path>` gives any file as it was then.
- **Checked when archived:** after `uv sync` in this folder, the fast test suite
  passes from here (`uv run pytest -m "not slow and not gpu and not network"`:
  235 passed, 1 skipped because torch was not installed), and so do the `slow`
  end-to-end EPSR runs (`-m "slow and not gpu and not network"`: 12 passed).

Start with [README.md](README.md) (usage), [docs/architecture.md](docs/architecture.md)
(design) and [docs/integration.md](docs/integration.md) (the Python API and run-folder
contract that the web app used).

---

## 1. What Mellor is

An LLM agent that takes measured structure-factor data (Gudrun `.mint01` files), proposes
an atomistic model as a JSON `SystemSpec`, generates the EPSR inputs and a standalone
`run_<name>.py`, runs that script through the `epsr` wrapper, analyses the result with
`epsrvis` and writes a PDF report. Afterwards it answers questions about the run through
tools. Design rules worth keeping in any reimplementation:

1. The model only produces the specification and prose. Every number comes from the engine.
2. The generated script is the single source of truth. Mellor executes it, so a user can
   re-run a result with no LLM.
3. The library of validated systems (water, silica, DMSO/water, benzene, Cu-BTC) is
   used whenever the data match one of them. Anything else is flagged as unvalidated.
4. A `scripted` backend runs the whole pipeline with no model.

LLM backends (`src/mellor/llm/backend.py`):

| backend | needs | default model |
|---|---|---|
| `transformers` (local) | `uv sync --extra local` (torch, transformers), ideally a GPU | `Qwen/Qwen3-4B-Instruct-2507` |
| `hf-api` | `HF_TOKEN` | set with `--model`, e.g. `Qwen/Qwen3-235B-A22B-Instruct-2507` |
| `scripted` | nothing | none |
| `auto` | picks one of the above; `MELLOR_BACKEND` and `MELLOR_MODEL` override | |

## 2. Layout

```
archive/mellor/
├── ARCHIVED.md            this file
├── README.md              usage
├── pyproject.toml, uv.lock
├── src/mellor/
│   ├── spec.py            SystemSpec (the JSON the model fills in)
│   ├── library.py         validated system templates (numbers copied from validation/ and EPSR-test/)
│   ├── inputs.py          SystemSpec -> EPSR input files (.mol, .NWTS.dat, .EPSR.inp, *.stdin)
│   ├── script.py          SystemSpec -> standalone run_<name>.py (calls epsr + epsrvis)
│   ├── runner.py          executes the script; run-folder manifest
│   ├── data.py            .mint01 reader and validation
│   ├── analysis.py, figures.py, report.py, export.py   metrics, figures, PDF report, tidy tables
│   ├── sources.py         local / Dissolve / Zenodo search for comparable data
│   ├── cli.py             `mellor run | chat | ...`
│   └── llm/               backends, prompts, agent loop, tool definitions
├── tests/                 unit tests; `slow` ones run EPSR end to end
└── docs/
    ├── architecture.md, integration.md
    ├── DATA_ACQUISITION_FEASIBILITY.md   spec for a data-fetching agent (ISIS, ORNL, ...); never built
    └── EPSR_Mellor_logo.png
```

## 3. How it connects to the rest of the repository

These are the points to check first when bringing Mellor back, because they are where
the rest of the repository can drift from it.

**Packages, by relative path** (`pyproject.toml` `[tool.uv.sources]` and `uv.lock`):

- `epsr-py` → `../../python`
- `epsrvis` → `../../visualisation/python_visualisation_toolkit`

**Wrapper API used** (almost all of it from the generated script, `src/mellor/script.py`):

- `epsr.Simulation(folder, inp_name, overrides=..., verbose=...)` and its accessors
- `epsr.tools.run(folder, command, args, answers, verbose=...)` for the setup programs
- `epsr.tools.keywords(name)` (`src/mellor/inputs.py`)
- `epsr.__version__`
- `epsrvis`: `Trajectory`, `record_run`, `make_report`, `metrics.fit_metrics`,
  `metrics.first_shell`, `plots.fq_fit`, `plots.final_outputs`, `outputs.OutputFolder`,
  `systems.get_case`, `report.FIGURES`

The wrapper now also has `epsr.vis`, which grew from the same code as `epsrvis`. A
reimplementation could depend on `epsr.vis` alone and drop the `epsrvis` dependency.

**Repository data** (found through `repository_root()` in `library.py` and `sources.py`):

- `validation/{water,silica_bulk,cu-btc}/`: library templates. Cu-BTC reads
  `validation/cu-btc/CuBTC.pdb` at run time.
- `EPSR-test/`: the DMSO/water template
- `experimental_data/<system>/data/*.mint01`: test data and local data search
  (`LOCAL_DATA_FOLDERS` in `sources.py`)

## 4. What changed when it was archived

| before | after |
|---|---|
| `mellor/` | `archive/mellor/` |
| `mellor/mellor_unseen_verification/{C6H6,C6D6,5050}.mint01` (NIMROD benzene data, for testing Mellor on unseen data) | `experimental_data/benzene/data/` |
| `experimental_data/<system>/{input.txt,CMakeLists.txt}` (Dissolve inputs; Cu-BTC also `7108574.cif`) | `experimental_data/<system>/dissolve/` |
| `docs/EPSR_Mellor_logo.png` | `archive/mellor/docs/` |
| `DATA_ACQUISITION_FEASIBILITY.md` | `archive/mellor/docs/` |
| `webapp/` (FastAPI front end around Mellor) | deleted |

Path edits made so that Mellor keeps working from its new location:

- `pyproject.toml`, `uv.lock`: `../python` → `../../python`, `../visualisation/...` → `../../visualisation/...`
- `src/mellor/library.py` `repository_root()`: `parents[3]` → `parents[4]`
- `src/mellor/sources.py` `repository_root()`: now returns the nearest ancestor that holds
  the data folders, instead of the parent of `mellor/`
- `tests/{specs,example_run,test_data,test_export,test_llm}.py`: `parents[2]` → `parents[3]`
- `README.md`, `docs/integration.md`: relative paths updated

**The web app.** `webapp/` was deleted, not archived. Recover it with
`git checkout f1fe713 -- webapp`. The parts that matter for an agent are:

- `webapp/src/epsrweb/mellor_adapter.py`: the only module that called Mellor
- `webapp/tests/test_mellor_contract.py`, `test_mellor_source.py`: tests that pinned
  every assumption the app made about Mellor's API
- `webapp/scripts/make_logo_assets.py`: generator for the logo

## 5. Reinitialising

**Run it as it is, from the archive:**

```bash
cd archive/mellor
uv sync                     # builds the epsr Rust extension on first sync
uv run pytest -m "not slow and not gpu and not network"    # fast suite
uv run pytest -m slow       # end-to-end EPSR runs (minutes)
uv run mellor run ../../experimental_data/water/data/*.mint01 --out work/water --system water --backend scripted
```

**Move it back to the top level:** `git mv archive/mellor mellor`, then undo the path
edits listed in §4 (one directory level fewer in each).

**Update it against a changed wrapper:**

1. Run the fast suite. A change to `epsr.Simulation`, `epsr.tools` or `epsrvis` shows up
   first in `tests/test_script.py` and `tests/test_inputs.py`, and then in the `slow` tests.
2. Update `src/mellor/script.py`. It is a template for the generated `run_<name>.py`, so
   changes to the wrapper API almost all land there.
3. If the data folders were renamed, update `LOCAL_DATA_FOLDERS` in `sources.py` and the
   paths in `library.py` (`validation/cu-btc/CuBTC.pdb`).
4. To use `epsr.vis` instead of `epsrvis`, change the imports in `script.py`,
   `analysis.py` and `figures.py`, and drop `epsrvis` from `pyproject.toml`.

**Open items when it was archived:**

- Mellor's choices (proposed system, potentials, schedule) had been checked only on its
  five library systems. The benzene data in `experimental_data/benzene/` were added to
  test it on unseen data, but no results of that test are recorded.
- The data-acquisition agent in `docs/DATA_ACQUISITION_FEASIBILITY.md` was specified but
  never built.
