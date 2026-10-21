# Plan: EPSR_python_translation

A line-by-line Python translation of the EPSR26 Fortran suite, in a new repository
`EPSR_python_translation`. The point is readability: being able to modify the code and see how and
when each routine is called. Speed is not a goal. The user-facing contract is `epsrshell`. It should
look and behave like the Fortran shell, with two fixes: `ls`/`pwd`/`cd` work, and the shell shuts
down cleanly with no residual processes.

The repository is created by [`setup_epsr_python_translation.sh`](setup_epsr_python_translation.sh)
(section 9). This file is copied into the new repository as `PLAN.md`.

---

## 1. Goals and non-goals

**Goals**

1. Translate every Fortran routine that `epsrshell` can reach into Python. That is the `epsr` engine,
   the 35 setup and analysis programs, the shell itself, and `setup` and `plot`. One Python function
   per Fortran subroutine/function, with the same name and in the same order of operations.
2. Produce output **byte-identical** to the Fortran on the three real-data validation cases:
   water (H₂O/D₂O/HDO), amorphous SiO₂ and Cu-BTC. The same standard as the Rust port.
3. Keep the `epsrshell` experience: prompt, banner, help text, command names, script mode
   (`runflag.txt`, `ss`/`ps`/`es`) and the questions each program asks on stdin.
4. Make the call flow visible: optional tracing of every routine call (section 5).

**Non-goals**

- Speed. No vectorisation, no numba/Cython, no parallelism. A faster path can come later, but the
  literal translation stays as the reference.
- PGPLOT graphics (`plot2d`, `plot3d`, `plotato`). Not in the first pass. `plot` (gnuplot scripts)
  and `plot3djmol`/`runjmol` (Jmol) are in scope.
- Changing physics or file formats. The five bugs fixed in the Rust port (`docs/FORTRAN_BUGS.md` §1)
  are fixed the same way here and are flagged in the code. Nothing else diverges.

## 2. Translation rules

The Rust port shows that byte-identical output is possible, and that it comes from the rules below.
Translators should keep `fortran_src/<file>` and the matching Rust module (EPSRpy
`rust/epsr-core/src/...`, mapped in `docs/rust_port/FORTRAN_MAP.md`) open side by side. Every
mixed-precision cast has already been worked out in the Rust code.

**Structure**

- One Python module per Fortran source file. Files with a `program` go in `epsr26/programs/`; shared
  routine files go in `epsr26/library/`. Each function's docstring starts with the Fortran
  location, e.g. `"""molfit.f:molfit — one Monte Carlo cycle."""`.
- Variable names are the Fortran names. `goto` becomes loops and flags, with the original label kept
  as a comment (`# 3000 continue`).
- COMMON blocks, `.inc` files and Fortran module variables become module-level objects in
  `epsr26/common/`, one per include file (`common.molecom.xa`). This is global state, as in
  the Fortran. `common.reset()` recreates all of it before each program runs, because every Fortran
  program started in a fresh process.
- `stop` raises `FortranStop(code, message)`. The shell catches it.

**Numbers**

- `real` → `numpy.float32`, `real*8`/`double precision` → Python `float`, `integer` → Python `int`
  (`numpy.int32` only where 32-bit wrap-around matters). Requires NumPy ≥ 2: under NEP 50, a Python
  literal no longer promotes a `float32` expression to double.
- Single-precision literals and constant expressions must be explicit. Fortran folds `1.0/3.0` in
  single precision, Python in double. Use the helper `r4(...)`.
- Mixed precision is promoted exactly as gfortran does it. Follow the Rust casts.
- Arrays are NumPy arrays with one extra element per dimension. Index 0 is unused so that `a[i, j]`
  is Fortran `a(i,j)`. Use `order='F'`. Where the Fortran relies on sequence association (passing
  `a(1,j)` as a 1-D array), pass the view `a[:, j]` explicitly.
- **No NumPy reductions or vectorised expressions in translated code.** `np.sum` uses pairwise
  summation and SIMD transcendentals, which change the last bit. Loops stay loops, in the Fortran
  order.
- Single-precision intrinsics (`sinf`, `cosf`, `expf`, `logf`, `powf`, `acosf`, `atan2f`, …) call
  glibc `libm` through `ctypes` (`epsr26/fortran/fmath.py`). NumPy's float32 ufuncs do not use glibc.
  Double precision uses `math.*`, which already calls glibc. The gfortran helpers `__powisf2`
  (`x**i`), `nint`, `sign`, `aint`, `anint` and constant powers (`x**2` → `x*x`) are ported from
  `fortran/fmath.rs`.
- Integer semantics: Fortran `mod` and integer division truncate toward zero, while Python `%` and `//`
  floor. Use `fmod_i`/`idiv` helpers and never the bare operators on values that can be negative.
- Reference build: gfortran `-O2 -fno-tree-vectorize`, scalar math, `OMP_NUM_THREADS=1`. glibc
  matters: the Cu-BTC references had to be regenerated after a glibc update (EPSRpy
  `docs/plans/cu_btc_fortran_discrepancy.md`). Record the glibc version with every reference set.

**Text I/O**

- `epsr26/fortran/` holds a FORMAT interpreter (F, E, ES, G, I, A, X, /, Hollerith, `1p` scale
  factor), gfortran list-directed output (`write(6,*)` spacing), list-directed input, blank-padded
  fixed-length strings, `backspace`/`rewind`, and a unit registry including the default `fort.NN`
  files. These are ports of the Rust `fortran/format.rs`, `listread.rs`, `text.rs` and `output.rs`
  (≈1,600 lines), and everything else depends on them. Do them first and test them on their own.

## 3. Repository layout

```
EPSR_python_translation/
  PLAN.md  README.md  CLAUDE.md  LICENSE  pyproject.toml
  src/epsr26/
    fortran/      compatibility layer: fmath, fformat, listio, text, arrays, units, FortranStop
    common/       one module per .inc / COMMON block / Fortran module, plus reset()
    library/      shared routine files (parse.py, molfit.py, ato_read_and_save.py, ...)
    programs/     one module per Fortran program (epsr_standalone_rev3.py, nwts.py, mixato.py, ...)
    shell/        epsrshell.py (translation), builtins.py (ls/pwd/cd), processes.py, signals.py
    trace.py      call tracing (section 5)
  fortran_src/    the Fortran, unchanged; ACTIVE_SOURCES.txt lists the files being translated
  EPSR-test/      the original 200 DMSO + 400 water test case
  experimental_data/{water,silica_bulk,cu-btc}/
  validation/{water,silica_bulk,cu-btc}/ + run_validation.sh, set_keyword.py, plot_gr_agreement.py
  test_problems/  115 feature cases, each with commands.txt and fortran_reference.tar.gz
  tools/          build_fortran_reference.sh, fortran_link_lists.txt, compare_outputs.py, gen/
  docs/           FORTRAN_BUGS.md, rust_port/ (architecture, Fortran map, tools, test cases)
  tests/          pytest: compatibility layer, shell (pexpect), ex cases, validation (slow marker)
```

The package is named `epsr26`, not `epsr`, so it can be installed next to `epsr-py`. The console
script is `epsrshell`. If the Fortran `epsrshell` is on the same `PATH`, use `python -m epsr26`.

## 4. The shell

### 4.1 What is preserved

`epsr26/shell/epsrshell.py` is a translation of `epsrshell_standalone.f` (960 lines). It keeps:

- the banner (`EPSRshell> Welcome to EPSR version ...`), the `EPSRshell>` prompt, the `help`/`?` text,
  the "Command ... not recognised -- Will attempt to run command anyway" message;
- every command name and alias in `process_command` (`setup` → `upset`, `readgpr`/`readjmol` →
  `readmole`, `mapgr` → `mapgrquad`, `bye`/`e`, ...);
- script mode: `runflag.txt`, `ss`/`ps`/`es`, `#` comments, rewinding the script at the end, and
  pausing/resuming across shell restarts (`start_script`, `check_script_status`, `pause_script`,
  `end_script`);
- `system_commands.txt` for external tools (gnuplot, Jmol, editor);
- line editing and history, using the stdlib `readline` with the history file the Fortran
  `FCreadline` writes.

What changes underneath: programs no longer run as separate binaries through `call system(...)`.
`run_program` looks the command up in a table (`programs/__init__.py`), calls `common.reset()`, and
runs the program in the shell's own process. Its stdin is the terminal (or the script), as before.
Running in-process is what makes both fixes below possible.

### 4.2 Fix 1: `ls`, `pwd`, `cd` (and `md`, `del`, `edit`)

In the Fortran, these commands build a command line from `system_commands.txt` (`system_ls`,
`system_cd`, `system_termination`) and pass it to `call system`. The template shipped with EPSR26
(`EPSR-test/system_commands.txt`) holds the Windows values (`dir/b`, `cd/d`, `\`). `cd` only changes
the shell's internal `directory_name` string. I have not reproduced the Linux failure. These
Windows values are the likely cause, so check it against the Fortran shell in phase 3.

The translation implements them directly, without `system_commands.txt`:

| command | behaviour |
|---|---|
| `pwd` | prints `os.getcwd()` |
| `cd [dir]` | `os.chdir`. `cd`, `cd home` and `cd .` return to the start folder (Fortran behaviour); `cd ..` and `cd -` work as in bash |
| `ls [args]` | with no arguments, lists the folder like `ls`. With arguments, runs `/bin/ls args` in the foreground |
| `md`, `del` | `os.makedirs`, `os.remove` (with glob) |
| `edit file` | `$EDITOR`, then `system_edit`, then `nano` |
| `system cmd ...` or `!cmd ...` | runs `cmd` through `/bin/sh` in the foreground |

Unknown commands still print the Fortran message. They then run through `/bin/sh` instead of from
`system_binaries`.

### 4.3 Fix 2: clean shutdown and no residual processes

| situation | behaviour |
|---|---|
| `bye`, `e`, `exit`, `quit`, Ctrl+D at the prompt | runs the shutdown sequence below, exit status 0 |
| Ctrl+C at the prompt | discards the line and shows a new prompt (as bash does) |
| Ctrl+C during `epsr` | the current iteration becomes the last: `write_results` runs and all files are saved consistently, then back to the prompt. A message says so |
| second Ctrl+C within 3 s, or Ctrl+C in a setup/analysis program | aborts the program with `KeyboardInterrupt` and returns to the prompt with a warning that its files may be incomplete |
| Ctrl+C while a script runs | as above, then the script is paused (as `ps`), so `runflag.txt` is never left saying `running` |
| SIGTERM, SIGHUP (terminal closed) | the shutdown sequence |

Shutdown sequence (`shell/signals.py`, registered with `atexit` and the signal handlers):

1. Save the readline history.
2. If a script is running, pause it in `runflag.txt`.
3. Terminate every child process the shell started. Every external process goes through
   `processes.spawn()`, which starts it in its own process group and records it. Each group gets
   SIGTERM, then SIGKILL after 3 s. Viewers (gnuplot, Jmol) are stopped too, unless
   `system_keep_viewers yes` is set in `system_commands.txt`.
4. Remove the lock file `.epsrshell.lock`.

The lock file holds the shell's PID and is created in the start folder. At startup, a lock owned by
a live process gives a warning ("another EPSRshell is running in this folder: PID n"). Two shells
running the same script are a known way to get stray runs. A `jobs` built-in lists the
live children.

Tests: `tests/shell/` drives the shell with `pexpect` and checks the prompt, `ls`/`pwd`/`cd`, Ctrl+C at the
prompt, during a long `epsr` and during a script. After each test, no process from the shell's
process group may remain (`psutil`).

## 5. Seeing how and when routines are called

- `@routine` decorator (`trace.py`) on every translated function. It has no effect unless tracing
  is on.
- `EPSR_TRACE=calls epsrshell`, or `trace on [pattern]` / `trace off` in the shell, writes an
  indented call tree to `trace.log` with iteration numbers, e.g. `iter 3 > molfit > movem >
  en_change`.
- `trace summary` prints the call count and cumulative time of each routine for the last program
  run.
- `tools/gen_call_graph.py` builds the static call graph from the Fortran `call` statements, into
  `docs/call_graph.md`, to compare against the traced one.

## 6. Validation

All comparisons use `tools/compare_outputs.py` (from EPSRpy). It masks only the date on line 2 of
`.out` and the timing columns of `fort.35`. Everything else must be byte-identical.

**References.** `tools/build_fortran_reference.sh` builds the Fortran programs (gfortran is
installed). `validation/make_reference.sh CASE` runs each `steps.txt` with the Fortran only and
saves the folder after every step as `references/CASE/stepNN/`, plus `references/CASE/glibc.txt`.
These references are generated locally, not committed (a few hundred MB). A SHA-256 manifest of each
step is committed, so a mismatch is visible without the Fortran.

**Two modes**, because the Python will be orders of magnitude slower than the Fortran:

- *Resume mode* (the main test): each step starts from the Fortran folder after the previous step,
  runs one step in Python, and compares. Each step is tested on its own, and a 100-iteration step can
  be cut to `niter 2` (with a matching Fortran reference) without running the 99 before it. This
  makes the EP stage testable without hours of reference-potential iterations.
- *Chain mode* (final acceptance): the full `steps.txt` in Python from the raw inputs, compared after
  every step, exactly as `run_validation.sh` does for Rust. Run once per case before calling the
  translation complete, overnight if needed.

Before writing `steps_short.txt` for each case, measure the time per iteration in phase 1.

**Coverage**

| case | exercises |
|---|---|
| water | `makemole`, `mixato`, `randomise`, `fmole`, `nwts` ×3, `epsr` scratch start / reference potential / EP / restart |
| silica_bulk | the same for an atomic network (two atom templates, `mixato` with two components) |
| cu-btc | `readpdb`, `changeato`, `nwts`, `epsr` with Bragg peaks (`.bragg`, `.inten`) on every iteration |
| EPSR-test | the original DMSO/water example, as a quick smoke test |
| test_problems ex1–ex115 | one case per feature and per tool (X-ray, triclinic boxes, three-body, coordination constraints, all analysis programs) |

The known Rust workaround applies here too: `fmole` with the `.mol` file present runs `makemole`
internally (validation README, issue 1). The Python should port that path properly. The water case
can then drop its `delete Water.mol` step, with a new Fortran reference.

## 7. Phases

Each phase ends with its tests passing on CI (GitHub Actions: ruff + pytest, without slow tests).

| # | phase | content | done when |
|---|---|---|---|
| 0 | repository | run the setup script; pyproject, CI, lint; build the Fortran reference; `make_reference.sh` + resume-mode runner | references for the three cases exist and the runner reports "not implemented" for every step |
| 1 | compatibility layer | `fortran/`: fmath (libm via ctypes), FORMAT, list-directed I/O, strings, units, arrays, `ran1`; time a float32 inner loop to estimate run times | unit tests against strings and values produced by small gfortran programs; `ran1` stream equals Fortran's |
| 2 | shared library | include files → `common/`; `parse`, `stringutility`, `file_and_directory_commands`, `system_commands`, `read_and_save_files`, `setup_symbols_generic`, `process_symbols`, `ato_read_and_save`, `ato_manipulation_rev4`, `minimagec`, `bondnumbers`, `get_molecules`, `reallocation_routines`, `lngamma_mod` | `.ato`, `.EPSR.inp` and `.pcof` files read and rewritten byte-identically for all reference folders |
| 3 | shell | translate `epsrshell_standalone.f`; built-ins; signal handling; process registry; program table with stubs | pexpect tests of section 4 pass; no residual processes; the user tries it on a real folder |
| 4 | setup programs on the validation path | `makemole`, `mixato`, `randomise` (introtcluster.f), `fmole`, `nwts`, `readpdb`, `changeato` | resume mode: all setup steps of water, silica and Cu-BTC identical |
| 5 | engine | `epsr_standalone_rev3`, `molfit`, `cellneighbourlist`, `ergpairreference`, `setup_epsr_rev1`, `pfroutines_mod`, energies, structure, potential, weights, Bragg (`bragginten` path), then `threebody_mod`, `coordconstraint_mod`, `containerforces_mod` | resume mode: every `epsr` step of the three cases identical (shortened `niter`); then chain mode once per case |
| 6 | remaining programs | the other 26 tools (`partials`, `coord`, `triangles`, `torangles`, `clusters`, `chains`, `rings`, `lsi`, `sharm`, `sdf`, `sdfcube`, `mapgr`, ...), `setup` (setup_generic.f), `plot` (gnuplot), `plot3djmol`/`runjmol`, `writexyz`/`readxyz`, `xwts`/`epsrwts`/`epsrwtsx` | all 115 ex cases identical |
| 7 | documentation | `docs/call_graph.md`, a per-routine index (Fortran file:line ↔ Python module), user README | — |

Phases 3 and 4 can run in parallel once phase 2 is done. The shell needs only the program table.

## 8. Risks and open questions

- **Run time.** A rough guess is 100–1000× slower than the Fortran, which is minutes per EPSR
  iteration for the water box. Phase 1 measures it. Resume mode is the mitigation. If an
  interactive run is unusably slow, a later fast path (numba on hot loops, checked against this
  translation) could help, but that is outside this plan.
- **glibc dependence.** Byte identity of `sinf`/`expf` holds only on the glibc used for the
  references. Regenerate the references on each machine and do not compare across machines.
- **gfortran list-directed and G-format output** have edge cases (field widths around 0.1 and 1e7,
  `-0.0`). The Rust `format.rs` has already solved them, so port it rather than reinvent it.
- **Single process.** The EPSRpy validation notes that the setup programs depend on the current
  directory. The shell `chdir`s, so this is fine for interactive use, but two folders cannot run in
  one process at the same time. That is acceptable for a shell.
- **`fmole` → `makemole` inside one process.** `fmole` calls `makemole`, which then shares the RNG
  state. Running in-process matches the Fortran here, because the Fortran also did this within one
  program.
- To decide: should the new repository be public? The setup script creates it private.

## 9. Files copied from EPSRpy

Copied from the committed `HEAD` with `git archive`, so uncommitted files are not included.

| EPSRpy | new repository | why |
|---|---|---|
| `fortran_src/` (all of it) | `fortran_src/` | the source being translated. Older revisions are kept for reference; `ACTIVE_SOURCES.txt` (generated) lists the ~75 files in scope |
| `EPSR-test/` | `EPSR-test/` | original test case, smoke test |
| `experimental_data/{README.md,water,silica_bulk,cu-btc}` | same | measured data and the Dissolve inputs the force fields came from |
| `validation/{water,silica_bulk,cu-btc}/`, `run_validation.sh`, `set_keyword.py`, `plot_gr_agreement.py` | `validation/` | the three validation cases (inputs, `steps.txt`, stdin answers, Cu-BTC `fortran_reference.tar.gz`) and the harness |
| `validation/README.md` | `validation/RUST_PORT_RESULTS.md` | expected results (R-factors, energies) and known issues |
| `test_problems/` (185 MB, largest file < 20 MB) | `test_problems/` | 115 feature cases with Fortran references, for phase 6. Skip with `INCLUDE_EX_CASES=0` |
| `rust/tools/{compare_outputs.py, build_fortran_reference.sh, build_patched_readxyz_readpdb.sh, fortran_link_lists.txt, run_ex_test.sh, gen/gen_symtabs.py}` | `tools/` | reference build, file comparison, ex-case runner, keyword-table generator (paths fixed by the script) |
| `rust/docs/{ARCHITECTURE,FORTRAN_MAP,TOOLS,TEST_CASES,VALIDATION}.md` | `docs/rust_port/` | routine map and the lessons of the Rust port |
| `FORTRAN_BUGS.md`, `LICENSE` | `docs/FORTRAN_BUGS.md`, `LICENSE` | known Fortran bugs and how they are handled; MIT licence |

Not copied: `rust/` sources (consult them in EPSRpy), `python/`, `archive/`, `visualisation/`,
`internal_play/`, `validation/api_driver/`, `experimental_data/benzene/`.
