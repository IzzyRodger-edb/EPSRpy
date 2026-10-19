# Briefing: Cu-BTC results differ from the archived Fortran reference

**Status:** open, 2026-10-20. Diagnosed but not resolved. Nothing in the engine, the
archive or `validation/` has been changed.
**For:** the agent that investigates and resolves it. Read this first; it records what
is established, with the commands to reproduce each fact.

## The symptom

`python/tests/test_cu_btc.py::test_matches_fortran_after_every_run` (an existing test, on
`main` before the setup-API work) and the new
`python/tests/test_framework.py::test_cu_btc_equals_the_fortran_after_every_run` both fail
on this machine. Both compare the run folder, after each of the four EPSR runs of
`validation/cu-btc/steps.txt`, with `validation/cu-btc/fortran_reference.tar.gz`.

| after step | files that differ from the archive |
|---|---|
| 4 (scratch start, 1 iteration) | `CuBTC.EPSR.o01`: 133 of 441 lines, in the 7th significant digit |
| 6 (+99 iterations, reference potential) | the same, plus last digits of `.erg` (4 lines), `.inpa`, `.out`, `fort.35` |
| 9 (+100 iterations, `ereq` 10) | as above; one `.out` value differs by one unit in its last printed digit (`Av. intermolecular pressure 0.787812E+01` vs `0.787813E+01`), which `epsr_cases.compare_file` classes as **different**, so the test fails |
| 10 (+100 iterations) | last digits of `.erg` (40 lines), `.inpa`, `.o01`, `.terg`, `fort.35`; all within tolerance |

Every other file (`.ato`, `.pcof`, `.bragg`, `.inten`, the g(r) and F(Q) tables...) is
identical after every step.

## What is established

1. **The Rust engine equals the Fortran built on this machine, byte for byte, after
   every EPSR run.** Both differ from the archive identically, line by line. Run with the
   Fortran build of 2026-10-17 at `~/epsr-work/fortran-build` (made by
   `rust/tools/build_fortran_reference.sh`: `-O2 -fno-tree-vectorize`, "scalar math";
   gfortran 15.2.0, glibc 2.43; the three files that differ from `fortran_src/` are the
   script's documented patches). The script and its output:
   `~/epsr-work/cubtc-diagnosis/` (folders `fortran/`, `rust/`, `archive/`), made by
   running `validation/cu-btc/steps.txt` with both programs:

   ```
   step  4/6/9/10  rust vs fortran-now: identical=True
   step  4/6/9/10  fortran-now vs archive: the differences in the table above
   ```

2. **The Python layer adds nothing.** The setup API (`Project.create` + `equilibrate` +
   `refine`), the old wrapper route (`tests/epsr_workflow.run_steps_with_wrapper`) and
   the `epsr-rs` CLI give byte-identical folders after every run
   (`tests/test_framework.py` checks the first two at each run; the CLI was compared
   after the full schedule).

3. **Same-machine comparisons are identical for all three cases; only cross-machine
   comparisons differ.** On this machine, on 2026-10-17, the local Fortran and Rust gave
   identical files at every step for water and silica (`~/epsr-work/val-water`,
   `~/epsr-work/val-silica`), and the setup API reproduces those Fortran folders exactly
   (water through the committed example 07 output; silica re-run on 2026-10-20). The
   Cu-BTC archive is the only reference that was **made on another machine**: it was
   committed in `d1130b2` (2026-10-17, author BradleyMak) by
   `validation/cu-btc/make_fortran_reference.py`.

4. **The first difference is in the reference-potential table.** In `.o01` after the
   scratch start, values differ in the 7th significant digit across the whole r range
   (e.g. CU-related columns at r = 0.15 Å: `0.9574267E+13` vs `0.9574269E+13`; at
   r = 11.91 Å: `-0.1404794E-05` vs `-0.1404793E-05`). The table is computed by
   `generate_rp` (`rust/epsr-core/src/potential.rs`), whose Lennard-Jones term uses
   single-precision `expf`/`sinhf`/`powf`. `python/docs/testing.md` already records the
   same effect for the EPSR-test `.o01` hash: "their 7th significant digit depends on the
   libm rounding of `powf`, which differs between glibc versions". The Monte Carlo is
   chaotic, so a last-digit difference in the potential appears in the energies after the
   first iterations; it stays in the last digits here.

## Working hypothesis

The archive was produced with a Fortran build whose single-precision libm results differ
from this machine's in the last bit for some arguments: a different glibc, a different
gfortran, or a vectorised build (`build_fortran_reference.sh --vectorize`, which makes
gfortran call glibc's SIMD math). The Rust port is faithful: it reproduces the local
Fortran exactly. The archive is therefore a platform-specific reference, not a portable
one.

## How to confirm or refute it

1. Ask for, or find, the environment `d1130b2` was made on (OS, glibc, gfortran, flags of
   `build_fortran_reference.sh`). If it is glibc ≠ 2.43 or a vectorised build, the
   hypothesis is confirmed.
2. Rebuild the local Fortran with `--vectorize` (into a new folder) and re-run the
   comparison script of the appendix with that build. If it reproduces the archive, the
   cause is the build mode.
3. Isolate the function: evaluate `expf`, `sinhf` and `powf` for the arguments
   `generate_rp` uses at the first differing r of `.o01` (r = 0.15 Å, the CU-O pair) in a
   tiny C program and compare the bits with the archive's value. `rust/epsr-core` calls
   libm through `f32` methods, so its results follow the platform's libm, as the
   Fortran's do.

## How it could be resolved (each needs the user's agreement)

- **Make the reference portable by comparing on the same machine:** when a Fortran build
  is available (`EPSR_FORTRAN_BUILD`), run it and compare with it, as
  `test_experimental_data.py::test_full_schedule_matches_fortran` already does for water
  and silica; fall back to the archive with the tolerance otherwise.
- **Keep the archive, accept last-digit differences:** treat a one-unit difference in the
  last printed digit as "close" (the `.out` line of step 9 is the only "different"
  result), and keep requiring byte-identity between the Rust CLI and the Python wrapper.
- **Regenerate the archive on the machine the tests run on**, with
  `make_fortran_reference.py`, and record its environment in `validation/cu-btc/README.md`.
  This changes a validated input, so it needs the user's explicit approval.

Do not change `rust/epsr-core` to match the archive: it already matches the Fortran on
this machine, and the engine's byte-identity with the Fortran is the project's main
guarantee.

## Appendix: the comparison script

```python
"""Run validation/cu-btc/steps.txt with the local Fortran build and with the epsr-rs CLI;
compare each after every EPSR run with fortran_reference.tar.gz and with each other.
usage (from python/): uv run python fortran_cubtc.py ~/epsr-work/cubtc-diagnosis"""
import shutil, sys
from pathlib import Path
sys.path.insert(0, "tests")
from epsr_workflow import (CU_BTC, load_steps, prepare_case_folder, run_step_externally,
                           cli_command, fortran_command)
from epsr_cases import compare_dirs, extract_archive, cli_binary
root = Path(sys.argv[1]); shutil.rmtree(root, ignore_errors=True); root.mkdir(parents=True)
archive = extract_archive(CU_BTC.case_folder / "fortran_reference.tar.gz", root / "archive")
build = Path.home() / "epsr-work" / "fortran-build"
folders = {"fortran": prepare_case_folder(CU_BTC, root / "fortran"),
           "rust": prepare_case_folder(CU_BTC, root / "rust")}
commands = {"fortran": fortran_command(build), "rust": cli_command(cli_binary())}
for number, step in enumerate(load_steps(CU_BTC), start=1):
    for name, folder in folders.items():
        run_step_externally(step, folder, CU_BTC.case_folder, commands[name])
    if step.is_epsr_run:
        reference = archive / f"after_step_{number:02d}"
        for label, (a, b) in {"fortran-now vs archive": (reference, folders["fortran"]),
                              "rust vs archive": (reference, folders["rust"]),
                              "rust vs fortran-now": (folders["fortran"], folders["rust"])}.items():
            comparison = compare_dirs(a, b)
            print(f"step {number} {label}: identical={comparison.all_identical}")
            if not comparison.all_identical:
                print(comparison.report())
```
