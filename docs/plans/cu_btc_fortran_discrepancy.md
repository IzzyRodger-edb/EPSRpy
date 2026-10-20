# Briefing: Cu-BTC results differ from the archived Fortran reference

**Status:** resolved, 2026-10-20: the cause is glibc's `sinhf` (see Findings), and the
archive was regenerated on glibc 2.43 (see Resolution). The engine is unchanged.
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

## Findings (2026-10-20)

**Cause:** the archive was made on a machine whose glibc is older than 2.41. glibc 2.41
replaced the single-precision `sinhf` (the fdlibm code, up to about 1 ulp of error) with
CORE-MATH's correctly rounded version, and 2.42 fixed its last rounding bug (BZ 32627).
`generate_rp`'s Lennard-Jones term calls `sinhf`. Substituting the pre-2.41 `sinhf` for
this machine's (glibc 2.43) makes **both the local Fortran and the Rust engine identical
to `fortran_reference.tar.gz` after every EPSR run of `steps.txt`**. No other change is
needed. The build mode, the inputs, `powf`, `expf` and `logf` play no part.

Work files are in `~/epsr-work/cubtc-libm/`: `old_sinhf.c` (the preloadable `sinhf`),
`trace.c`, `mvec.c`, `cmp_args.c`, the comparison-script variants and the run logs.
The vectorised build is in `~/epsr-work/fortran-build-vectorize/`, its run in
`~/epsr-work/cubtc-vectorize/`. The old-`sinhf` wrappers are in
`~/epsr-work/fortran-build-oldsinhf/`, that run in `~/epsr-work/cubtc-oldsinhf/`.

### 1. A vectorised build does not reproduce the archive (confidence: high)

`rust/tools/build_fortran_reference.sh ~/epsr-work/fortran-build-vectorize --vectorize`
(the arguments are `BUILD_DIR [--vectorize]`). The vectorised build calls libmvec only
for `sinf` (`epsr_standalone_rev3.o`, `lsi.o`) and for `expf`/`exp`
(`auxiliary_routines.o`): `nm -u *.o | grep _ZGV`. `generate_rp`'s loop is not
vectorised: it calls `chsph`, and there is no `_ZGV…sinhf`/`logf`. Running the appendix
script with this build (`build = …/fortran-build-vectorize`):

```
step 4 fortran-now vs archive: identical=False
  CuBTC.EPSR.d01  close 142 line(s)   CuBTC.EPSR.o01  close 133 line(s)
  CuBTC.EPSR.q01  close 123 line(s)   CuBTC.EPSR.r01  close 228 line(s)   CuBTC.EPSR.t01 close 102 line(s)
step 9 fortran-now vs archive: identical=False
  CuBTC.EPSR.out  different line 21: 'Av. intermolecular energy -0.537743E+03…' vs '-0.538882E+03…'
  (and .d01 .erg .f01 .g01 .inpa .inten .l01 .p01 .q01 .r01 … "different")
```

The vectorised build has the same 133 `.o01` lines as the scalar build, plus differences
in the transforms (libmvec `sinf`). The Monte Carlo then diverges macroscopically. The
archive has none of these, so it was made with a scalar build.

### 2. The first difference is one `sinhf` call per value, 1 ulp apart (confidence: high)

Correction to "What is established" §4: the r = 0.15 Å entries `0.9574267E+13` /
`0.9574269E+13` are in the C1-C1, C1-C2, C1-C3, C2-C2, C2-C3 and C3-C3 columns, not CU
columns. The step-4 `.o01` differences are spread over every pair (4 to 41 rows each)
and go both ways. CU-O first differs at r = 0.18 Å.

The LJ term for a C-C pair (σ = 3.473 Å, ε = 0.39748 kJ/mol, `power` 12, `refpotfac` 1)
at x = 0.15 Å, in the Fortran's float operation order (`trace.c`):

```
ddif = dc/x = 23.1533356     logf -> 3.14213896 (404918ce)
pwr1 = 9*xx clamped to 20    expf(20)     = 485165184  (4de75844)
pwr2 = 3*xx = 9.4264164      exact sinh   = 6205.98364743
  sinhf glibc 2.43  = 6205.98389 (45c1efdf)    correctly rounded
  sinhf glibc<=2.40 = 6205.9834  (45c1efde)    -1 ulp
sumlj = 2*fac*sinhf(pwr2)*expf(pwr1):
  2.43  -> 550b52ee = 9574268534784 -> printed 0.9574269E+13   (this machine)
  <=2.40-> 550b52ed = 9574267486208 -> printed 0.9574267E+13   (archive)
```

The Coulomb term (order 10²) is below one float ulp at 10¹³ (about 10⁶), and
`wfac` = `wfacc` = 1 here. `logf`, `expf` and `powf` (the last only in `prefac`) give the
same bits on both glibc versions: their implementations have not changed since glibc
2.28 (`/usr/share/doc/libc6/NEWS.gz`: the CORE-MATH imports of 2.41 list `sinhf`, `coshf`,
`tanhf`, `expm1f`… but not `expf`/`logf`/`powf`). Scalar against libmvec on this machine
(`mvec.c`): `logf` and `expf` agree. libmvec's `sinhf` gives 45c1efde for this argument,
but across all 4400 distinct `sinhf` arguments of step 4 it differs from the old scalar
`sinhf` in 1215. Preloading it gives 149 `.o01` lines different from the archive.

**The decisive test.** `old_sinhf.c` is an `LD_PRELOAD` `sinhf` implementing glibc
2.27–2.40's `e_sinhf.c` + `s_expm1f.c` (the fdlibm code, written from the algorithm
because fetching the glibc source was not possible here). Its constants match their
documented bit patterns. Over 2·10⁶ random arguments it is off by 1 ulp for 10 %,
against 0 for glibc 2.43. The Rust CLI imports `sinhf@GLIBC_2.43` dynamically (`nm -D`),
so the preload replaces it in both programs. Step 4 alone:

```
system sinhf: 18480 calls, 4742 results differ from the old sinhf -> .o01 DIFFERS (133 lines)
old sinhf:                                                       -> .o01 identical to archive
```

The whole schedule, with the appendix script (`fortran_command(…/fortran-build-oldsinhf)`,
`cli_command(…/fortran-build-oldsinhf/epsr-rs)`, both wrappers exporting
`LD_PRELOAD=libold_sinhf.so`):

```
step 4 / 6 / 9 / 10   fortran-now vs archive: identical=True
step 4 / 6 / 9 / 10   rust vs archive:        identical=True
step 4 / 6 / 9 / 10   rust vs fortran-now:    identical=True
```

**Corroboration: the EPSR-test hash.** `python/docs/testing.md` (BradleyMak, 43d1c14,
2026-10-17) says `cargo test --release` reports `200DMSO400Water.EPSR.o01` as
differing and blames `powf`. On this machine the test **passes**. With the old `sinhf`
preloaded into the test binary it fails on exactly that file:

```
system sinhf: test result: ok. 1 passed
old sinhf:    1 of 26 files differ from the Fortran reference: ["200DMSO400Water.EPSR.o01"]
```

`reference_hashes.txt` was made by IzzyRodger-edb (c3e5322) with a glibc ≥ 2.41
`sinhf`. The Cu-BTC archive was made by BradleyMak with an older one. The same effect
explains both, and the `testing.md` paragraph names the wrong function and has the
direction backwards for this machine.

### 3. Nothing else contributes (confidence: high)

- `CuBTC.pdb`, `readpdb.stdin`, `changeato.stdin`, `Empty-CuBTC.mint01` and
  `vanderWaalsRadii.txt` in `after_step_04/` are byte-identical to `validation/cu-btc/`.
  `CuBTC.NWTS.dat` and `CuBTC.EPSR.inp` differ from the committed copies only because
  `nwts` and `epsr` rewrite them. The local run's rewritten files equal the archive's
  (`.NWTS.dat`, `.EPSR.inp`, `.pcof`, `.NWTStot.wts` after step 10 all `cmp` equal).
- Since d1130b2, `validation/cu-btc/` changed only `README.md`, `gr_agreement.png`,
  `make_pdb.py` (the `.pdb` it made is unchanged) and `make_fortran_reference.py`
  (refactoring and error reporting in b616baa; the same steps, environment and copy
  rule). `fortran_src/`, `build_fortran_reference.sh`,
  `build_patched_readxyz_readpdb.sh`, `fortran_link_lists.txt` and `set_keyword.py`
  are unchanged.
- The tar members belong to `bradleymak/bradleymak` (written 2026-10-17 13:26–13:55
  +0100), not `izzyr`, which confirms another machine. Nothing in the archive records
  its OS, glibc or compiler.
- The `sinhf` substitution alone reproduces every file at every step (§2), which leaves
  no room for another cause.

### What could not be determined

- The exact environment of the archive machine. The evidence fixes it to a glibc whose
  `sinhf` is the fdlibm code (glibc 2.27–2.40, e.g. Ubuntu 22.04/24.04 or Debian 12),
  with a scalar build. Not 2.41, whose CORE-MATH `sinhf` would differ from 2.43 in only a
  few arguments. Not the vectorised build (§1). Asking BradleyMak for `ldd --version`
  would confirm it. It is presumably x86-64; an aarch64 libm with the same `sinhf`
  results cannot be excluded.
- Whether `expf`/`logf`/`powf`, which are not correctly rounded, would differ on another
  architecture or another libm. On x86-64 glibc ≥ 2.28 they agree.

### 4. Recommended resolution

The archive is not wrong, and neither is the engine: the two Fortran results differ only
by glibc's `sinhf`. glibc ≥ 2.42's `sinhf` is **correctly rounded**, so any machine with
it gives the same bits, now and in future versions. That makes it the right definition
of the reference. Recommended, in order (each needs the user's agreement):

1. **Regenerate the archive on glibc ≥ 2.42** with the scalar build
   (`make_fortran_reference.py ~/epsr-work/fortran-build …`) and record in
   `validation/cu-btc/README.md` the environment (glibc ≥ 2.42, gfortran 15.2, scalar
   `build_fortran_reference.sh`) and why it matters (`sinhf`). This is the briefing's
   third option, with the reason it is portable rather than only local. It makes the
   archive consistent with `reference_hashes.txt` and with the water and silica
   comparisons, which already assume a correctly rounded `sinhf`. The `--vectorize`
   result shows the build mode must also stay scalar.
2. **Make the archive tests say why they fail on an old glibc** instead of failing on
   last digits: check `sinhf(9.4264164f)` (bits 0x4116d29a → expected 0x45c1efdf) through
   `ctypes` and skip (or xfail) with "needs glibc ≥ 2.42 (correctly rounded sinhf)". The
   same check could gate the `.o01` warning in `test_reference.py` and the Rust hash test.
3. **Optionally also compare with a local Fortran build when `EPSR_FORTRAN_BUILD` is
   set** (the first option), as `test_experimental_data.py` does. This works on any
   machine, but it is a complement, not a replacement: without a local build, the
   archive is still needed.
4. **Fix `python/docs/testing.md`** (the `powf` paragraph): the function is `sinhf`, the
   hashes match glibc ≥ 2.41/2.42, and older glibc is the case that reports the `.o01`
   mismatch. Also correct "CU-related columns" in §4 above.

Not recommended: **treating a last-digit difference as "close"**. It would hide the
step-9 `.out` line, but the difference grows through the chaotic Monte Carlo. Here it
stays in the last digits by step 10 only by luck; the vectorised run shows how far a
one-ulp change can go (intermolecular energy −537.7 vs −538.9 kJ/mol after step 9).
Relaxing the tolerance would also weaken the byte-identity check that is the project's
main guarantee, while the cause is fully explained and avoidable.

## Resolution (2026-10-20, approved by the user)

The archive machine's environment is now known: `python/docs/testing.md` (43d1c14,
BradleyMak) records its results as "Linux (WSL2), glibc 2.39, … 17 October 2026".
glibc 2.39's `sinhf` is the fdlibm code, as §2 inferred.

- `validation/cu-btc/fortran_reference.tar.gz` regenerated with
  `make_fortran_reference.py ~/epsr-work/fortran-build` (scalar build, glibc 2.43,
  gfortran 15.2.0). It has the same members as before. It differs from the old archive
  exactly as in "The symptom" above. Its `after_step_10/` is identical to the Rust CLI's
  folder from the diagnosis run. The environment is recorded in
  `validation/cu-btc/README.md`.
- `python/tests/epsr_cases.py`: `HAS_CORRECTLY_ROUNDED_SINHF` checks the platform's
  `sinhf` (looked up the way the engine resolves it) on the argument of §2.
  Without a correctly rounded `sinhf`, `test_cu_btc.py::test_matches_fortran_after_every_run`
  skips, `test_framework.py::test_cu_btc_equals_the_fortran_after_every_run` skips
  only its archive comparison (with a warning), and `test_reference.py` keeps
  reporting the EPSR-test `.o01` as a warning. With a correctly rounded `sinhf`, that
  `.o01` hash is now required to match.
- Docs: the `powf` explanation in `python/docs/testing.md`, `test_reference.py` and
  `test_framework.py` is replaced by the `sinhf` one (new section "libm and the
  Fortran references"). The top-level `README.md` paragraph on the discrepancy is
  replaced.
- Not changed: the engine, `rust/epsr-core/tests/test_case.rs` (it still fails on
  `.o01` with glibc ≤ 2.40, as documented) and the tolerance policy. The
  `EPSR_FORTRAN_BUILD` comparison for Cu-BTC was not added: with the archive made
  on a correctly rounded `sinhf`, it would add only a slower route to the same check.
