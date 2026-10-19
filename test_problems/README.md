# Test problems

Verification cases for the Rust port. Every case comes with an archived Fortran
reference, and the Rust output must match it byte for byte.

| path | what |
|---|---|
| [`../EPSR-test/`](../EPSR-test) | the original EPSR test case supplied with the Fortran: 200 DMSO + 400 water, fitted to three neutron data sets. It is the base of every case below, and `rust/epsr-core/tests` pins FNV hashes of its outputs. It stays at the top level because the case scripts and tests reach it by relative path. |
| [`200DMSO400Water/EPSR-test-ex*/`](200DMSO400Water) | 115 feature cases, one per engine feature or tool (X-ray data, Gaussian/Morse potentials, triclinic boxes, every setup and analysis program, error paths, ...). |

Each `EPSR-test-exN/` folder holds the inputs, `commands.txt` (one program call per
line; `repeat N` chains the list N times), a `README.md` describing what is tested, a
one-page summary `EPSR-test-exN.pdf`, and `fortran_reference.tar.gz` with the files the
Fortran wrote.

- Index of all cases: [`rust/docs/TEST_CASES.md`](../rust/docs/TEST_CASES.md), generated
  by `rust/tools/gen_test_index.py`
- Run them against the Fortran: `rust/tools/run_ex_test.sh` or
  `rust/tools/run_all_ex_tests.sh` (see [`rust/docs/VALIDATION.md`](../rust/docs/VALIDATION.md))
- Replay them through the Python wrapper: `python/tests/test_ex_cases.py`

Validation on real measured data (water, silica, Cu-BTC) is in
[`../validation/`](../validation/README.md), and the data themselves are in
[`../experimental_data/`](../experimental_data/README.md).
