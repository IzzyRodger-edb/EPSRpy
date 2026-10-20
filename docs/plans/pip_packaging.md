# Plan: distributing `epsr` with pip (status and hand-over)

**Status (2026-10-20):** a pip install from source works and is verified. The CI wheel
workflow is written but **has never run**. Nothing is published on PyPI. Local `main`
holds all of this as merge commit `266786b`, which was **not pushed** when this was
written: check with `git log origin/main..main`.
**For:** an agent picking this up later. Read the *User preferences* first, then *Open items*.

**Goal (from the user):** make the Python wrapper (`python/`, imported as `epsr`) easy to
install and use with `pip install`, and document the setup in `README.md`.

## User preferences and constraints

- **Never run `git push`** (`CLAUDE.md`). The user pushes. Local commits and merges are
  fine. For the same reason, don't merge pull requests or push tags; tell the user the
  exact command instead.
- **Don't publish to PyPI or TestPyPI**, and don't trigger releases. These are
  outward-facing and the user does them.
- The user is new to PyPI and GitHub Actions and asked for plain-language,
  step-by-step explanations (see *PyPI setup* below, written the way they understood it).
- When asked how to share the package with a few people, the agent recommended the
  GitHub-URL install (*Ways to give someone EPSRpy*) and leaving PyPI until there are
  outside users. The user hasn't decided yet.

## What was done

### Packaging (`python/pyproject.toml`, maturin + PyO3)

`python/` was already a maturin mixed Rust/Python project. The extension is
`epsr._core`, built for the stable ABI (`abi3-py310`). Changes:

| change | why |
|---|---|
| `dynamic = ["version"]` | the single source of the version is `python/Cargo.toml` (`version = "0.1.0"`); `epsr.__version__` comes from `CARGO_PKG_VERSION` too |
| `license = "MIT"` (SPDX) + `license-files = ["LICENSE"]`; MIT classifier removed; `maturin>=1.9` | PEP 639; maturin only packages licence files inside `python/`, so **`python/LICENSE` is a copy of the root `LICENSE`** (keep them in sync) |
| keywords, classifiers, `[project.urls]` (GitHub) | PyPI metadata |
| `[project.scripts] epsrvis = "epsr.vis.cli:main"` | console command for `epsr.vis` |
| `[tool.maturin] exclude` (`tests/`, `examples/`, `docs/`, `uv.lock`, `.python-version`, sdist only) | sdist went from 23 MB (example 07's committed outputs) to 441 KB |

The engine is the path dependency `epsr-core = { path = "../rust/epsr-core" }` in
`python/Cargo.toml`. This works for every route:
- **sdist:** maturin copies `rust/epsr-core` into it and rewrites the path.
- **Clone or git URL:** pip ≥ 21.3 builds in place, so `../rust` resolves.

The only runtime data file is `epsr/data/vanderWaalsRadii.txt`, loaded with
`importlib.resources`. It is in the wheel and the sdist. The runtime needs nothing else
from the repository; the examples do (`experimental_data/`).

### CI: `.github/workflows/wheels.yml` (new, untested)

- **Wheels** (`PyO3/maturin-action@v1`, `working-directory: python`,
  `--release --strip`): Linux manylinux2014 x86_64 (`ubuntu-latest`) and aarch64
  (`ubuntu-24.04-arm`), macOS x86_64 (`macos-13`) and arm64 (`macos-14`), and Windows
  x64. Each job installs its wheel and imports `epsr`.
- **sdist** job.
- **publish** job: runs only for tags `v*` and uses `environment: pypi` with PyPI
  trusted publishing (`id-token: write`, `pypa/gh-action-pypi-publish@release/v1`).
- **Triggers:** `workflow_dispatch` (build only) and tag push `v*` (build and publish).
- GitHub lists a workflow, and shows its *Run workflow* button, only once the file is on
  the **default branch (`main`)**. The user didn't see it at first for this reason.

### Documentation

- **`README.md` §1 "Quick start" became "1. Installation":**
  - Linux dependencies: a table of wheel vs source route, commands for Debian/Ubuntu/WSL2,
    Fedora/RHEL and Arch, Rust via rustup (distribution packages are too old; Ubuntu
    24.04 has 1.75, and ≥ 1.85 is needed), a newer Python via `uv python install`,
    commands to check each dependency, the PEP 668 "externally-managed-environment"
    note, and the `viewer` extra's system libraries.
  - Install routes (GitHub URL, clone), the extras table, pre-built wheels, the release
    workflow, an "Updating" table, running the examples, and the uv development setup.
  - §2 intro and §7 repository map updated (`.github/workflows/` row).
- **`python/README.md`** "Install": pip first, uv for development.
- **Example 07** (`water_from_scratch.py` docstring and its `README.md`): how to run it
  with a pip-installed package.

### Git state at hand-over

- **`feature/python-setup-api`** has commit `e110cc1` (the user's commit of the work
  above), pushed to `origin`.
- **Local `main`:**
  - `0e99a7c`, the user's own commit with the same tree as `e110cc1`.
  - `266786b`, the agent's merge of `origin/main`. It brought in PR #1's merge commit
    `2a5b289`, which local `main` lacked, so the push fast-forwards. The tree is
    unchanged.
  - **Waiting for the user's `git push origin main`.** The agent couldn't
    `git fetch` (no GitHub credentials in the shell), so `origin/main` was as of the
    user's last push.
- `docs/plans/pip_packaging.md` (this file) is uncommitted on `feature/python-setup-api`.

## Verification done (all in fresh venvs in the session scratchpad)

| test | result |
|---|---|
| `maturin sdist`, then install from the sdist with a **cold** Rust build, Python 3.10 (the minimum) | ✔ about 8 min; `import epsr` gives `0.1.0 26: 2017-03-10`; `epsrvis --help` works |
| `maturin build --release --strip`, install the wheel | ✔ wheel 1.6 MB, `.so` 3.4 MB (unstripped 23 MB); contents: the package, data file, licence and entry point only |
| example 07 `--quick --no-3d` against the installed wheel, run from outside `python/`, `EPSR_EXAMPLES_WORK` pointing at scratch | ✔ finishes, about 75 s; R-factor with the empirical potential 9.8e-4 |
| `pytest -m "not slow" -n 6` against the installed wheel (venv's `pytest`, not `python -m pytest`; the run imported site-packages `epsr`) | ✔ 227 passed |
| `epsr-py @ git+file:///…/EPSRpy@feature/python-setup-api#subdirectory=python` | ✔ builds and installs (tested on commit `89c5fd1`, before the packaging changes; the mechanism is the same) |
| `uv sync --all-extras` (dev setup) with the dynamic version | ✔ `python/uv.lock` updated: `version` line dropped, extras reordered |
| black and ruff on the edited example | ✔ |

**pip's reinstall behaviour at an unchanged version** (tested; it matters for updates):

| installed from | result |
|---|---|
| local folder (`pip install ./python`) | always rebuilt and reinstalled |
| local wheel file | skipped ("already installed with the same version"); needs `--force-reinstall` |
| git URL (also with `-U`) | clones, reads the metadata, **installs nothing**; needs `--force-reinstall` or a version bump |
| PyPI | a version can never be re-uploaded |

So every release must bump `version` in `python/Cargo.toml`.

### Sandbox notes for the next agent

- Plain `pip`'s network access was blocked in the agent's shell: isolated builds hung
  while fetching maturin. `uv pip` had network access. To test pip itself, install
  maturin into the venv, put the venv's `bin/` on `PATH`, and use
  `--no-build-isolation --no-index`.
- `curl` to pypi.org was denied, so **the availability of the name `epsr-py` on PyPI
  wasn't checked**.
- `rm -rf` was denied; use new output folders instead.
- A cold release build of the engine takes about 8 minutes. Point `CARGO_TARGET_DIR` at
  `python/target` to reuse the warm cache.
- A wheel built on this machine is tagged `linux_x86_64`, not manylinux, because the
  host glibc is 2.43. It installs only on systems with glibc ≥ 2.43. Portable wheels
  come from CI (manylinux2014 container). Locally, `maturin build --zig` with
  `pip install ziglang` might work, but it is untested.

## Ways to give someone EPSRpy (as explained to the user)

1. **GitHub URL (works now, no setup).** The other person needs Rust and a C compiler
   (README §1) and waits about 8 minutes for the build:
   `pip install "epsr-py[plot] @ git+https://github.com/IzzyRodger-edb/EPSRpy#subdirectory=python"`.
   Until `main` is pushed, add `@feature/python-setup-api` after `EPSRpy`. If the
   repository is private, they also need access to it.
2. **A wheel file.** It installs in seconds with no Rust, but must be portable. Run the
   workflow by hand (Actions → wheels → Run workflow) and download the wheels from the
   run's *Artifacts*. Nothing is published.
3. **PyPI.** Users run `pip install epsr-py`. It needs the setup below.

## PyPI setup (steps for the user; none done yet)

1. Make a pypi.org account with two-factor authentication. Check that
   pypi.org/project/epsr-py is free; if not, change `name` in `python/pyproject.toml`
   and the docs.
2. On PyPI, go to Account → Publishing → *Add a new pending publisher*:
   - project: `epsr-py`
   - owner: `IzzyRodger-edb`
   - repository: `EPSRpy`
   - workflow: `wheels.yml`
   - environment: `pypi`
3. **Optional:** the GitHub environment `pypi` doesn't need creating (GitHub creates it
   on first use). Create it by hand (Settings → Environments) only to add *Required
   reviewers*, so each publish waits for approval.
4. The user pushes `main`, then does a trial run from the Actions tab (build only).
5. To release:
   1. Bump `version` in `python/Cargo.toml`.
   2. Run `uv sync` (updates the lock files) and commit.
   3. On `main`, run `git tag vX.Y.Z` and `git push origin vX.Y.Z`; the user does this.
   4. Users update with `pip install --upgrade epsr-py`.
6. **Optional:** a dry run on TestPyPI (separate account and publisher; add
   `repository-url: https://test.pypi.org/legacy/` to the publish step).

## Open items

- [ ] **The user pushes `main`** (`git push origin main`). If it is rejected, run
      `git pull origin main` first.
- [ ] **First CI run.** Watch Windows especially: the engine has never been compiled
      there. `setup_common.rs` shells out to `sh`, but fails gracefully. Also check
      the ARM runners.
- [ ] **Offered, not done:** a workflow step that fails the publish when the tag
      (`vX.Y.Z`) doesn't equal the package version. Without it, a mismatched tag
      publishes whatever `Cargo.toml` says, and PyPI versions can't be reused. The
      user hasn't answered.
- [ ] Check the PyPI name `epsr-py` (see above).
- [ ] Source installs ship the unstripped 23 MB `.so`, because the release profile has
      `debug = 1`, kept for profiling. Only CI strips. Setting `strip = true` under
      `[tool.maturin]` would strip every build, including `uv sync` dev builds; it was
      left as it is deliberately.
- [ ] **`epsrvis` name clash:** `visualisation/python_visualisation_toolkit` (package
      `epsrvis`) also installs an `epsrvis` command with the same CLI. They collide
      only if both are installed in one environment.
- [ ] `python/README.md` is the PyPI long description; its relative links (`../rust`,
      `../README.md`) will be broken on PyPI. Use absolute GitHub URLs if that matters.
- [ ] Harmless, pre-existing, and not packaging-related: the setup tools print
      `sh: 1: ./Water.ato: Permission denied`. With no `system_commands.txt`,
      `search_for_files` in `rust/epsr-core/src/tools/setup_common.rs` runs an empty
      `system_ls`, faithfully to the Fortran.
