#!/usr/bin/env bash
# Create the EPSR_python_translation repository next to EPSRpy, copy in the files the
# translation needs (see python_translation.md section 9), make the first commit and create
# the GitHub repository. It does not push.
#
# usage: bash docs/plans/setup_epsr_python_translation.sh
#   SRC=...  EPSRpy checkout (default ~/EPSRpy); files come from its committed HEAD
#   DST=...  new repository (default ~/EPSR_python_translation); must not exist
#   GH_OWNER=...         GitHub account (default IzzyRodger-edb)
#   INCLUDE_EX_CASES=0   leave out test_problems/ (185 MB)
set -euo pipefail
SRC=${SRC:-$HOME/EPSRpy}
DST=${DST:-$HOME/EPSR_python_translation}
GH_OWNER=${GH_OWNER:-IzzyRodger-edb}
INCLUDE_EX_CASES=${INCLUDE_EX_CASES:-1}
[ -e "$DST" ] && { echo "$DST already exists"; exit 1; }
mkdir -p "$DST"
git init -q -b main "$DST"
cd "$DST"

# --- 1. copy from EPSRpy's committed HEAD
copy() { git -C "$SRC" archive HEAD "$@" | tar -x -C "$DST"; }
copy fortran_src EPSR-test LICENSE FORTRAN_BUGS.md
copy experimental_data/README.md experimental_data/water experimental_data/silica_bulk experimental_data/cu-btc
copy validation/water validation/silica_bulk validation/cu-btc \
     validation/run_validation.sh validation/set_keyword.py validation/plot_gr_agreement.py validation/README.md
[ "$INCLUDE_EX_CASES" = 1 ] && copy test_problems
copy rust/tools/compare_outputs.py rust/tools/build_fortran_reference.sh rust/tools/build_patched_readxyz_readpdb.sh \
     rust/tools/fortran_link_lists.txt rust/tools/run_ex_test.sh rust/tools/gen/gen_symtabs.py \
     rust/docs/ARCHITECTURE.md rust/docs/FORTRAN_MAP.md rust/docs/TOOLS.md rust/docs/TEST_CASES.md rust/docs/VALIDATION.md

# --- 2. rearrange
mkdir -p tools docs/rust_port
mv rust/tools/* tools/
mv rust/docs/* docs/rust_port/
rm -rf rust
mv FORTRAN_BUGS.md docs/
mv validation/README.md validation/RUST_PORT_RESULTS.md
rm -f fortran_src/.DS_Store
cp "$SRC/docs/plans/python_translation.md" PLAN.md
sed -i 's#(setup_epsr_python_translation.sh)#(not copied; it lives in EPSRpy docs/plans/)#' PLAN.md

# --- 3. fix paths that pointed into EPSRpy's rust/ folder
sed -i 's#\$HERE/\.\./\.\./fortran_src#$HERE/../fortran_src#' tools/build_fortran_reference.sh tools/build_patched_readxyz_readpdb.sh
sed -i 's#\.\./\.\./\.\./fortran_src#../../fortran_src#' tools/gen/gen_symtabs.py
sed -i 's#\.\./rust/tools/compare_outputs\.py#../tools/compare_outputs.py#' validation/run_validation.sh
grep -q '/../fortran_src' tools/build_fortran_reference.sh
grep -q '../tools/compare_outputs.py' validation/run_validation.sh

# --- 4. list the Fortran files in scope
{
  echo "# Fortran files covered by the Python translation (generated from tools/fortran_link_lists.txt"
  echo "# plus the shell, setup, plot and Jmol programs). Every other file in fortran_src/ is an"
  echo "# older revision or an out-of-scope program, kept for reference only."
  for b in $( (grep -v '^#' tools/fortran_link_lists.txt | cut -d: -f2 | tr ' ' '\n' | grep '\.o$' | sed 's/\.o$//'
              echo molecule_coordinates_mod) | sort -u); do
    if [ -f "fortran_src/$b.f90" ]; then echo "$b.f90"; else echo "$b.f"; fi
  done
  printf '%s\n' epsrshell_standalone.f setup_generic.f plot.f runjmol.f plot3djmol.f
  (cd fortran_src && ls *.inc)
} > fortran_src/ACTIVE_SOURCES.txt

# --- 5. package skeleton
for d in fortran common library programs shell; do
  mkdir -p "src/epsr26/$d"; : > "src/epsr26/$d/__init__.py"
done
mkdir -p tests
cat > src/epsr26/__init__.py <<'EOF'
"""Line-by-line Python translation of EPSR26 (A. K. Soper). See PLAN.md."""
__version__ = "0.0.1"
EOF
cat > src/epsr26/__main__.py <<'EOF'
from epsr26.shell.epsrshell import main

main()
EOF
cat > src/epsr26/shell/epsrshell.py <<'EOF'
"""epsrshell_standalone.f: the EPSR command shell. Not translated yet (PLAN.md phase 3)."""


def main():
    raise SystemExit("epsrshell: not translated yet")
EOF
cat > pyproject.toml <<'EOF'
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "epsr26"
version = "0.0.1"
description = "Line-by-line Python translation of the EPSR26 Fortran suite"
readme = "README.md"
license = "MIT"
requires-python = ">=3.10"
dependencies = ["numpy>=2"]

[project.optional-dependencies]
dev = ["pytest", "pexpect", "psutil", "ruff"]

[project.scripts]
epsrshell = "epsr26.shell.epsrshell:main"

[tool.hatch.build.targets.wheel]
packages = ["src/epsr26"]

[tool.pytest.ini_options]
testpaths = ["tests"]
markers = ["slow: full validation runs (chain mode)"]
EOF
cat > .gitignore <<'EOF'
__pycache__/
*.egg-info/
.venv/
build/
dist/
references/
work/
fortran_build/
EOF
cat > CLAUDE.md <<'EOF'
Never run git push. I handle all git pushes manually.

This repository translates the EPSR26 Fortran (fortran_src/) into Python line by line.
Follow the translation rules in PLAN.md section 2. Fidelity to the Fortran comes before
idiomatic Python, and output must stay byte-identical to the Fortran references.
EOF
cat > README.md <<'EOF'
# EPSR_python_translation

A line-by-line Python translation of EPSR26 (Empirical Potential Structure Refinement,
A. K. Soper), written to be read and modified rather than to be fast. It keeps the
`epsrshell` interface and is validated byte for byte against the Fortran on measured
neutron data for water, amorphous silica and Cu-BTC.

Status: setting up. See [PLAN.md](PLAN.md).

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e '.[dev]'
```
EOF

# --- 6. first commit and GitHub repository (private; no push)
git add -A
git commit -q -m "Initial import: Fortran source, validation cases and tools from EPSRpy; plan and package skeleton"
gh repo create "$GH_OWNER/EPSR_python_translation" --private --source=. --remote=origin \
  --description "Line-by-line Python translation of the EPSR26 Fortran suite"
echo
echo "Created $DST. Push when ready:  cd $DST && git push -u origin main"
