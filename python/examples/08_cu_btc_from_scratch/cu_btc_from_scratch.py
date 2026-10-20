"""Refine the Cu-BTC (HKUST-1) framework against neutron data, with Bragg peaks.

Run from the ``python/`` folder:

    uv run python examples/08_cu_btc_from_scratch/cu_btc_from_scratch.py            # ~4 minutes
    uv run python examples/08_cu_btc_from_scratch/cu_btc_from_scratch.py --quick    # ~30 s

One cubic unit cell of the metal-organic framework (a = 26.3336 Å, 624 atoms) is read
from a PDB file (``CuBTC.pdb`` next to this script) and fitted to the NIMROD data of the empty framework, with Bragg peaks
calculated up to Q = 5 1/Å. Unlike example 07 the configuration is not built from
molecules: ``readpdb`` finds the molecules (24 Cu2 paddlewheels and 32 benzene-
tricarboxylate linkers) by bond detection, and the crystal positions are kept. The
schedule is that of ``validation/cu-btc/steps.txt``: equilibrate for 100 iterations on
the reference potential, then two runs of 100 iterations with the empirical potential
(EReq 10 kJ/mol). Everything goes to ``output/`` next to this script.
"""

import argparse
import os
from pathlib import Path

import epsr
from epsr.samples import NeutronSample
from epsr.system import AtomType, Framework, System

EXAMPLE_FOLDER = Path(__file__).resolve().parent
REPOSITORY = EXAMPLE_FOLDER.parents[2]
DATA = REPOSITORY / "experimental_data" / "cu-btc" / "data"
OUTPUT = Path(os.environ.get("EPSR_EXAMPLES_WORK", EXAMPLE_FOLDER)) / "output"

parser = argparse.ArgumentParser()
parser.add_argument("--quick", action="store_true", help="10 iterations per run")
options = parser.parse_args()
EQUILIBRATION = 11 if options.quick else 100  # includes the 1-iteration scratch start
REFINEMENT = 10 if options.quick else 100

# --- 1. Atom types and the framework ---------------------------------------------------
# Atom types in readpdb's order (first appearance in the PDB); labels = PDB atom names.
# The Lennard-Jones parameters and charges are those of the Dissolve example.
CU = AtomType("CU", element="Cu", epsilon=0.02092, sigma=3.144, charge=1.098)
OX = AtomType("O", element="O", epsilon=0.401664, sigma=3.033, charge=-0.665)
C1 = AtomType("C1", element="C", epsilon=0.39748, sigma=3.473, charge=0.778)
C2 = AtomType("C2", element="C", epsilon=0.39748, sigma=3.473, charge=-0.092)
C3 = AtomType("C3", element="C", epsilon=0.39748, sigma=3.473, charge=-0.014)
H = AtomType("H", element="H", epsilon=0.06276, sigma=2.846, charge=0.109, mass=1.008)

framework = Framework.from_pdb(  # readpdb; checked against the PDB label order
    EXAMPLE_FOLDER / "CuBTC.pdb",  # from validation/cu-btc/ (made by make_pdb.py)
    types=[CU, OX, C1, C2, C3, H],
    vdw_factor=0.68,  # bond detection: finds the Dissolve bonds + the Cu-Cu pairs
    dihedrals=False,
)
print(framework)  # cubic cell 26.3336 Å, 624 atoms, 56 residues (24 × 2, 32 × 18 atoms)
# The 192 Cu-O bonds of the crystal join paddlewheels to linkers, so they are
# intermolecular here: Lennard-Jones + Coulomb only, the main approximation of the model.

system = System(
    "CuBTC", temperature=300.0, structure=framework
)  # density from the cell
system.monte_carlo(
    tol=0.0, stepmi=0.1, stepri=0.3, stepr=0.3, stepm=1.0
)  # changeato names
print(system.pair_parameters())

# --- 2. The data ----------------------------------------------------------------------
samples = [NeutronSample("CuBTC", DATA / "Empty-CuBTC.mint01")]  # natural, no exchange

# --- 3. Build the EPSR project --------------------------------------------------------
project = (
    epsr.Project.create(  # readpdb -> changeato -> nwts -> .EPSR.inp; no EPSR run yet
        OUTPUT / "epsr_run",
        system,
        samples,
        settings={
            "qwidthqmax": (0.2, 5.0),  # Bragg peaks up to 5 1/Å, Gaussian width 0.2 1/Å
            "hklqmin": (
                0.0,
                10.0,
                0.0,
            ),  # minimum Q of hkl, radius for the Bragg g(r), DW
            "qmin": (0.05, 0.0),
            "rshmin": 0.7,  # below the 0.93 Å C-H bond
        },
        overwrite=True,
    )
)
print(project.status())

# --- 4. Refine ------------------------------------------------------------------------
equilibration = project.equilibrate(
    iterations=EQUILIBRATION
)  # scratch + 99, Bragg peaks
history = project.refine(ereq=10.0, runs=2, iterations=REFINEMENT)

# --- 5. Results -----------------------------------------------------------------------
result = project.results()
print()
print(f"R-factor: {equilibration.r_factor[0]:.3e} after the scratch start")
print(f"          {equilibration.r_factor[-1]:.3e} on the reference potential")
print(f"          {history.r_factor[-1]:.3e} with the empirical potential")
bragg = epsr.io.read_bragg(project.output("bragg"))
print(f"{bragg.n_reflections} Bragg reflections at {bragg.n_unique} unique Q values")
copper_copper = epsr.analysis.first_peak(
    result.r, result.gr_intra("CU", "CU"), (2.0, 3.2)
)
copper_oxygen = epsr.analysis.first_peak(
    result.r, result.gr_inter("CU", "O"), (1.5, 2.8)
)
print(f"Cu-Cu (paddlewheel) at {copper_copper:.2f} Å")
print(
    f"Cu-O (intermolecular in this model) at {copper_oxygen:.2f} Å (1.95 Å in the crystal)"
)
print(result.potential.pair_params())  # CU-O row: LJ + Coulomb (the Cu-O "bond")
