"""The validated systems written with the setup API (``epsr.System``, ``NeutronSample``).

Every number is that of ``validation/water`` and ``validation/silica_bulk``; a
project built from these definitions must give the validated files.
"""

from __future__ import annotations

from epsr import AtomType, Molecule, NeutronSample, System

from epsr_cases import REPOSITORY

WATER_DATA = REPOSITORY / "experimental_data" / "water" / "data"
SILICA_DATA = REPOSITORY / "experimental_data" / "silica_bulk" / "data"
#: The committed example 07 run, identical to validation/water (plan constraint 1).
WATER_RUN = (
    REPOSITORY / "python" / "examples" / "07_water_from_scratch" / "output" / "epsr_run"
)

OW = AtomType("OW", "O", 0.6503, 3.165492, -0.82, coarse_grain={"rmax": 1.55})
HW = AtomType("HW", "H", 0.0, 0.0, 0.41, 2.0, exchangeable=True)
SI = AtomType("Si", "Si", 0.175, 1.03, 2.0, coarse_grain={"rmax": 1.11})
OX = AtomType("O", "O", 0.165, 3.5, -1.0, coarse_grain={"rmax": 1.55})

#: settings of validation/water and validation/silica_bulk that differ from EPSR's defaults
WATER_SETTINGS = {"rshmin": 0.7}
SILICA_SETTINGS = {"qmin": (0.1, 0.0), "rshmin": 1.2}


def water_molecule() -> Molecule:
    return Molecule.from_geometry(
        "Water",
        centre=("O", OW),
        arms=[("H1", HW, 0.976), ("H2", HW, 0.976)],
        angle=107.134,
    )


def water_system() -> System:
    return System("waterbox", temperature=300.0, density=0.1).add(
        water_molecule(), 1000
    )


def water_samples() -> list[NeutronSample]:
    return [
        NeutronSample(
            "D2O", WATER_DATA / "SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}
        ),
        NeutronSample(
            "HDO",
            WATER_DATA / "SLS18500-HDO5050.mint01",
            isotopes={"HW": {1: 0.5, 2: 0.5}},
        ),
        NeutronSample("H2O", WATER_DATA / "SLS18498-H2O.mint01"),
    ]


def silica_system() -> System:
    system = System("silica", temperature=298.0, density=0.0664)
    return system.add(Molecule.atom(SI), 1500).add(Molecule.atom(OX), 3000)


def silica_samples() -> list[NeutronSample]:
    return [NeutronSample("SiO2", SILICA_DATA / "SiO2_NIMROD.mint01")]
