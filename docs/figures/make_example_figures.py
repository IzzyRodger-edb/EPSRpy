"""Figures of the three example systems as EPSR starts them: configuration, reference
potential and measured data.

Run from the ``python/`` folder (needs the ``plot`` extra):

    uv run python ../docs/figures/make_example_figures.py

For each case (water, amorphous silica, Cu-BTC) the system is built with the setup API
exactly as in its example (``Project.create``; no EPSR run), and the folder is read with
``project.results()``, which changes no file:

* **left** - the initial configuration the Monte Carlo starts from: the atoms of a slab
  through the middle of the box, seen along z;
* **middle** - the reference potential U(r) the engine uses for selected atom-type pairs
  (Lennard-Jones + Coulomb, ``sim.potential.reference()``);
* **right** - the measured interference differential cross section F(Q) being fitted
  (the ``.mint01`` data, with their errors as a band).

Writes ``water.png``, ``silica.png`` and ``cu_btc.png`` next to this script.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import epsr  # noqa: E402
from epsr import AtomType, Framework, Molecule, NeutronSample, System  # noqa: E402

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
DATA = REPOSITORY / "experimental_data"

# Reference palette (dataviz skill): three categorical slots validated all-pairs in
# light mode, one muted ink for a fourth group, chart chrome.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")
MUTED = "#898781"
INK, INK_2 = "#0b0b0b", "#52514e"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
SLAB = 8.0  # Å, thickness of the configuration slab


@dataclass
class Case:
    name: str
    title: str
    system: System
    samples: list[NeutronSample]
    settings: dict[str, object]
    #: element -> colour of the configuration panel (in legend order)
    colours: dict[str, str]
    #: atom-type pairs drawn in the potential panel, and the y range (kJ/mol)
    pairs: list[tuple[str, str]]
    potential_range: tuple[float, float]


def water() -> Case:
    ow = AtomType("OW", "O", 0.6503, 3.165492, -0.82, coarse_grain={"rmax": 1.55})
    hw = AtomType("HW", "H", 0.0, 0.0, 0.41, 2.0, exchangeable=True)
    molecule = Molecule.from_geometry(
        "Water", centre=("O", ow), arms=[("H1", hw, 0.976), ("H2", hw, 0.976)],
        angle=107.134,
    )
    data = DATA / "water" / "data"
    return Case(
        "water",
        "Water: 1000 H2O, 0.1 atoms/Å³, 300 K (SANDALS)",
        System("waterbox", 300.0, 0.1).add(molecule, 1000),
        [
            NeutronSample("D2O", data / "SLS18502-D2O.mint01", isotopes={"HW": {2: 1.0}}),
            NeutronSample(
                "HDO", data / "SLS18500-HDO5050.mint01", isotopes={"HW": {1: 0.5, 2: 0.5}}
            ),
            NeutronSample("H2O", data / "SLS18498-H2O.mint01"),
        ],
        {"rshmin": 0.7},
        {"O": SERIES[0], "H": SERIES[1]},
        [("OW", "OW"), ("OW", "HW"), ("HW", "HW")],
        (-400.0, 600.0),
    )


def silica() -> Case:
    si = AtomType("Si", "Si", 0.175, 1.03, 2.0, coarse_grain={"rmax": 1.11})
    oxygen = AtomType("O", "O", 0.165, 3.5, -1.0, coarse_grain={"rmax": 1.55})
    system = System("silica", 298.0, 0.0664)
    system.add(Molecule.atom(si), 1500).add(Molecule.atom(oxygen), 3000)
    return Case(
        "silica",
        "Amorphous silica: 1500 Si + 3000 O, 0.0664 atoms/Å³, 298 K (NIMROD)",
        system,
        [NeutronSample("SiO2", DATA / "silica_bulk" / "data" / "SiO2_NIMROD.mint01")],
        {"qmin": (0.1, 0.0), "rshmin": 1.2},
        {"Si": SERIES[0], "O": SERIES[1]},
        [("Si", "Si"), ("Si", "O"), ("O", "O")],
        (-2500.0, 2500.0),
    )


def cu_btc() -> Case:
    types = [
        AtomType("CU", "Cu", 0.02092, 3.144, 1.098),
        AtomType("O", "O", 0.401664, 3.033, -0.665),
        AtomType("C1", "C", 0.39748, 3.473, 0.778),
        AtomType("C2", "C", 0.39748, 3.473, -0.092),
        AtomType("C3", "C", 0.39748, 3.473, -0.014),
        AtomType("H", "H", 0.06276, 2.846, 0.109, 1.008),
    ]
    framework = Framework.from_pdb(REPOSITORY / "validation" / "cu-btc" / "CuBTC.pdb", types)
    system = System("CuBTC", 300.0, structure=framework)
    system.monte_carlo(tol=0.0, stepmi=0.1, stepri=0.3, stepr=0.3, stepm=1.0)
    return Case(
        "cu_btc",
        "Cu-BTC (HKUST-1): one 624-atom unit cell, 300 K (NIMROD)",
        system,
        [NeutronSample("CuBTC", DATA / "cu-btc" / "data" / "Empty-CuBTC.mint01")],
        {"qwidthqmax": (0.2, 5.0), "hklqmin": (0.0, 10.0, 0.0), "qmin": (0.05, 0.0),
         "rshmin": 0.7},
        {"Cu": SERIES[0], "O": SERIES[1], "C": SERIES[2], "H": MUTED},
        [("CU", "O"), ("CU", "CU"), ("O", "H")],
        (-400.0, 400.0),
    )


def _style(axes: plt.Axes) -> None:
    axes.set_facecolor(SURFACE)
    axes.grid(True, color=GRID, linewidth=0.6)
    axes.tick_params(colors=INK_2, labelsize=8)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(AXIS)
    axes.xaxis.label.set_color(INK_2)
    axes.yaxis.label.set_color(INK_2)


def _configuration(axes: plt.Axes, result: epsr.Simulation, case: Case) -> None:
    positions = result.positions
    side = result.box_size
    elements = {t.label: t.element for t in case.system.atom_types}
    element_of_atom = np.array([elements.get(label, "?") for label in result.atom_labels])
    centre = 0.5 * side if positions[:, 2].min() >= 0 else 0.0
    in_slab = np.abs(positions[:, 2] - centre) < SLAB / 2
    for element, colour in case.colours.items():
        selected = in_slab & (element_of_atom == element)
        size = 6 if element == "H" else 14
        axes.scatter(positions[selected, 0], positions[selected, 1], s=size, c=colour,
                     edgecolors=SURFACE, linewidths=0.3, label=element)
    low = positions[:, 0].min() if positions[:, 0].min() < 0 else 0.0
    axes.add_patch(plt.Rectangle((low, low), side, side, fill=False, edgecolor=AXIS))
    axes.set_aspect("equal")
    axes.set_xlabel("x (Å)")
    axes.set_ylabel("y (Å)")
    axes.set_title(f"Initial configuration ({SLAB:g} Å slab, box {side:.2f} Å)",
                   fontsize=9, color=INK, loc="left")
    axes.legend(fontsize=7, frameon=False, loc="upper right", markerscale=1.5,
                labelcolor=INK_2, bbox_to_anchor=(1.0, 1.12), ncol=len(case.colours))


def _potential(axes: plt.Axes, result: epsr.Simulation, case: Case) -> None:
    r = result.r
    for (first, second), colour in zip(case.pairs, SERIES, strict=False):
        axes.plot(r, result.potential.reference(first, second), color=colour,
                  linewidth=1.5, label=f"{first}-{second}")
    axes.axhline(0.0, color=AXIS, linewidth=0.8)
    axes.set_xlim(0.0, 12.5)
    axes.set_ylim(*case.potential_range)
    axes.set_xlabel("r (Å)")
    axes.set_ylabel("U(r) (kJ/mol)")
    axes.set_title("Reference potential (Lennard-Jones + Coulomb)", fontsize=9,
                   color=INK, loc="left")
    axes.legend(fontsize=7, frameon=False, labelcolor=INK_2)


def _data(axes: plt.Axes, case: Case) -> None:
    offset_step = 0.5 if len(case.samples) > 1 else 0.0
    for index, (sample, colour) in enumerate(zip(case.samples, SERIES, strict=False)):
        data = sample.data
        assert data is not None
        valid = data.valid
        q, f, errors = data.q[valid], data.f[valid] + index * offset_step, data.errors[valid]
        axes.fill_between(q, f - errors, f + errors, color=colour, alpha=0.25, linewidth=0)
        label = sample.name + (f" (+{index * offset_step:g})" if index and offset_step else "")
        axes.plot(q, f, color=colour, linewidth=1.0, label=label)
    axes.set_xscale("log")
    axes.set_xlabel("Q (Å⁻¹)")
    axes.set_ylabel("F(Q) (barn sr⁻¹ atom⁻¹)")
    axes.set_title("Measured data being fitted", fontsize=9, color=INK, loc="left")
    axes.legend(fontsize=7, frameon=False, labelcolor=INK_2)


def draw(case: Case, work: Path) -> Path:
    project = epsr.Project.create(work / case.name, case.system, case.samples, case.settings)
    result = project.results()
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True,
                                facecolor=SURFACE)
    for panel in axes:
        _style(panel)
    _configuration(axes[0], result, case)
    _potential(axes[1], result, case)
    _data(axes[2], case)
    figure.suptitle(case.title, fontsize=11, color=INK, x=0.01, ha="left")
    path = HERE / f"{case.name}.png"
    figure.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(figure)
    return path


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as work:
        for make in (water, silica, cu_btc):
            print(draw(make(), Path(work)))
