# Experimental data

Measured neutron total-scattering data, one folder per system. Each folder has the same
layout:

| subfolder | contents |
|---|---|
| `data/` | reduced data from Gudrun (`.mint01`: Q, F(Q), error), the files EPSR fits |
| `dissolve/` | the [Dissolve](https://www.projectdissolve.com/) example inputs the system was taken from (`input.txt`, `CMakeLists.txt`). These are Dissolve files, not EPSR inputs. |
| `reference_partials/` | published partial g(r), where there are any (see the `SOURCE.md` inside) |

| system | data (ISIS instrument) | EPSR setup | used by |
|---|---|---|---|
| [water/](water) | SANDALS: H₂O, D₂O, HDO 50:50 (`SLS18498/18500/18502.mint01`); Soper (2000) reference partials | [`validation/water`](../validation/water) | `python/examples/05`, `07`; `python/tests/test_experimental_data.py` |
| [silica_bulk/](silica_bulk) | NIMROD: amorphous SiO₂ (`SiO2_NIMROD.mint01`) | [`validation/silica_bulk`](../validation/silica_bulk) | `python/examples/06`, `09`; `python/tests/test_experimental_data.py` |
| [cu-btc/](cu-btc) | NIMROD: empty Cu-BTC (HKUST-1) framework (`Empty-CuBTC.mint01`); the crystal structure `dissolve/7108574.cif` | [`validation/cu-btc`](../validation/cu-btc) | `python/examples/08`; `python/tests/test_cu_btc.py` |
| [benzene/](benzene) | NIMROD: liquid C₆H₆, C₆D₆ and a 50:50 mixture (runs 27749, 14901, 27731) | none (set up in Python only) | `python/examples/10`; first gathered to test the archived Mellor agent on data it had not seen |

The EPSR inputs built from these systems (molecules, weights files, `.EPSR.inp`,
`steps.txt`) are in [`../validation/`](../validation/README.md). Each validation case
folder holds its own copy of the `.mint01` files, because the validation scripts copy
the whole case folder.

The original EPSR test case (200 DMSO + 400 water, three neutron data sets) is in
[`../EPSR-test/`](../EPSR-test), and the feature tests built from it are in
[`../test_problems/`](../test_problems/README.md).
