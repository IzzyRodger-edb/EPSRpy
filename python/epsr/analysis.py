"""Common post-processing of partial g(r): peak positions and coordination numbers.

Plain numpy, so these work without the ``plot`` extra. They take the arrays of a
:class:`epsr.Simulation` (``result.r``, ``result.gr_inter("OW", "HW")``):

>>> peak = first_peak(result.r, result.gr_inter("OW", "OW"), (2.0, 4.0))
>>> n = coordination_number(result.r, result.gr_inter("OW", "OW"),
...                         partner_density=result.density / 3, cutoff=3.3)
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["first_peak", "coordination_number"]


def first_peak(
    r: npt.ArrayLike, g: npt.ArrayLike, window: tuple[float, float]
) -> float:
    """Position (Å) of the highest point of ``g`` with ``window[0] < r < window[1]``.

    The window keeps the search on the peak of interest (e.g. ``(1.5, 2.3)`` for the
    water hydrogen bond, which a search over all r would miss).
    """
    r_values, g_values = np.asarray(r, dtype=float), np.asarray(g, dtype=float)
    if r_values.shape != g_values.shape:
        raise ValueError(
            f"r and g differ in shape: {r_values.shape} vs {g_values.shape}"
        )
    low, high = window
    inside = (r_values > low) & (r_values < high)
    if not inside.any():
        raise ValueError(f"no r values inside the window {window}")
    return float(r_values[inside][np.argmax(g_values[inside])])


def coordination_number(
    r: npt.ArrayLike,
    g: npt.ArrayLike,
    partner_density: float,
    cutoff: float,
    r_min: float = 0.0,
) -> float:
    """Average number of partner atoms between ``r_min`` and ``cutoff`` (Å).

    N = 4π ρ ∫ r² g(r) dr, with ``partner_density`` ρ the number density of the
    partner atoms (atoms/Å³; e.g. ``result.density`` times the atom fraction of the
    partner type). The integral uses the trapezoidal rule, because the EPSR r grid
    is not uniform (its first interval is a quarter of the others).
    """
    r_values, g_values = np.asarray(r, dtype=float), np.asarray(g, dtype=float)
    if r_values.shape != g_values.shape:
        raise ValueError(
            f"r and g differ in shape: {r_values.shape} vs {g_values.shape}"
        )
    inside = (r_values >= r_min) & (r_values <= cutoff)
    r_inside = r_values[inside]
    integrand = r_inside**2 * g_values[inside]
    integral = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r_inside))
    return float(4 * np.pi * partner_density * integral)
