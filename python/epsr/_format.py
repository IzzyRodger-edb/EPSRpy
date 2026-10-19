"""Numbers and keyword lines in the layouts the EPSR programs write.

The generated input files use the layout the programs themselves write when they
rewrite a file (Fortran edit descriptors, reals rounded to single precision), so a
generated file is the file EPSR would produce from the same numbers.
"""

from __future__ import annotations

import math
import struct

__all__ = [
    "f32",
    "fixed",
    "exponent",
    "list_real",
    "list_real_field",
    "number_text",
    "keyword_line",
    "KEYWORD_COLUMN_WIDTH",
    "DESCRIPTION_GAP",
    "TITLE_MARKER",
]

#: Keyword files: keyword padded to 11, a blank, the value, 15 blanks, the description
#: (Fortran format ``(a11,1x,a,15x,a)``).
KEYWORD_COLUMN_WIDTH = 11
DESCRIPTION_GAP = " " * 15
TITLE_MARKER = "Title of this file"


def f32(value: float) -> float:
    """``value`` rounded to single precision, as EPSR stores every real."""
    rounded: float = struct.unpack("f", struct.pack("f", float(value)))[0]
    return rounded


def fixed(value: float, width: int, decimals: int) -> str:
    """Fortran ``Fw.d`` of a single-precision value."""
    return f"{f32(value):{width}.{decimals}f}"


def exponent(value: float, width: int = 13, decimals: int = 6) -> str:
    """Fortran ``Ew.d``: ``0.100000E+01`` style (mantissa below 1), right-justified."""
    value = f32(value)
    if value == 0.0:
        mantissa, power = 0.0, 0
    else:
        power = math.floor(math.log10(abs(value))) + 1
        mantissa = value / 10.0**power
        if abs(round(mantissa, decimals)) >= 1.0:  # rounding carried into a new digit
            mantissa /= 10.0
            power += 1
    text = f"{mantissa:.{decimals}f}E{power:+03d}"
    return f"{text:>{width}}"


def list_real(value: float) -> str:
    """A single-precision real as Fortran list-directed output prints it.

    ``0.976000011``, ``15.9994001``, ``4.99999989E-02``.
    """
    value = f32(value)
    if value == 0.0:
        return "0.00000000"
    power = math.floor(math.log10(abs(value)))
    if -1 <= power <= 8:
        return f"{value:.{8 - power}f}"
    return f"{value:.8E}"


def list_real_field(value: float) -> str:
    """A single-precision real as one item of a Fortran list-directed record, without
    the separating blank: a fixed-form value right-justified in 12 columns plus four
    trailing blanks (where the exponent would be), an E-form value in 16 columns."""
    text = list_real(value)
    return f"{text:>16}" if "E" in text else f"{text:>12}    "


def number_text(value: float) -> str:
    """A number as EPSR input text: ``1`` for whole numbers, else the shortest exact form."""
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"cannot write {value!r} to an EPSR file")
    return str(int(value)) if value.is_integer() else repr(value)


def keyword_line(keyword: str, value: str, description: str) -> str:
    """One line of a keyword file: ``keyword value               description``."""
    return f"{keyword:<{KEYWORD_COLUMN_WIDTH}} {value}{DESCRIPTION_GAP}{description}"
