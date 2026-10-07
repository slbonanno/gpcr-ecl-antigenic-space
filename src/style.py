"""Shared plotting style. Imported by every figure notebook.

Arial everywhere, with fallbacks so it still renders on a machine that does not
have it installed. Nothing here draws anything.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

SANS = ["Arial", "Helvetica", "Helvetica Neue", "Liberation Sans",
        "Nimbus Sans", "DejaVu Sans"]

PALETTE = {
    "N-term": "#8172B3",
    "ECL1": "#4C72B0",
    "ECL2": "#DD8452",
    "ECL3": "#55A868",
    "ECF": "#937860",
}

FORMAT_PALETTE = {
    "linear peptide": "#4C72B0",
    "cyclic / disulfide-constrained": "#DD8452",
    "folded ectodomain": "#8172B3",
    "cell-surface only": "#C44E52",
}


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": SANS,
        "mathtext.fontset": "dejavusans",
        "pdf.fonttype": 42,          # editable text in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.dpi": 130,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 9,
        "axes.titlesize": 10,
        "legend.frameon": False,
    })
    return plt


def available() -> bool:
    """True if a real Arial-like face resolved, False if matplotlib fell back."""
    from matplotlib.font_manager import findfont, FontProperties
    # family must be a list: a bare string with a hyphen is parsed as a
    # fontconfig pattern and blows up on "sans-serif"
    f = findfont(FontProperties(family=SANS), fallback_to_default=True)
    return "DejaVu" not in f
