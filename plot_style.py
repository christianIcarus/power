"""House plot style: seaborn whitegrid at cmr10.

Mirrors _style() and the palette from Icarus-Matrix/vehicle-simulation's
plots.py so solar figures and endurance figures are visually the same family.
Kept as a copy rather than an import because the two repos are not installed
into one another; if the endurance house style changes, change it here too.

Three things cmr10 (Computer Modern Roman) needs and does not announce:

  * It is a SINGLE-WEIGHT face. Anything asking for bold silently falls back
    to regular and matplotlib logs `findfont: Failed to find font weight bold`
    once per glyph, which buries real warnings. The fallback is what we want,
    so the logger is quieted rather than the styling changed.
  * It has no U+2212 MINUS SIGN, so `axes.unicode_minus` must be False or
    every negative tick renders as a missing-glyph box.
  * It has no en/em dash either -- figure text uses `--`. Docstrings are free
    to use whatever, since they are never drawn.
"""
from __future__ import annotations

import logging

# Ink / structure, then the semantic series colours.
INK, MUTED, GRID = "#1a1a1f", "#6b6b76", "#dcdce4"
SOLAR, DEMAND, BATT, ACCENT, WARN = (
    "#e0a020", "#c0392b", "#2e8b57", "#6a4ca8", "#b03030")

_applied = False


def apply(force: bool = False) -> None:
    """Install the house style into matplotlib's global rcParams.

    Idempotent: figure-drawing code paths call this freely without stacking
    seaborn theme applications.
    """
    global _applied
    if _applied and not force:
        return

    import matplotlib.pyplot as plt
    import seaborn as sns

    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams.update(
        {
            "font.family": "cmr10",
            "mathtext.fontset": "cm",
            "axes.formatter.use_mathtext": True,
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfe",
            "axes.edgecolor": "#2a2a33",
            "axes.linewidth": 1.1,
            "axes.titlesize": 15.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 13.5,
            "font.size": 13,
            "xtick.labelsize": 11.5,
            "ytick.labelsize": 11.5,
            "grid.color": "#dadae3",
            "grid.linewidth": 0.7,
            "legend.frameon": True,
            "legend.framealpha": 0.92,
            "legend.edgecolor": "#c8c8d0",
            "legend.fontsize": 11.5,
            "figure.dpi": 130,
            "savefig.dpi": 175,
            "savefig.bbox": "tight",
        }
    )
    _applied = True


def heat_cmap():
    """The house sequential colormap (deep blue -> teal -> pale yellow)."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "heat", ["#12263f", "#2f6fb0", "#63c9a8", "#f2e58a"])


def label(text: str) -> str:
    """Make text safe to DRAW in cmr10.

    cmr10 has no underscore glyph -- matplotlib silently substitutes an
    apostrophe or a raised dot, so `America/Los_Angeles` renders as
    `America/Los'Angeles`. Identifiers that reach a figure (timezone names,
    log-file stems, column names) go through here first.
    """
    return text.replace("_", " ")
