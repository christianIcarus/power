#!/usr/bin/env python3
"""One typeset image of the full derivation of the spectral mismatch factor M.

Pure exposition -- it recomputes the handful of numbers it quotes from the real
EQE + reference spectra so the sheet can never drift from the code, but it
reads no flight log. Rendered with matplotlib mathtext in the cm fontset so it
matches the cmr10 house style.

Layout note: tall mathtext (\\int, \\frac) is much taller than a text line, so
every element is DRAWN first, its rendered height MEASURED via the renderer,
and the cursor advanced past it. Guessing the advance is what makes equations
collide with prose.

Usage:
    python plot_m_math.py [--out outputs/M_math.png]
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
from pvlib import spectrum

import plot_style
import spectral_mismatch as sm

CUTOFF_NM = 1190.0        # Max7 EQE is < 1% here: the silicon bandgap edge
DAYOFYEAR = 236           # 2026-08-24, the flight this was derived on
HOLD_ALT_M, HOLD_ZENITH = 16765.0, 38.2


def compute_numbers(eqe_csv: pathlib.Path) -> dict:
    wl, sr = sm.load_spectral_response(eqe_csv)
    refs = spectrum.get_reference_spectra()
    w0 = refs.index.to_numpy(float)

    def band(W, E):
        k = (W >= 300) & (W <= 4000)
        W, E = W[k], E[k]
        G = np.trapezoid(E, x=W)
        J = np.trapezoid(E * np.interp(W, wl, sr, left=0.0, right=0.0), x=W)
        unus = np.trapezoid(E[W > CUTOFF_NM], x=W[W > CUTOFF_NM]) / G
        return G, J, J / G, unus

    g15, j15, u15, x15 = band(w0, refs["global"].to_numpy(float))
    w2, d2 = sm.spectrl2_direct(HOLD_ZENITH, HOLD_ALT_M, DAYOFYEAR)
    g55, j55, u55, x55 = band(w2, d2[:, 0])

    # Curves for the graphical derivation, on one common grid so the two
    # spectra and the spectral response can be multiplied pointwise.
    grid = np.arange(300.0, 2500.1, 2.0)
    e15 = np.interp(grid, w0, refs["global"].to_numpy(float))
    e55 = np.interp(grid, w2, d2[:, 0])
    srg = np.interp(grid, wl, sr, left=0.0, right=0.0)
    # How much energy the 2500 nm plot limit hides, so the shaded areas can be
    # described honestly against integrals that actually run to 4000 nm.
    k = (w0 >= 300) & (w0 <= 4000)
    beyond = (np.trapezoid(refs["global"].to_numpy(float)[k & (w0 > 2500)],
                           x=w0[k & (w0 > 2500)])
              / np.trapezoid(refs["global"].to_numpy(float)[k], x=w0[k]))
    return dict(u15=u15, u55=u55, m=u55 / u15, g15=g15, g55=g55,
                j15=j15, j55=j55, x15=x15, x55=x55,
                grid=grid, e15=e15, e55=e55, srg=srg, beyond=beyond,
                eqe_wl=wl, eqe=np.where(wl > 0, sr * 1239.84 / wl, 0.0))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eqe", default="max7_eqe.csv")
    ap.add_argument("--out", default="outputs/M_math.png")
    args = ap.parse_args()

    n = compute_numbers(pathlib.Path(args.eqe))
    print(f"U(AM1.5G) = {n['u15']:.4f} A/W;  U(55 kft) = {n['u55']:.4f} A/W;  "
          f"M = {n['m']:.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plot_style.apply()
    plt.rcParams.update({"axes.grid": False})   # this sheet is typography

    # Pads below are tuned in figure-fraction at PAD_REF_H inches; scaling
    # them by PAD_REF_H/H keeps their PHYSICAL size fixed when the sheet
    # height changes, so the layout stays stable instead of re-tuning.
    FIG_H, PAD_REF_H = 42.0, 20.0
    K = PAD_REF_H / FIG_H
    fig = plt.figure(figsize=(13.0, FIG_H))
    fig.patch.set_facecolor("white")
    R = fig.canvas.get_renderer()
    FH = fig.bbox.height
    INK, MUTED, ACC = plot_style.INK, plot_style.MUTED, plot_style.ACCENT
    RED, GRN = plot_style.DEMAND, plot_style.BATT
    LEFT, RIGHT = 0.055, 0.945

    y = [0.985]

    def measure(t):
        return t.get_window_extent(renderer=R).height / FH

    def txt(s, size=13.5, color=INK, x=LEFT, ha="left", weight=None,
            pad=0.0055 * 1.0):
        """Draw at the cursor, then advance past the element's real height."""
        t = fig.text(x, y[0], s, fontsize=size, color=color, va="top", ha=ha,
                     fontweight=weight)
        y[0] -= measure(t) + pad * K
        return t

    def eq(s, size=18, color=ACC, pad=0.016, pre=0.012):
        y[0] -= pre * K
        t = fig.text(0.5, y[0], s, fontsize=size, color=color, va="top",
                     ha="center")
        y[0] -= measure(t) + pad * K
        return t

    def head(s, num):
        y[0] -= 0.014 * K
        t = fig.text(LEFT, y[0], f"{num}  {s}", fontsize=15.5, color=INK,
                     va="top", fontweight="bold")
        fig.text(LEFT, y[0], f"{num}", fontsize=15.5, color=RED, va="top",
                 fontweight="bold")
        y[0] -= measure(t) + 0.010 * K

    def rule(pad=0.011):
        y[0] -= pad * K
        fig.add_artist(plt.Line2D([LEFT, RIGHT], [y[0], y[0]],
                                  color="#dcdce4", lw=1.0))
        y[0] -= pad * K

    # ---------------- header ----------------
    t = fig.text(0.5, y[0], "The Spectral Mismatch Factor  $M$", fontsize=26,
                 ha="center", va="top", color=INK, fontweight="bold")
    y[0] -= measure(t) + 0.008 * K
    txt("standard name: spectral mismatch factor (IEC 60904-7) -- also called "
        "the spectral correction factor, or spectral factor",
        size=12.5, color=MUTED, x=0.5, ha="center")
    rule()

    # ---------------- 1 ----------------
    head("The assumption being made", "1.")
    txt("The model turns irradiance into power with one broadband number and "
        "one datasheet efficiency:")
    eq(r"$P_{avail} \;=\; G \cdot A \cdot \eta_{STC} \cdot T_{enc}$")
    txt(r"$\eta_{STC} = 25.4\%$ is measured under AM1.5G. Applying it to a "
        r"different spectrum is only valid if that spectrum")
    txt("produces the same current per watt as AM1.5G does. At altitude it "
        "does not.")
    rule()

    # ---------------- 2 ----------------
    head("A cell counts photons, not watts", "2.")
    txt(r"Spectral irradiance $E(\lambda)$ [W m$^{-2}$ nm$^{-1}$] is the power "
        r"arriving per unit wavelength. Integrating it throws")
    txt("away the colour -- and that integral is exactly what POA / GHI / DNI "
        "are:")
    eq(r"$G \;=\; \int_0^{\infty}\;E(\lambda)\, d\lambda"
       r"\qquad \mathrm{[W\,m^{-2}]}$")
    txt(r"A photon at $\lambda$ carries $hc/\lambda$ joules, so the arriving "
        r"photon flux is")
    eq(r"$\Phi(\lambda) \;=\; E(\lambda)\,\frac{\lambda}{hc}"
       r"\qquad \mathrm{[photons\;s^{-1}\,m^{-2}\,nm^{-1}]}$")
    txt(r"EQE$(\lambda)$ -- the external quantum efficiency, i.e. your "
        r"spreadsheet -- is the fraction of those photons that")
    txt("become collected electrons. Photocurrent is therefore a weighted "
        "photon count:")
    eq(r"$J_{sc} \;=\; q \int\;\Phi(\lambda)\,\mathrm{EQE}(\lambda)\,"
       r"d\lambda \;=\; \int\;E(\lambda)\,SR(\lambda)\,d\lambda$")
    txt("where the spectral response folds the constants together:")
    eq(r"$SR(\lambda) \;=\; \frac{q\lambda}{hc}\,\mathrm{EQE}(\lambda)"
       r"\;=\; \frac{\mathrm{EQE}(\lambda)\cdot\lambda_{nm}}{1239.84}"
       r"\qquad \mathrm{[A\,W^{-1}]}$")
    rule()

    # ---------------- 3 ----------------
    head("The usable fraction: amps per broadband watt", "3.")
    eq(r"$U(E) \;=\; \frac{\int\;E(\lambda)\,SR(\lambda)\,d\lambda}"
       r"{\int\;E(\lambda)\,d\lambda} \;=\; \frac{J_{sc}}{G}"
       r"\qquad \mathrm{[A\,W^{-1}]}$")
    txt(r"$U$ is scale-invariant, $U(cE) = U(E)$, so it describes the SHAPE of "
        r"the spectrum and nothing about how bright")
    txt(r"it is. Two skies delivering identical W m$^{-2}$ but different "
        r"colour have different $U$.")
    rule()

    # ---------------- 4 ----------------
    head("$M$ is the ratio of two usable fractions", "4.")
    y[0] -= 0.014 * K
    t = fig.text(0.5, y[0], r"$M \;=\; \frac{U(E_{flight})}{U(E_{AM1.5G})}$",
                 fontsize=21, color=RED, va="top", ha="center")
    h = measure(t)
    fig.patches.append(plt.Rectangle((0.33, y[0] - h - 0.008 * K), 0.34,
                                     h + 0.016 * K, transform=fig.transFigure,
                                     zorder=0, facecolor="#fdf3f1",
                                     edgecolor=RED, lw=1.6))
    y[0] -= h + 0.020 * K
    txt("Substituting it back shows exactly what the model was leaving out:")
    eq(r"$J_{sc} \;=\; G\cdot U(E_{flight}) \;=\; "
       r"\left[\, G\cdot U(E_{AM1.5G}) \,\right] \;\cdot\; M$", pad=0.006)
    txt(r"the bracket is what the model already computes; $M$ is the factor it "
        r"omits", size=12, color=MUTED, x=0.5, ha="center")
    y[0] -= 0.008 * K
    txt(r"$M = 1$ at AM1.5G by construction, so the whole correction is one "
        r"multiplication:")
    eq(r"$P_{corrected} \;=\; G \cdot A \cdot \eta_{STC} \cdot T_{enc} "
       r"\cdot M$", color=RED)
    txt(r"Equivalently: the cell behaves as though it saw $G \cdot M$ of "
        r"AM1.5G-coloured light.")
    rule()

    # ---------------- 5: the graphical derivation ----------------
    head("Where the EQE enters, and where $M$ comes from -- drawn", "5.")
    txt("Read this top to bottom. Only the middle panel comes from your "
        "spreadsheet; everything else is the two")
    txt("spectra and their pointwise product with it.")
    y[0] -= 0.010 * K

    AX_H, GAP = 0.052, 0.013           # per-panel height, figure fraction
    xlim = (300, 2500)

    def add_panel(hfrac=AX_H):
        bottom = y[0] - hfrac
        ax = fig.add_axes([0.115, bottom, 0.80, hfrac])
        y[0] = bottom - GAP * K
        ax.set_xlim(*xlim)
        ax.grid(True, color="#e6e6ee", lw=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=10)
        return ax

    # (a) the two spectra -- area under each IS G
    axA = add_panel()
    axA.fill_between(n["grid"], n["e15"], color=plot_style.SOLAR, alpha=0.30)
    axA.plot(n["grid"], n["e15"], color=plot_style.SOLAR, lw=1.3,
             label=f"$E_{{AM1.5G}}(\\lambda)$   area $= G = {n['g15']:.0f}$ W/m$^2$")
    axA.plot(n["grid"], n["e55"], color=plot_style.DEMAND, lw=1.3,
             label=f"$E_{{55kft}}(\\lambda)$   area $= G = {n['g55']:.0f}$ W/m$^2$")
    axA.set_ylabel("$E$\n[W/m$^2$/nm]", fontsize=10)
    axA.legend(fontsize=9.5, loc="upper right")
    axA.set_xticklabels([])
    axA.set_title("(a) the two spectra -- the area under each is the broadband "
                  "irradiance $G$", fontsize=11.5, pad=4)

    # (b) SR from the EQE sheet -- the ONLY input from the spreadsheet
    axB = add_panel()
    axB.fill_between(n["eqe_wl"], n["eqe"], color=plot_style.BATT, alpha=0.22)
    axB.plot(n["eqe_wl"], n["eqe"], color=plot_style.BATT, lw=1.6,
             label="measured EQE (bare cell)")
    axB.axvline(CUTOFF_NM, color=INK, ls=":", lw=1.1)
    axB.annotate(f"{CUTOFF_NM:.0f} nm: EQE $\\rightarrow 0$\n(silicon bandgap)",
                 xy=(CUTOFF_NM, 0.5), xytext=(1500, 0.62), fontsize=9.5,
                 color=INK, arrowprops=dict(arrowstyle="->", lw=0.9, color=INK))
    axB.set_ylabel("EQE\n[-]", fontsize=10)
    axB.set_ylim(0, 1.15)
    axB.legend(fontsize=9.5, loc="upper left")
    axB.set_xticklabels([])
    axB.set_title(r"(b) YOUR SPREADSHEET: EQE$(\lambda)$, converted to "
                  r"$SR(\lambda) = \mathrm{EQE}\cdot\lambda/1239.84$",
                  fontsize=11.5, pad=4)

    # (c) the pointwise product -- area under each IS Jsc
    axC = add_panel()
    axC.fill_between(n["grid"], n["e15"] * n["srg"], color=plot_style.SOLAR,
                     alpha=0.30)
    axC.plot(n["grid"], n["e15"] * n["srg"], color=plot_style.SOLAR, lw=1.3,
             label=f"$E_{{AM1.5G}}\\cdot SR$   area $= J_{{sc}} = "
                   f"{n['j15']:.1f}$ A/m$^2$")
    axC.plot(n["grid"], n["e55"] * n["srg"], color=plot_style.DEMAND, lw=1.3,
             label=f"$E_{{55kft}}\\cdot SR$   area $= J_{{sc}} = "
                   f"{n['j55']:.1f}$ A/m$^2$")
    axC.axvspan(CUTOFF_NM, xlim[1], color=INK, alpha=0.06)
    axC.text(1840, axC.get_ylim()[1] * 0.55,
             "everything here is zero:\nthe cell cannot convert it,\n"
             "but panel (a) still counts it in $G$",
             fontsize=9.5, ha="center", va="center", color=INK)
    axC.set_ylabel("$E\\cdot SR$\n[A/m$^2$/nm]", fontsize=10)
    axC.legend(fontsize=9.5, loc="upper right")
    axC.set_xlabel("Wavelength (nm)", fontsize=11)
    axC.set_title(r"(c) product $E(\lambda)\,SR(\lambda)$ -- the area under "
                  r"each is the photocurrent $J_{sc}$", fontsize=11.5, pad=4)

    y[0] -= 0.048 * K      # clear panel (c)'s tick labels and xlabel
    txt(f"Shown to 2500 nm; the integrals run to 4000 nm, where a further "
        f"{100*n['beyond']:.1f}% of AM1.5G energy sits (all of it",
        size=12, color=MUTED)
    txt("unusable). Divide panel (c) by panel (a) to get $U$, then divide the "
        "two $U$ values to get $M$:", size=12, color=MUTED)
    rule()

    # ---------------- 6 ----------------
    head("The four numbers, and $M$", "6.")
    y[0] -= 0.006 * K
    box_top = y[0]
    rows = [("", "AM1.5G (lab)", "55 kft (flight)"),
            (r"broadband  $G$", f"{n['g15']:.0f} W m$^{{-2}}$",
             f"{n['g55']:.0f} W m$^{{-2}}$"),
            (r"photocurrent  $J_{sc}$", f"{n['j15']:.1f} A m$^{{-2}}$",
             f"{n['j55']:.1f} A m$^{{-2}}$"),
            (r"usable fraction  $U$", f"{n['u15']:.4f} A/W",
             f"{n['u55']:.4f} A/W"),
            (f"energy beyond {CUTOFF_NM:.0f} nm (unusable)",
             f"{100*n['x15']:.1f}\\%", f"{100*n['x55']:.1f}\\%")]
    for i, (a, b, c) in enumerate(rows):
        y[0] -= 0.014 * K
        w = "bold" if i == 0 else None
        col = MUTED if i == 0 else INK
        t = fig.text(0.105, y[0], a, fontsize=13.5, color=col, va="top",
                     fontweight=w)
        fig.text(0.60, y[0], f"${b}$" if "\\%" in b else b, fontsize=13.5,
                 color=col, va="top", ha="right", fontweight=w)
        fig.text(0.895, y[0], f"${c}$" if "\\%" in c else c, fontsize=13.5,
                 color=col, va="top", ha="right", fontweight=w)
        y[0] -= measure(t)
    y[0] -= 0.010 * K
    fig.patches.append(plt.Rectangle((0.085, y[0]), 0.83, box_top - y[0],
                                     transform=fig.transFigure, zorder=0,
                                     facecolor="#f5f5fb", edgecolor="#dcdce4"))
    eq(rf"$M \;=\; {n['u55']:.4f} \,/\, {n['u15']:.4f} \;=\; "
       rf"\mathbf{{{n['m']:.3f}}}$", size=20, color=RED, pre=0.018)
    txt(f"The mechanism in one line: the share of arriving energy the cell "
        f"cannot convert grows from {100*n['x15']:.1f}% to "
        f"{100*n['x55']:.1f}%.")
    txt(f"Beyond {CUTOFF_NM:.0f} nm a photon lacks the energy to lift an "
        f"electron across the silicon bandgap, so that irradiance")
    txt(r"adds to the DENOMINATOR of $U$ and contributes nothing to the "
        r"numerator.", color=GRN)
    rule()

    # ---------------- 6 ----------------
    head("Why this is exactly the IEC 60904-7 factor", "7.")
    txt("The standard defines mismatch between a test device and a reference "
        "device across two spectra:")
    eq(r"$MM \;=\; \frac{\int\;E_t\,SR_t}{\int\;E_r\,SR_t} \;\cdot\; "
       r"\frac{\int\;E_r\,SR_r}{\int\;E_t\,SR_r}$")
    txt(r"Here the 'reference device' is the irradiance model itself, which is "
        r"spectrally flat -- $SR_r = c$, an ideal")
    txt(r"thermopile. Setting $SR_r = c$ cancels $c$ and collapses the four "
        r"integrals to two:")
    eq(r"$MM \;\rightarrow\; \frac{\int\;E_t SR_t \,/\, \int\;E_t}"
       r"{\int\;E_r SR_t \,/\, \int\;E_r} \;=\; \frac{U(E_t)}{U(E_r)} \;=\; M$")
    txt(r"So $M$ is not an ad-hoc fudge factor: it is the standard spectral "
        r"mismatch factor, with the broadband")
    txt("irradiance model playing the role of the reference cell.")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=165, facecolor="white")
    print(f"Saved -> {out}   (content ends at y={y[0]:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
