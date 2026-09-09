#!/usr/bin/env python3
"""Two panels, nothing else: altitude vs time, and a zoom on the high-altitude
hold showing that the spectral mismatch factor M closes most of the gap between
the modeled clear-sky irradiance and the irradiance implied by measured current.

Usage:
    python plot_m_simple.py --csv outputs/<log>_solar_efficiency.csv
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

import plot_style
import spectral_mismatch as sm

ENCAPSULATION = 0.93 * 0.92     # ETFE x POE, as solar_efficiency.py applies it


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--eqe", default="max7_eqe.csv")
    args = ap.parse_args()

    csv_path = pathlib.Path(args.csv)
    df = pd.read_csv(csv_path, parse_dates=["time"], low_memory=False).set_index("time")
    doy = int(df.index[0].dayofyear)

    # M on an (altitude, zenith) grid, then looked up per sample.
    sr_wl, sr = sm.load_spectral_response(pathlib.Path(args.eqe))
    u_ref, _ = sm.reference_usable_fraction(sr_wl, sr, 300.0, 4000.0)
    alt_grid, zen_grid, m_grid, _ = sm.build_mismatch_grid(sr_wl, sr, u_ref, doy)
    alt = df["alt_msl_m"].to_numpy(float)
    zen = 90.0 - df["sun_elevation_deg"].to_numpy(float)
    m_spectral = sm.interp_mismatch(alt_grid, zen_grid, m_grid, alt, zen)

    modeled = pd.to_numeric(df["poa_string1_w_m2"], errors="coerce") * ENCAPSULATION
    corrected = modeled * m_spectral
    measured = pd.to_numeric(df["g_cell_string1_w_m2"], errors="coerce")

    hold = df["flight_phase"] == "holding_high"
    if not hold.any():
        print("no high-altitude hold in this flight")
        return 1
    t0, t1 = df.index[hold][0], df.index[hold][-1]

    # Percent error of each modeled curve against the measurement.
    # error = model/measured - 1 ; positive means the model reads HIGH.
    err_raw = (modeled / measured - 1.0) * 100.0
    err_corr = (corrected / measured - 1.0) * 100.0

    # Gap statistics inside the hold, on samples where all three exist.
    k = hold & modeled.notna() & measured.notna() & np.isfinite(m_spectral) & (modeled > 150)
    m_hold = float(np.nanmedian(m_spectral[k]))
    before, after = float(np.nanmedian(err_raw[k])), float(np.nanmedian(err_corr[k]))
    print(f"High hold {t0:%H:%M}-{t1:%H:%M}, {int(k.sum())} samples, "
          f"{alt[hold].mean()*3.28084:.0f} ft,  M = {m_hold:.3f}\n")
    print(f"  {'':28s} {'median':>9s} {'mean':>8s} {'RMS':>8s} {'p10..p90':>16s}")
    for tag, e in (("modeled broadband vs meas.", err_raw[k]),
                   ("modeled x M vs measured", err_corr[k])):
        e = e.to_numpy(float)
        e = e[np.isfinite(e)]
        print(f"  {tag:28s} {np.median(e):+8.1f}% {e.mean():+7.1f}% "
              f"{np.sqrt((e**2).mean()):7.1f}% "
              f"{np.percentile(e,10):+7.1f}..{np.percentile(e,90):+.1f}%")
    print(f"\n  HEADLINE: {before:+.1f}%  ->  {after:+.1f}%   "
          f"(M removes {before-after:.1f} of the {before:.1f} points of bias)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    plot_style.apply()

    fig, (ax0, ax1, ax2) = plt.subplots(
        3, 1, figsize=(12, 12.5),
        gridspec_kw=dict(height_ratios=[1.0, 1.3, 1.1]))
    tz = df.index.tz

    # ---- (a) altitude vs time -------------------------------------------
    ax0.plot(df.index, alt / 1000.0, color=plot_style.ACCENT, lw=2.2)
    ax0.axvspan(t0, t1, color=plot_style.SOLAR, alpha=0.25)
    ax0.text(t0 + (t1 - t0) / 2, alt[hold].mean() / 1000.0 - 3.0,
             f"high hold\n{alt[hold].mean()*3.28084:.0f} ft",
             ha="center", fontsize=12, color=plot_style.INK)
    ax0.set_ylabel("Altitude (km)")
    ax0.set_title("(a) Altitude")
    ax0.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=tz))

    # ---- (b) zoom on the hold -------------------------------------------
    z = df.loc[t0:t1]
    zc = corrected.loc[t0:t1]
    ax1.plot(z.index, modeled.loc[t0:t1], color=plot_style.DEMAND, lw=2.0,
             label=f"Modeled, broadband  (reads {before:+.0f}% high)")
    ax1.plot(z.index, zc, color=plot_style.ACCENT, lw=2.2, ls="--",
             label=f"Modeled $\\times$ M = {m_hold:.3f}  (reads {after:+.0f}% high)")
    ax1.plot(z.index, measured.loc[t0:t1], color=plot_style.BATT, lw=1.2,
             label="Measured (implied by cell current)")
    ax1.set_ylabel("Cell-plane irradiance (W/m$^2$)")
    ax1.set_xlim(t0, t1)
    ax1.set_title(f"(b) Zoom on the hold: M closes {before-after:.0f} of the "
                  f"{before:.0f} point gap")
    ax1.legend(loc="lower left", fontsize=11)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=tz))
    ax1.tick_params(labelbottom=False)

    # ---- (c) percent error vs time, same window --------------------------
    ax2.plot(z.index, err_raw.loc[t0:t1], color=plot_style.DEMAND, lw=1.2,
             label=f"Modeled broadband  (median {before:+.1f}%)")
    ax2.plot(z.index, err_corr.loc[t0:t1], color=plot_style.ACCENT, lw=1.2,
             label=f"Modeled $\\times$ M  (median {after:+.1f}%)")
    ax2.axhline(0.0, color=plot_style.INK, lw=1.2)
    ax2.axhline(before, color=plot_style.DEMAND, ls=":", lw=1.2)
    ax2.axhline(after, color=plot_style.ACCENT, ls=":", lw=1.4)
    ax2.set_ylabel("Model error (%)")
    ax2.set_ylim(-15, 45)
    ax2.set_xlim(t0, t1)
    ax2.set_title("(c) The same thing as percent error: M shifts the whole "
                  "trace down, the oscillation is the racetrack legs")
    ax2.legend(loc="upper right", fontsize=11)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=tz))
    ax2.set_xlabel(f"Local time ({df.index[0].date()})")

    fig.tight_layout()
    out = csv_path.parent / f"{csv_path.name.replace('_solar_efficiency.csv','')}_M_simple.png"
    fig.savefig(out)
    print(f"\nSaved plot -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
