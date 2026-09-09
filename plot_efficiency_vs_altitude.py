#!/usr/bin/env python3
"""True array aperture efficiency as a function of altitude.

APERTURE EFFICIENCY here means DC watts measured at the panel side (pre-MPPT)
per broadband watt of sunlight landing on the cell aperture (1.116 m^2 of
cells, not airframe planform):

    eta(h) = P_measured / (POA_true * A_cells)

POA_true divides the model's POA by r, the gray-Beer-Lambert airmass bias
(see spectral_mismatch.py) -- otherwise the denominator is wrong by up to 25%
low down and the "efficiency" absorbs a modelling error.

Alongside the measurement it plots the predicted chain

    eta(h) = eta_cell * T_enc * M(h, z) * L * f_temp(T_cell(h))

so the two can be compared with nothing fitted per-altitude except the cell
temperature, which is read from the flight's own Vmpp-derived estimate.

IMPORTANT: below ~7 km the measurement is NOT an array property. The
denominator is a CLEAR-SKY model that cannot see cloud, so any cloud shows up
as apparent inefficiency. That band is shaded and excluded from conclusions.

Usage:
    python plot_efficiency_vs_altitude.py --csv outputs/<log>_solar_efficiency.csv
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

import plot_style
import spectral_mismatch as sm

CELL_AREA_M2 = 155.0 / 1e4
N_STRING0, N_STRING1 = 14, 58
ETA_CELL, ENCAP = 0.254, 0.93 * 0.92
L_RESIDUAL = 0.918          # measured, flat above 7 km (am_calibration.py)
GAMMA_PCT_PER_C = -0.27
STC_TEMP_C = 25.0
CLOUD_FLOOR_M = 7000.0      # below this the clear-sky denominator is unreliable


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--eqe", default="max7_eqe.csv")
    ap.add_argument("--bin-km", type=float, default=1.5)
    args = ap.parse_args()

    csv_path = pathlib.Path(args.csv)
    df = pd.read_csv(csv_path, parse_dates=["time"], low_memory=False).set_index("time")
    doy = int(df.index[0].dayofyear)

    sr_wl, sr = sm.load_spectral_response(pathlib.Path(args.eqe))
    u_ref, _ = sm.reference_usable_fraction(sr_wl, sr, 300.0, 4000.0)
    alt_grid, zen_grid, m_grid, r_grid = sm.build_mismatch_grid(
        sr_wl, sr, u_ref, doy)

    alt = df["alt_msl_m"].to_numpy(float)
    zen = 90.0 - df["sun_elevation_deg"].to_numpy(float)
    M = sm.interp_mismatch(alt_grid, zen_grid, m_grid, alt, zen)
    Rb = sm.interp_mismatch(alt_grid, zen_grid, r_grid, alt, zen)

    a0, a1 = N_STRING0 * CELL_AREA_M2, N_STRING1 * CELL_AREA_M2
    p = (pd.to_numeric(df["pv_power_w_0"], errors="coerce")
         + pd.to_numeric(df["pv_power_w_1"], errors="coerce"))
    sun_model = (pd.to_numeric(df["poa_string0_w_m2"], errors="coerce") * a0
                 + pd.to_numeric(df["poa_string1_w_m2"], errors="coerce") * a1)
    sun_true = sun_model / Rb
    tcell = pd.to_numeric(df["tc_string1_est_c"], errors="coerce")

    eta = (p / sun_true).to_numpy(float)
    f_temp = 1.0 + (GAMMA_PCT_PER_C / 100.0) * (tcell.to_numpy(float) - STC_TEMP_C)
    eta_pred = ETA_CELL * ENCAP * M * L_RESIDUAL * f_temp

    ok = (np.isfinite(eta) & np.isfinite(Rb) & (sun_model.to_numpy(float) > 50)
          & (p.to_numpy(float) > 5))

    print(f"cell aperture: {N_STRING0}+{N_STRING1} cells = {a0 + a1:.3f} m2")
    print("APERTURE EFFICIENCY vs ALTITUDE  (DC pre-MPPT / broadband on cells)\n")
    print(f"{'band':13s} {'n':>7s} {'measured':>10s} {'predicted':>10s} "
          f"{'M':>7s} {'f temp':>7s} {'T cell':>8s}")
    edges = np.arange(0.0, 18.0 + args.bin_km, args.bin_km)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = ok & (alt >= lo * 1000) & (alt < hi * 1000)
        if k.sum() < 200:
            continue
        e_m = float(np.nanmedian(eta[k]))
        e_p = float(np.nanmedian(eta_pred[k]))
        flag = "" if lo * 1000 >= CLOUD_FLOOR_M else "   <- cloud-contaminated"
        print(f"{lo:4.1f}-{hi:4.1f} km {int(k.sum()):7d} {100*e_m:9.1f}% "
              f"{100*e_p:9.1f}% {np.nanmedian(M[k]):7.3f} "
              f"{np.nanmedian(f_temp[k]):7.3f} "
              f"{np.nanmedian(tcell.to_numpy(float)[k]):+7.1f}C{flag}")
        rows.append((0.5 * (lo + hi), e_m, e_p, k.sum()))

    hi_k = ok & (alt > CLOUD_FLOOR_M)
    print(f"\n  ABOVE {CLOUD_FLOOR_M/1000:.0f} km (trustworthy): measured "
          f"{100*np.nanmedian(eta[hi_k]):.1f}%   predicted "
          f"{100*np.nanmedian(eta_pred[hi_k]):.1f}%")
    print(f"  installed STC ceiling (cell x encapsulation) = "
          f"{100*ETA_CELL*ENCAP:.1f}%")

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plot_style.apply()

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(12, 10.5),
                                  gridspec_kw=dict(height_ratios=[1.5, 1.0]),
                                  sharex=True)
    ctr = np.array([r[0] for r in rows])
    meas = 100 * np.array([r[1] for r in rows])
    pred = 100 * np.array([r[2] for r in rows])

    ax.scatter(alt[ok] / 1000.0, 100 * eta[ok], s=1.5, alpha=0.06,
               color=plot_style.INK)
    ax.plot(ctr, meas, "o-", color=plot_style.DEMAND, lw=2.2, ms=7,
            label="measured aperture efficiency (binned median)")
    ax.plot(ctr, pred, "s--", color=plot_style.ACCENT, lw=2.0, ms=6,
            label=r"predicted: $\eta_{cell}\cdot T_{enc}\cdot M \cdot L "
                  r"\cdot f_{temp}$")
    ax.axhline(100 * ETA_CELL * ENCAP, color=plot_style.SOLAR, ls=":", lw=1.8,
               label=f"installed STC ceiling {100*ETA_CELL*ENCAP:.1f}% "
                     f"(cell $\\times$ encapsulation)")
    ax.axhline(100 * ETA_CELL, color=plot_style.MUTED, ls=":", lw=1.4,
               label=f"bare cell datasheet {100*ETA_CELL:.1f}%")
    ax.axvspan(0, CLOUD_FLOOR_M / 1000.0, color=plot_style.WARN, alpha=0.10)
    ax.text(CLOUD_FLOOR_M / 2000.0, 6.0,
            "below 7 km the denominator is a CLEAR-SKY model:\n"
            "cloud shows up here as apparent inefficiency,\n"
            "so these points are not an array property",
            fontsize=10.5, ha="center", color=plot_style.INK)
    ax.set_ylabel("Aperture efficiency (%)")
    ax.set_ylim(0, 28)
    ax.legend(loc="lower right", fontsize=10.5)
    ax.set_title("True array efficiency vs altitude -- flat above the weather, "
                 "because M and temperature cancel")

    # why it is flat
    axb = ax2
    axb.plot(ctr, [np.nanmedian(M[ok & (np.abs(alt / 1000.0 - c) <= args.bin_km / 2)])
                   for c in ctr], "o-", color=plot_style.BATT, lw=2.0,
             label="$M$ (spectral) -- falls with altitude")
    axb.plot(ctr, [np.nanmedian(f_temp[ok & (np.abs(alt / 1000.0 - c) <= args.bin_km / 2)])
                   for c in ctr], "s-", color=plot_style.DEMAND, lw=2.0,
             label="$f_{temp}$ (cold cells) -- rises with altitude")
    prod = [np.nanmedian((M * f_temp)[ok & (np.abs(alt / 1000.0 - c) <= args.bin_km / 2)])
            for c in ctr]
    axb.plot(ctr, prod, "^-", color=plot_style.INK, lw=2.4,
             label="product -- nearly flat")
    axb.axhline(1.0, color=plot_style.MUTED, ls=":", lw=1.2)
    axb.axvspan(0, CLOUD_FLOOR_M / 1000.0, color=plot_style.WARN, alpha=0.10)
    axb.set_xlabel("Altitude (km)")
    axb.set_ylabel("Correction factor")
    axb.legend(loc="lower left", fontsize=10.5)
    axb.set_title("The two altitude-dependent terms pull opposite ways")

    fig.tight_layout()
    out = csv_path.parent / (csv_path.name.replace("_solar_efficiency.csv", "")
                             + "_efficiency_vs_altitude.png")
    fig.savefig(out)
    print(f"\nSaved plot -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
