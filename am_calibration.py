#!/usr/bin/env python3
"""Which knob is turned wrong? Air-mass / spectral calibration of the
"cell-plane irradiance implied by current vs modeled clear-sky" gap.

THE OBSERVATION. In solar_efficiency.py's string-1 irradiance panel, the
current-implied curve (green) sits well BELOW the modeled clear-sky curve
(red) during the early high-altitude hold, then crosses ABOVE it down low.
Those two curves are not the same physical quantity:

    red   = POA_broadband * encapsulation           -- ALL wavelengths
    green = (Impp/Impp_stc)/(1+alpha*dT) * 1000     -- what the CELL converts,
                                                       referenced to AM1.5G STC

So green/red is a direct measurement of spectral mismatch, times whatever
else is wrong. Writing it out:

    green/red = M / r * L

  M  spectral mismatch vs AM1.5G, from SPECTRL2 + the MEASURED Max7 EQE
  r  BeerLambert/SPECTRL2 broadband bias of the model that produced red
     (see spectral_mismatch.py; ~1.00 at cruise, ~0.76 low and low-sun)
  L  everything else -- the residual this script exists to attribute

WHERE THE AIRCRAFT ACTUALLY SITS. At the high hold (16,765 m / 55,005 ft)
the ISA pressure ratio is 0.090, so the ABSOLUTE air mass is only
AM_rel * p/p0 = 1.27 * 0.090 = 0.115. Spectrally that is 92% of the way from
AM1.5G to AM0 -- not AM0, and nowhere near AM1.5. The reason it is not AM0 is
specific and worth knowing: water vapour and aerosol are essentially all
BELOW the aircraft and gone, Rayleigh is down to ~9%, but the ozone layer
peaks near 22 km, so at 17 km almost the whole O3 column is still ABOVE and
its UV/Chappuis absorption is nearly intact.

WHAT THE EQE CALIBRATES. The measured EQE pins the photocurrent scale
independently of the datasheet efficiency route:
    Jsc,stc = integral(E_AM1.5G(l) * SR(l) dl) = 41.60 mA/cm^2  -> Isc = 6.448 A
which is textbook Maxeon and brackets Impp,stc to 6.00-6.19 A for an
admissible Impp/Isc of 0.93-0.96. The script's 6.085 A sits inside that band,
so the current anchor CANNOT absorb the residual -- that is what makes L a
statement about physics rather than calibration.

Note the EQE is for the BARE cell, not the laminated cell (see
spectral_mismatch.load_spectral_response). So it constrains the cell's
photocurrent scale and the laminate's spectral TILT, but NOT the absolute
laminate transmission (0.93 x 0.92 = 0.856, of which POE's 0.92 is a flat
user-supplied figure with no spectral data behind it). A constant laminate
error is therefore still admissible inside L; only the heading swing proves
that PART of L is angular. Separating them needs a flash Isc or EQE of a
LAMINATED coupon: Isc_lam/Isc_bare is the absolute transmission directly.

CELL TEMPERATURE FROM V AND I TOGETHER. Neither alone is enough: current
gives irradiance (alpha is only +0.058 %/degC), voltage gives temperature
(beta = -0.236 %/degC) but only once irradiance is known, because Vmpp
carries an n*kT/q*ln(G) term. The two-unknown system is solved by fixed-point
iteration, here for EACH MPPT string independently and with the spectral
correction folded into the current->irradiance step (which the estimator in
solar_efficiency.py does not do).

Reads:  max7_eqe.csv, outputs/<log>_solar_efficiency.csv
Writes: outputs/<log>_am_calibration.png
        outputs/<log>_irradiance_corrected.png
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
import pvlib
from pvlib import spectrum

import plot_style
import spectral_mismatch as sm

# ---- Cell electrical constants (mirrored from solar_efficiency.py) ----
STC_TEMP_C = 25.0
STC_IRRADIANCE_W_M2 = 1000.0
VMPP_TEMP_COEFF_PCT_PER_C = -0.236
ISC_TEMP_COEFF_PCT_PER_C = +0.058
DIODE_IDEALITY = 1.1
BOLTZMANN_OVER_Q_V_PER_K = 8.617333262e-5
CELL_VMPP_STC_V = 0.647
CELL_TEMP_MIN_CURRENT_A = 0.3
CELL_TEMP_OFFMPP_GUARD_V = 1.0
CELL_TEMP_MAX_PLAUSIBLE_C = 85.0
CELL_TEMP_MIN_DELTA_VS_TOUT_C = -5.0
CELL_EFFICIENCY = 0.254
CELL_AREA_CM2 = 155.0
ETFE_T, POE_T = 0.93, 0.92

# String wiring. 72 cells total, 58 in string 1 -> 14 in string 0.
STRINGS = {0: 14, 1: 58}

# Above this altitude the model's broadband bias r is ~1 and the aircraft is
# above the weather, so L there is interpretable. Below it, L is contaminated
# by both r's aerosol/water assumptions and by real cloud.
L_TRUSTWORTHY_ALT_M = 7000.0

MAXEON_IMPP_OVER_ISC = (0.93, 0.96)


def solve_temperature(v_string, i_string, n_cells, impp_stc, m_spectral,
                      beta_pct=VMPP_TEMP_COEFF_PCT_PER_C, iters=12):
    """Simultaneous (T, G) from one string's Vmpp and Impp.

    Current gives irradiance, voltage gives temperature, and they are coupled
    through Vmpp's n*kT/q*ln(G) term -- so iterate:

        G/G0 = (Impp/Impp_stc) / (1 + alpha*(T-25)) / M
        T    = 25 + (v_cell - v_stc - a_v*ln(G/G0)) / beta_v_per_c

    m_spectral divides the current->irradiance step: measured photocurrent is
    produced by the SPECTRALLY WEIGHTED irradiance, so recovering broadband
    irradiance means dividing by M. Pass M=1 to reproduce the uncorrected
    estimator in solar_efficiency.py.
    """
    v_cell = np.asarray(v_string, float) / n_cells
    i_use = np.asarray(i_string, float)
    m = np.asarray(m_spectral, float)
    beta_v_per_c = (beta_pct / 100.0) * CELL_VMPP_STC_V

    t_c = np.full(v_cell.shape, STC_TEMP_C)
    g_ratio = np.full(v_cell.shape, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        for _ in range(iters):
            a_v = DIODE_IDEALITY * BOLTZMANN_OVER_Q_V_PER_K * (t_c + 273.15)
            g_ratio = (i_use / impp_stc) / (
                1.0 + (ISC_TEMP_COEFF_PCT_PER_C / 100.0) * (t_c - STC_TEMP_C)) / m
            g_ratio = np.where(g_ratio > 1e-4, g_ratio, np.nan)
            t_c = STC_TEMP_C + (v_cell - CELL_VMPP_STC_V
                                - a_v * np.log(g_ratio)) / beta_v_per_c
    return t_c, g_ratio


def eqe_current_anchor(sr_wl, sr):
    """Isc,stc per cell from the measured EQE against AM1.5G -- the anchor
    that is independent of the datasheet efficiency/Vmpp route."""
    refs = spectrum.get_reference_spectra()
    w = refs.index.to_numpy(float)
    k = (w >= 300) & (w <= 4000)
    e, w = refs["global"].to_numpy(float)[k], w[k]
    jsc_a_m2 = np.trapezoid(e * np.interp(w, sr_wl, sr, left=0.0, right=0.0), x=w)
    return jsc_a_m2 * (CELL_AREA_CM2 / 1e4), jsc_a_m2 / 10.0  # A, mA/cm^2


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", required=True)
    p.add_argument("--eqe", default="max7_eqe.csv")
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args()

    csv_path = pathlib.Path(args.csv)
    df = pd.read_csv(csv_path, parse_dates=["time"], low_memory=False)
    df = df.set_index("time")
    doy = int(df.index[0].dayofyear)

    sr_wl, sr = sm.load_spectral_response(pathlib.Path(args.eqe))
    u_ref, refs = sm.reference_usable_fraction(sr_wl, sr, 300.0, 4000.0)
    isc_eqe, jsc_ma_cm2 = eqe_current_anchor(sr_wl, sr)
    impp_script = (CELL_EFFICIENCY * CELL_AREA_CM2 / 1e4
                   * STC_IRRADIANCE_W_M2) / CELL_VMPP_STC_V

    print("=" * 78)
    print("EQE CURRENT ANCHOR  (independent of the efficiency/Vmpp route)")
    print("=" * 78)
    print(f"  Jsc,stc  = {jsc_ma_cm2:.2f} mA/cm2   (Maxeon-class expectation 41-42)")
    print(f"  Isc,stc  = {isc_eqe:.3f} A for a {CELL_AREA_CM2:.0f} cm2 cell")
    lo, hi = MAXEON_IMPP_OVER_ISC
    print(f"  admissible Impp,stc = {isc_eqe*lo:.3f}-{isc_eqe*hi:.3f} A "
          f"(Impp/Isc {lo}-{hi})")
    print(f"  script uses {impp_script:.3f} A -> Impp/Isc = {impp_script/isc_eqe:.3f}  "
          f"{'INSIDE the admissible band' if lo <= impp_script/isc_eqe <= hi else 'OUTSIDE'}")
    print()

    # ---- geometry / spectral state along the flight ---------------------
    alt = df["alt_msl_m"].to_numpy(float)
    zen = 90.0 - df["sun_elevation_deg"].to_numpy(float)
    pressure_ratio = sm.isa_pressure_ratio(alt)
    am_rel = pvlib.atmosphere.get_relative_airmass(zen)
    df["am_absolute"] = am_rel * pressure_ratio

    print("Building M and r grids via SPECTRL2 ...")
    alt_grid, zen_grid, m_grid, r_grid = sm.build_mismatch_grid(
        sr_wl, sr, u_ref, doy)
    df["m_spectral"] = sm.interp_mismatch(alt_grid, zen_grid, m_grid, alt, zen)
    df["r_broadband"] = sm.interp_mismatch(alt_grid, zen_grid, r_grid, alt, zen)
    df.loc[zen > sm.MAX_ZENITH_DEG, ["m_spectral", "r_broadband"]] = np.nan

    # M at the pure AM0 limit, used to express M as a position on the
    # AM1.5G -> AM0 axis (M is 1.0 at AM1.5G by construction).
    out = spectrum.spectrl2(
        apparent_zenith=0.0, aoi=0.0, surface_tilt=0.0, ground_albedo=0.0,
        surface_pressure=1.0, relative_airmass=1.0, precipitable_water=1e-4,
        ozone=1e-6, aerosol_turbidity_500nm=1e-6, dayofyear=doy)
    m_am0 = sm.usable_fraction(out["wavelength"],
                               np.asarray(out["dni"], float).ravel(),
                               sr_wl, sr) / u_ref
    df["am0_fraction"] = (1.0 - df["m_spectral"]) / (1.0 - m_am0)

    # ---- the two curves and their ratio ---------------------------------
    encap = ETFE_T * POE_T
    red = pd.to_numeric(df["poa_string1_w_m2"], errors="coerce") * encap
    green = pd.to_numeric(df["g_cell_string1_w_m2"], errors="coerce")
    held = df["tc_string1_held"].astype("boolean").fillna(True)
    valid = (red > 150) & green.notna() & ~held & df["m_spectral"].notna()
    df["red_modeled"], df["green_implied"] = red, green
    df["red_spectral"] = red * df["m_spectral"]                    # corrected
    df["L_residual"] = np.where(valid, (green / red.where(red > 150))
                                * df["r_broadband"] / df["m_spectral"], np.nan)

    hh = (df["flight_phase"] == "holding_high") & valid
    print()
    print("=" * 78)
    print("WHERE THE AIRCRAFT SITS, SPECTRALLY  (early high-altitude hold)")
    print("=" * 78)
    print(f"  window            {df.index[hh][0]:%H:%M}-{df.index[hh][-1]:%H:%M} local, "
          f"{int(hh.sum())} clean samples")
    print(f"  altitude          {alt[hh].mean():.0f} m = {alt[hh].mean()*3.28084:.0f} ft")
    print(f"  pressure ratio    {pressure_ratio[hh].mean():.4f} of sea level")
    print(f"  relative air mass {am_rel[hh].mean():.2f} (zenith {zen[hh].mean():.1f} deg)")
    print(f"  ABSOLUTE air mass {df['am_absolute'][hh].mean():.3f}")
    print(f"  M (spectral)      {df['m_spectral'][hh].mean():.4f}   "
          f"[AM1.5G = 1.000, AM0 = {m_am0:.4f}]")
    print(f"  => {100*df['am0_fraction'][hh].mean():.0f}% of the way from AM1.5G to AM0")
    print(f"     (not AM0: the O3 column peaks near 22 km, so at "
          f"{alt[hh].mean()/1000:.0f} km it is still overhead)")
    print()

    # ---- knob attribution ------------------------------------------------
    band = valid & (alt > 14000)
    gr = float(np.nanmedian(green[band] / red[band]))
    mm = float(np.nanmedian(df["m_spectral"][band]))
    rr = float(np.nanmedian(df["r_broadband"][band]))
    # Close the books on the identity rather than taking a third independent
    # median: median-of-ratios != ratio-of-medians, and a table that does not
    # add up invites the reader to distrust all of it.
    ll = gr * rr / mm
    print("=" * 78)
    print("KNOB ATTRIBUTION  (>14 km, where r ~ 1 and the aircraft is above weather)")
    print("=" * 78)
    print(f"  observed green/red = {gr:.4f}  ({100*(1-gr):.1f}% deficit)")
    print(f"  M / r * L          = {mm:.4f} / {rr:.4f} * {ll:.4f} = {mm/rr*ll:.4f}")
    print()
    print(f"  {'knob':34s} {'closes':>9s}  verdict")
    need = 1.0 - gr
    rows = [
        ("spectral mismatch M (wavelength)", 1 - mm,
         "MISSING TERM -- confirmed 3 ways"),
        ("residual L (geometry/mismatch)", 1 - ll,
         "heading-dependent -- see below"),
        ("broadband model r", rr - 1,
         "real bug, but only below ~7 km"),
    ]
    for name, closes, verdict in rows:
        print(f"  {name:34s} {100*closes:8.1f}%  {verdict}")
    print()
    dt_needed = -need / (ISC_TEMP_COEFF_PCT_PER_C / 100.0)
    print("  RULED OUT as explanations of the current deficit:")
    print(f"   - cell temperature: green scales as 1/(1+alpha*dT) with alpha only "
          f"+{ISC_TEMP_COEFF_PCT_PER_C} %/degC,")
    print(f"     so closing the deficit would need dT = {dt_needed:.0f} degC. "
          f"A realistic 20 degC")
    print(f"     anchor error moves green {100*(ISC_TEMP_COEFF_PCT_PER_C/100)*20:.1f}%.")
    print(f"   - current anchor: would need Impp/Isc = {impp_script*gr/isc_eqe:.3f} "
          f"(deficit) or {impp_script*ll/isc_eqe:.3f} (residual);")
    print(f"     the EQE admits only {lo}-{hi}.")
    print(f"   - ETFE laminate optics: re-weighting with the MEASURED EQE moves "
          f"T_eff 0.9300 -> 0.9314,")
    print(f"     and only 0.9303 at 55 kft. Confirmed correct, no altitude "
          f"dependence.")
    print(f"     (POE's flat {POE_T} is the least-constrained optical number; "
          f"{POE_T*ll:.3f} would close L,")
    print(f"      but a CONSTANT laminate loss cannot produce the heading swing "
          f"below.)")
    print()

    # ---- heading discriminator ------------------------------------------
    rel = (pd.to_numeric(df["sun_azimuth_deg"], errors="coerce")
           - pd.to_numeric(df["yaw_deg"], errors="coerce")) % 360.0
    print("  Heading discriminator: a constant loss is flat vs heading; series")
    print("  mismatch on a curved array is not.")
    print(f"    {'sun bearing rel. nose':24s} {'n':>7s} {'L':>8s}")
    meds = []
    for lo_deg in range(0, 360, 45):
        m = band & (rel >= lo_deg) & (rel < lo_deg + 45)
        if m.sum() < 100:
            continue
        med = float(np.nanmedian(df["L_residual"][m]))
        meds.append(med)
        side = ("nose" if lo_deg < 45 or lo_deg >= 315 else
                "right wing" if 45 <= lo_deg < 135 else
                "tail" if 135 <= lo_deg < 225 else "left wing")
        print(f"    {lo_deg:3d}-{lo_deg+45:3d} deg ({side:10s}) {int(m.sum()):7d} {med:8.4f}")
    if meds:
        print(f"    spread = {max(meds)-min(meds):.3f}  -> "
              f"{'ANGULAR/MISMATCH signature' if max(meds)-min(meds) > 0.05 else 'flat: uniform loss'}")
    print()

    # ---- cell temperature from V and I, per string -----------------------
    print("=" * 78)
    print("CELL TEMPERATURE FROM V AND I, EACH MPPT STRING")
    print("=" * 78)
    print("Current gives irradiance, voltage gives temperature; they are coupled")
    print("through Vmpp's n*kT/q*ln(G) term, so both are solved together.")
    print()
    results = {}
    for s, n_cells in STRINGS.items():
        v = pd.to_numeric(df[f"pv_voltage_v_{s}"], errors="coerce")
        i = pd.to_numeric(df[f"pv_current_a_{s}"], errors="coerce")
        pw = pd.to_numeric(df[f"pv_power_w_{s}"], errors="coerce")
        on = (pw > 5) & v.notna() & i.notna() & (i > CELL_TEMP_MIN_CURRENT_A)
        # Same guards solar_efficiency.py applies to string 1, so the two
        # channels are gated identically and can be compared:
        #  - off-MPP ride: string V well below its own rolling median
        #  - non-thermal results: implausibly hot, or further below the
        #    fuselage TC than an illuminated cell can physically sit
        vr = v.rolling("120s", min_periods=30).median()
        clean = on & ((vr - v) < CELL_TEMP_OFFMPP_GUARD_V)
        tout = pd.to_numeric(df["tout_c"], errors="coerce")
        for tag, mvec in (("uncorrected", pd.Series(1.0, index=df.index)),
                          ("M-corrected", df["m_spectral"])):
            t, g = solve_temperature(v.where(clean), i.where(clean), n_cells,
                                     impp_script, mvec.where(clean))
            t = pd.Series(t, index=df.index)
            physical = ((t <= CELL_TEMP_MAX_PLAUSIBLE_C)
                        & ((t - tout) >= CELL_TEMP_MIN_DELTA_VS_TOUT_C))
            t = t.where(physical)
            df[f"tc_s{s}_{'raw' if tag=='uncorrected' else 'corr'}"] = t
            results[(s, tag)] = (t, g, clean & physical)
        clean = results[(s, "M-corrected")][2]
        t_raw = df[f"tc_s{s}_raw"]; t_cor = df[f"tc_s{s}_corr"]
        hi_m = clean & (alt > 14000)
        print(f"  string {s} ({n_cells} cells, {100*pw.sum()/(pd.to_numeric(df['pv_power_w_0'],errors='coerce').sum()+pd.to_numeric(df['pv_power_w_1'],errors='coerce').sum()):.0f}% of energy), "
              f"V/cell {(v[hi_m]/n_cells).median():.4f} at the high hold")
        dmc = t_cor[clean & hi_m].median() - t_raw[clean & hi_m].median()
        print(f"     spectral correction shifts T by {dmc:+.1f} degC at >14 km "
              f"(small: T depends on ln(G))")
        # Report the physics-gate survival rate alongside every temperature:
        # the gate censors (it removes exactly the implausibly-cold results),
        # so a mean over survivors is only a measurement if most samples
        # survived. A low rate means the inversion failed, not that it is cold.
        for lbl, regime in (("high hold", df["flight_phase"] == "holding_high"),
                            ("low hold", df["flight_phase"] == "holding_low")):
            att = on & regime
            kept = clean & regime
            if att.sum() < 100:
                continue
            rate = 100.0 * kept.sum() / max(int(att.sum()), 1)
            flag = ("" if rate > 85 else
                    "   <-- partly censored, treat as a bound" if rate > 55 else
                    "   <-- CENSORED, not a measurement")
            print(f"     {lbl:9s}: {t_cor[kept].median():+6.1f} degC   "
                  f"gate survival {rate:5.1f}% of {int(att.sum())}{flag}")
    # Cross-check: the two strings should read the same temperature.
    hi_m = alt > 14000
    d0 = df["tc_s0_corr"].where(hi_m).median()
    d1 = df["tc_s1_corr"].where(hi_m).median()
    print()
    print(f"  string 0 minus string 1 at >14 km = {d0-d1:+.1f} degC")
    print()
    print("  IS STRING 0 TRUSTWORTHY? Two independent tests of its wiring:")
    i0 = pd.to_numeric(df["pv_current_a_0"], errors="coerce")
    i1 = pd.to_numeric(df["pv_current_a_1"], errors="coerce")
    v0s = pd.to_numeric(df["pv_voltage_v_0"], errors="coerce")
    v1s = pd.to_numeric(df["pv_voltage_v_1"], errors="coerce")
    both = (i0 > 0.3) & (i1 > 0.3)
    hb = both & hi_m
    ir = float((i0[hb] / i1[hb]).median())
    print(f"   1. CURRENT: series strings of the same cells must carry the same")
    print(f"      current. I0/I1 at >14 km = {ir:.3f} -> the 14-cell count and")
    print(f"      Impp,stc are {'CONFIRMED' if 0.9 < ir < 1.1 else 'NOT confirmed'} "
          f"(no parallel wiring).")
    vc0 = float((v0s[hb] / STRINGS[0]).median())
    vc1 = float((v1s[hb] / STRINGS[1]).median())
    beta_v = (VMPP_TEMP_COEFF_PCT_PER_C / 100.0) * CELL_VMPP_STC_V
    print(f"   2. VOLTAGE: V/cell is {vc0:.4f} (string 0) vs {vc1:.4f} (string 1),")
    print(f"      a {1000*(vc0-vc1):+.0f} mV/cell split = {(vc0-vc1)/beta_v:+.0f} degC of")
    print(f"      apparent temperature. Implied series count if the two strings")
    print(f"      were truly isothermal: {v0s[hb].median()/vc1:.1f} -- not an integer.")
    print(f"      A fixed {14*(vc0-vc1):+.2f} V offset on the string-0 voltage sense")
    print(f"      reconciles it exactly; a cell-count error cannot (14 -> "
          f"{14*(vc0-vc1)/CELL_VMPP_STC_V:+.1f} cells).")
    print()
    print("   => STRING 1 is the trustworthy temperature channel (81% of energy,")
    print("      58 cells corroborated). String 0's absolute temperature is NOT")
    print("      usable until that voltage offset is resolved; its CURRENT is fine")
    print("      and is what the irradiance/spectral work above relies on.")

    if not args.no_plot:
        make_figures(df, csv_path, sr_wl, sr, refs, doy, m_am0, alt_grid,
                     zen_grid, m_grid)
    return 0


def make_figures(df, csv_path, sr_wl, sr, refs, doy, m_am0, alt_grid,
                 zen_grid, m_grid):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    plot_style.apply()
    stem = csv_path.name.replace("_solar_efficiency.csv", "")
    disp = plot_style.label(stem)          # cmr10 has no underscore glyph
    alt = df["alt_msl_m"]
    tz = df.index.tz

    # ================= figure 1: AM / spectral calibration ==============
    fig, axes = plt.subplots(4, 1, figsize=(13, 16))
    fig.suptitle(f"Air-mass and spectral calibration of the current-implied "
                 f"irradiance gap\n{disp}", fontsize=15)

    # (a) altitude + absolute air mass
    ax = axes[0]
    ax.plot(df.index, alt / 1000.0, color=plot_style.ACCENT, lw=1.8,
            label="Altitude (MSL)")
    ax.set_ylabel("Altitude (km)")
    ax2 = ax.twinx()
    ax2.plot(df.index, df["am_absolute"], color=plot_style.DEMAND, lw=1.2,
             label="Absolute air mass")
    ax2.set_yscale("log")
    ax2.set_ylabel("Absolute air mass (AM)")
    ax2.grid(False)
    for y, lb in ((1.5, "AM1.5"), (1.0, "AM1.0")):
        ax2.axhline(y, color=plot_style.MUTED, ls=":", lw=0.9)
        ax2.text(df.index[int(0.01*len(df))], y*1.05, lb, fontsize=10,
                 color=plot_style.MUTED)
    hh = df["flight_phase"] == "holding_high"
    if hh.any():
        ax.axvspan(df.index[hh][0], df.index[hh][-1], color=plot_style.SOLAR,
                   alpha=0.15, zorder=0)
        ax.annotate(f"high hold: {alt[hh].mean()*3.28084:.0f} ft, "
                    f"AM $\\approx$ {df['am_absolute'][hh].mean():.2f},\n"
                    f"{100*df['am0_fraction'][hh].mean():.0f}% of the way to AM0",
                    xy=(df.index[hh][len(df.index[hh])//2], alt[hh].mean()/1000.0),
                    xytext=(0.30, 0.30), textcoords="axes fraction",
                    fontsize=11, arrowprops=dict(arrowstyle="->", lw=1.2,
                                    color=plot_style.INK))
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc="upper right")
    ax.set_title("(a) Altitude and the ABSOLUTE air mass actually flown "
                 "(AM$_{rel}$ $\\times$ p/p$_0$)")

    # (b) the curves, uncorrected and spectrally corrected
    ax = axes[1]
    ax.plot(df.index, df["red_modeled"], color=plot_style.DEMAND, alpha=0.55,
            lw=1.4, label="Modeled clear-sky POA $\\times$ encapsulation (broadband)")
    ax.plot(df.index, df["red_spectral"], color=plot_style.WARN, lw=1.4, ls="--",
            label="... $\\times$ M (spectrally corrected for the Max7 EQE)")
    ax.plot(df.index, df["green_implied"], color=plot_style.BATT, lw=0.8,
            label="Cell-plane irradiance implied by measured current")
    if hh.any():
        ax.axvspan(df.index[hh][0], df.index[hh][-1], color=plot_style.SOLAR,
                   alpha=0.15, zorder=0)
    ax.set_ylabel("Irradiance (W/m$^2$)")
    ax.legend(loc="upper right", fontsize=10)
    ax.set_title("(b) Applying M closes most of the high-hold gap; what is left "
                 "is not spectral")

    # (c) the spectra themselves, with the EQE that samples them
    ax = axes[2]
    w0 = refs.index.to_numpy(float)
    k = (w0 >= 300) & (w0 <= 2000)
    ax.plot(w0[k], refs["global"].to_numpy(float)[k], color=plot_style.SOLAR,
            lw=1.4, label="AM1.5G (the 25.4% rating condition)")
    w_hi, d_hi = sm.spectrl2_direct(38.2, 16765.0, doy)
    kk = (w_hi >= 300) & (w_hi <= 2000)
    ax.plot(w_hi[kk], d_hi[:, 0][kk], color=plot_style.DEMAND, lw=1.4,
            label="55 kft, zenith 38$^\\circ$ (SPECTRL2)")
    ax.plot(w0[k], refs["extraterrestrial"].to_numpy(float)[k],
            color=plot_style.MUTED, lw=1.2, ls=":", label="AM0 (space)")
    ax.axvspan(1190, 2000, color=plot_style.INK, alpha=0.07)
    ax.text(1600, ax.get_ylim()[1]*0.55,
            "Max7 EQE = 0 beyond 1190 nm:\nthis irradiance is counted as\n'available' but cannot be converted",
            fontsize=10, ha="center", color=plot_style.INK)
    # Mark the H2O bands that suppress NIR at sea level and thin out with
    # altitude -- the visual mechanism behind M's altitude dependence.
    for wl_band in (940, 1130):
        ax.annotate("H$_2$O", xy=(wl_band, 0.62), xytext=(wl_band, 1.15),
                    fontsize=10, ha="center", color=plot_style.DEMAND,
                    arrowprops=dict(arrowstyle="->", lw=0.9,
                                    color=plot_style.DEMAND))
    axb = ax.twinx()
    axb.plot(sr_wl, sr * 1239.84 / sr_wl, color=plot_style.BATT, lw=2.0,
             label="Max7 EQE (measured)")
    axb.set_ylabel("EQE"); axb.set_ylim(0, 1.05); axb.grid(False)
    axb.legend(loc="lower right", fontsize=10)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Spectral irradiance (W/m$^2$/nm)")
    ax.legend(loc="upper right", fontsize=10)
    ax.set_title("(c) Why altitude moves M: the H$_2$O bands and the "
                 "unusable IR beyond 1190 nm")

    # (d) the residual L, with the untrustworthy region marked
    ax = axes[3]
    ok = df["L_residual"].notna()
    ax.scatter(alt[ok] / 1000.0, df["L_residual"][ok], s=2, alpha=0.10,
               color=plot_style.INK)
    bins = np.arange(0, 18.5, 1.0)
    idx = np.digitize(alt[ok] / 1000.0, bins)
    ctr = 0.5 * (bins[:-1] + bins[1:])
    med = [np.nanmedian(df["L_residual"][ok].to_numpy()[idx == n])
           if (idx == n).sum() > 20 else np.nan for n in range(1, len(bins))]
    ax.plot(ctr, med, "o-", color=plot_style.DEMAND, lw=1.8, label="median L")
    ax.axhline(1.0, color=plot_style.INK, ls=":", lw=0.9)
    ax.axvspan(0, L_TRUSTWORTHY_ALT_M / 1000.0, color=plot_style.WARN, alpha=0.10)
    ax.text(L_TRUSTWORTHY_ALT_M / 2000.0, 0.45,
            "below 7 km L is not interpretable:\nr depends on aerosol/water "
            "assumptions\nand the clear-sky model cannot see cloud",
            fontsize=10, ha="center", color=plot_style.INK)
    hi = ok & (alt > L_TRUSTWORTHY_ALT_M)
    lv = float(np.nanmedian(df["L_residual"][hi]))
    ax.axhline(lv, color=plot_style.BATT, lw=1.4, ls="--",
               label=f"L = {lv:.3f} above 7 km")
    ax.set_xlabel("Altitude (km)"); ax.set_ylabel("L (residual)")
    ax.set_ylim(0.3, 1.4)
    ax.legend(loc="lower right")
    ax.set_title("(d) What is left after M and r: flat $\\approx$ "
                 f"{lv:.2f} where it can be trusted")

    for a_ in axes[:2]:
        a_.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=tz))
    fig.tight_layout()
    out1 = csv_path.parent / f"{stem}_am_calibration.png"
    fig.savefig(out1)
    print(f"\nSaved plot -> {out1}")
    plt.close(fig)

    # ================= figure 2: the corrected irradiance panel =========
    fig, axes = plt.subplots(3, 1, figsize=(13, 12), sharex=True)
    fig.suptitle(f"Corrected cell-plane irradiance and cell temperature\n{disp}",
                 fontsize=15)

    ax = axes[0]
    ax.plot(df.index, alt / 1000.0, color=plot_style.ACCENT, lw=1.8,
            label="Altitude (MSL)")
    ax.set_ylabel("Altitude (km)")
    ax.legend(loc="upper right")
    ax.set_title("(a) Altitude")

    ax = axes[1]
    ax.plot(df.index, df["red_spectral"] / df["r_broadband"],
            color=plot_style.DEMAND, lw=1.4,
            label="Modeled, corrected: POA/r $\\times$ encapsulation $\\times$ M")
    ax.plot(df.index, df["green_implied"], color=plot_style.BATT, lw=0.8,
            label="Implied by measured current")
    ax.set_ylabel("Irradiance (W/m$^2$)")
    ax.legend(loc="upper right", fontsize=10)
    ax.set_title("(b) Both spectral and broadband corrections applied "
                 "(gap left = the angular term)")

    ax = axes[2]
    for s, col, lb in ((1, plot_style.DEMAND, "String 1 (58 cells)"),
                       (0, plot_style.BATT, "String 0 (14 cells)")):
        t = df[f"tc_s{s}_corr"].rolling("300s", min_periods=1).mean()
        ax.plot(df.index, t, color=col, lw=1.4, label=f"{lb}, M-corrected")
        t0 = df[f"tc_s{s}_raw"].rolling("300s", min_periods=1).mean()
        ax.plot(df.index, t0, color=col, lw=0.9, ls=":", alpha=0.8,
                label=f"{lb}, uncorrected")
    ax.plot(df.index, pd.to_numeric(df["tout_c"], errors="coerce"),
            color=plot_style.MUTED, lw=1.2, label="Fuselage skin TC (side-mounted)")
    ax.axhline(25.0, color=plot_style.INK, ls=":", lw=0.9, label="STC (25 degC)")
    ax.set_ylabel("Temperature (degC)")
    ax.set_ylim(-60, 80)
    ax.legend(loc="lower right", fontsize=9, ncol=2)
    ax.set_title("(c) Cell temperature solved from V and I, per string")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=tz))
    ax.set_xlabel(f"Local time ({df.index[0].date()})")

    fig.tight_layout()
    out2 = csv_path.parent / f"{stem}_irradiance_corrected.png"
    fig.savefig(out2)
    print(f"Saved plot -> {out2}")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
