#!/usr/bin/env python3
"""Does spectral mismatch (+ cell temperature) explain the modeled-vs-measured gap?

solar_efficiency.py models available power as

    P_avail = POA_broadband * area * eta_25.4% * T_encap

which embeds two assumptions that break at stratospheric altitude:

  1. SPECTRAL. POA is broadband (all wavelengths), and the Beer-Lambert
     attenuation SEA_LEVEL_TRANSMITTANCE**pressure_ratio is GRAY -- the same
     transmittance at every wavelength. The real atmosphere attenuates
     selectively (Rayleigh ~lambda^-4, H2O bands at 940/1130/1380/1870 nm,
     O3 in the UV), so climbing does not just scale the spectrum up, it
     changes its SHAPE. Max7 EQE is ~unity over 600-900 nm and hard-zero
     past ~1190 nm, while AM1.5G carries a large fraction of its energy
     beyond 1200 nm. Broadband irradiance the cell cannot convert is
     therefore counted as "available" -- and the unusable share is exactly
     what changes with altitude, because the H2O bands that suppress it at
     sea level thin out with pressure and humidity.

     Correction: the spectral mismatch factor
         U(E) = integral(E(l) * SR(l) dl) / integral(E(l) dl)
         M    = U(actual) / U(AM1.5G)
     where SR(l) = EQE(l) * q * l / (h*c) is the spectral response [A/W].
     M multiplies photocurrent exactly and power to first order.

  2. TEMPERATURE. eta = 25.4% is a 25 degC rating and the model does not
     derate by default (POWER_TEMP_COEFF_PCT_PER_C is "informational"). Real
     cells run ~41 degC at low hold (model over-predicts) and ~0 degC at high
     hold (model UNDER-predicts). The two effects fight each other at
     altitude, so they have to be evaluated together, not one at a time.

On the user's "use current to decipher what temperature everything is at":
the current/voltage roles are actually inverted. Impp carries alpha_Isc =
+0.058 %/degC -- far too weak to be a thermometer -- but it is very nearly
proportional to SPECTRALLY-WEIGHTED irradiance, which makes it an excellent
independent probe of M. Vmpp carries beta = -0.236 %/degC and is the real
thermometer, which is what estimate_string1_cell_temperature() already uses.
So this script uses Vmpp-derived cell temperature for the temperature leg and
measured string-1 CURRENT as the empirical check on the spectral leg.

Modeled spectra come from pvlib's SPECTRL2 (Bird & Riordan) driven by ISA
pressure at the aircraft's own altitude. Two validation anchors are checked
before any result is reported (see validate()).

Reads:  max7_eqe.csv, outputs/<log>_solar_efficiency.csv
Writes: outputs/<log>_spectral_mismatch.png + a decomposition table on stdout
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pvlib
from pvlib import spectrum

# ---- Constants mirrored from solar_efficiency.py (single source of truth
# there; duplicated here so this script stays standalone/read-only). ----
STC_TEMP_C = 25.0
POWER_TEMP_COEFF_PCT_PER_C = -0.27    # Pmpp, 546209 Rev C / EQE workbook Sheet1
ISC_TEMP_COEFF_PCT_PER_C = +0.058     # Isc
CELL_EFFICIENCY = 0.254
CELL_AREA_CM2 = 155.0
CELL_VMPP_STC_V = 0.647
STRING1_CELL_COUNT = 58

# Atmospheric profile assumptions for SPECTRL2. Water vapour is the dominant
# driver of the effect under test, so its profile is the most load-bearing
# assumption in this script and is reported explicitly.
PW_SEA_LEVEL_CM = 1.5          # precipitable water at sea level
PW_SCALE_HEIGHT_M = 2000.0     # H2O scale height (troposphere, ~2 km)
AOD500_SEA_LEVEL = 0.08        # aerosol optical depth at 500 nm
AOD_SCALE_HEIGHT_M = 1500.0    # aerosol scale height
OZONE_ATM_CM = 0.30            # column O3; ~flat below the 22 km O3 peak
SEA_LEVEL_PRESSURE_PA = 101325.0

MAX_ZENITH_DEG = 80.0          # Bird's model is not trusted at grazing sun

# SPECTRL2 spans 300-4000 nm; ~1.3% of the AM0 integral lies outside that
# window, so its integrated DNI is divided by this to compare against a
# true broadband model.
SPECTRL2_BAND_COVERAGE = 0.987


def isa_pressure_ratio(alt_m):
    """ISA static pressure as a fraction of sea level. Ported verbatim from
    solar_efficiency.py so both scripts see the same atmosphere."""
    alt_m = np.asarray(alt_m, dtype=float)
    density = np.where(
        alt_m <= 11000.0,
        1.225 * ((288.15 - 0.0065 * alt_m) / 288.15) ** 4.25587,
        np.where(
            alt_m <= 20000.0,
            0.363918 * np.exp(-0.000157688 * (alt_m - 11000.0)),
            0.088035 * ((216.65 + 0.001 * (alt_m - 20000.0)) / 216.65) ** -35.1632,
        ),
    )
    temperature_k = np.where(
        alt_m <= 11000.0,
        288.15 - 0.0065 * alt_m,
        np.where(alt_m <= 20000.0, 216.65, 216.65 + 0.001 * (alt_m - 20000.0)),
    )
    return (density / 1.225) * (temperature_k / 288.15)


# --------------------------------------------------------------------------
# Spectral response and the mismatch integral
# --------------------------------------------------------------------------
def load_spectral_response(path: Path):
    """SR(lambda) [A/W] from measured external QE of the BARE CELL.

    Provenance matters here (user, 2026-09-02): Max7_EQE_Data_Only.xlsx is the
    bare cell, NOT the laminated cell. That is the right input for both uses:

      * M is a ratio with the same SR in numerator and denominator, so a
        laminate transmission that is FLAT in wavelength cancels out of it
        EXACTLY -- injecting a flat 0.8586 into SR moves M by 0.00%. Only the
        laminate's spectral TILT survives, and measured against the digitized
        curves that is bare 0.9114 vs 0.9098 (ETFE x POE uv-through) or 0.9042
        (uv-cut) at 55 kft, i.e. -0.18% / -0.78%. (An earlier revision of this
        note said 0.07%; that figure was ETFE-only and predates
        poe_transmittance.csv, so it understated the uv-cut case ~10x.)
        Either way the ~14% the laminate ABSORBS is absent from M, which is
        why solar_efficiency.py stacks M on top of the encapsulation factors
        rather than substituting it for them.
      * The absolute Isc it yields must pair with a bare-cell Impp_stc, which
        is what solar_efficiency.py derives from the bare-cell datasheet
        efficiency. Against bare Isc the implied Impp/Isc is 0.944 (Maxeon
        0.93-0.96); against a laminated Isc it would be 1.097, non-physical.

    The laminate belongs on the MODELED side (POA x encapsulation), not here.
    Consequence: this sheet cannot validate the ABSOLUTE laminate
    transmission (0.93 x 0.92), so any laminate loss lands in the residual L.

    EQE is already external (the workbook's reflectance column is baked in),
    so it is used as-is."""
    eqe = pd.read_csv(path).dropna(subset=["wavelength_nm", "eqe"])
    wl = eqe["wavelength_nm"].to_numpy(float)
    # SR = EQE * q*lambda/(h*c); 1239.84 nm*eV is the hc/q conversion.
    sr = eqe["eqe"].to_numpy(float) * wl / 1239.84
    return wl, sr


def usable_fraction(wavelength_nm, spectral_irradiance, sr_wl, sr):
    """U = integral(E*SR) / integral(E). Non-uniform wavelength grids are the
    norm here (SPECTRL2 returns 122 unevenly spaced points), so integrate
    with an explicit x= rather than assuming constant spacing."""
    sr_on_grid = np.interp(wavelength_nm, sr_wl, sr, left=0.0, right=0.0)
    e = np.asarray(spectral_irradiance, dtype=float)
    total = np.trapezoid(e, x=wavelength_nm)
    if not np.isfinite(total) or total <= 0:
        return np.nan
    return float(np.trapezoid(e * sr_on_grid, x=wavelength_nm) / total)


def reference_usable_fraction(sr_wl, sr, lo_nm, hi_nm):
    """U for the datasheet's AM1.5G global-tilt reference (ASTM G173),
    truncated to the same bounds as the SPECTRL2 grid so numerator and
    denominator integration limits match."""
    refs = spectrum.get_reference_spectra()
    wl = refs.index.to_numpy(float)
    m = (wl >= lo_nm) & (wl <= hi_nm)
    return usable_fraction(wl[m], refs["global"].to_numpy(float)[m], sr_wl, sr), refs


def spectrl2_direct(zenith_deg, alt_m, dayofyear):
    """SPECTRL2 direct-normal spectral irradiance at the aircraft.

    Direct beam (not global) is the right comparison: solar_efficiency.py
    models POA from DNI with dhi forced to zero, so the flight sees an
    essentially pure direct beam. The datasheet reference stays AM1.5G
    global -- that lab-vs-flight difference is itself part of the mismatch.
    """
    zenith_deg = np.atleast_1d(np.asarray(zenith_deg, dtype=float))
    alt_m = np.atleast_1d(np.asarray(alt_m, dtype=float))
    pressure_pa = SEA_LEVEL_PRESSURE_PA * isa_pressure_ratio(alt_m)
    out = spectrum.spectrl2(
        apparent_zenith=zenith_deg,
        aoi=zenith_deg,                 # panel normal to the sun; only the
        surface_tilt=zenith_deg,        # direct component is used anyway
        ground_albedo=0.2,
        surface_pressure=pressure_pa,
        relative_airmass=pvlib.atmosphere.get_relative_airmass(zenith_deg),
        precipitable_water=PW_SEA_LEVEL_CM * np.exp(-alt_m / PW_SCALE_HEIGHT_M),
        ozone=OZONE_ATM_CM,
        aerosol_turbidity_500nm=AOD500_SEA_LEVEL * np.exp(-alt_m / AOD_SCALE_HEIGHT_M),
        dayofyear=dayofyear,
    )
    return out["wavelength"], np.atleast_2d(out["dni"].T).T


def build_mismatch_grid(sr_wl, sr, u_ref, dayofyear,
                        alt_grid=None, zen_grid=None):
    """M and the gray-model bias ratio on an (altitude x zenith) grid.

    SPECTRL2 per flight sample would be 76k spectra; the grid + bilinear
    interpolation is smooth in both axes and ~400 spectra instead.

    Returns (alt_grid, zen_grid, m_grid, r_grid) where r_grid is
    DNI_beer_lambert / DNI_spectrl2 -- the altitude/zenith-dependent bias of
    solar_efficiency.py's own broadband irradiance model. r<1 means that
    model UNDER-predicts irradiance, so it under-states available power and
    inflates the reported efficiency.
    """
    if alt_grid is None:
        alt_grid = np.concatenate([np.arange(0, 3000, 500.0),
                                   np.arange(3000, 19001, 1000.0)])
    if zen_grid is None:
        zen_grid = np.arange(0.0, MAX_ZENITH_DEG + 0.1, 5.0)
    aa, zz = np.meshgrid(alt_grid, zen_grid, indexing="ij")
    wl, dni = spectrl2_direct(zz.ravel(), aa.ravel(), dayofyear)
    m = np.array([usable_fraction(wl, dni[:, i], sr_wl, sr) / u_ref
                  for i in range(dni.shape[1])]).reshape(aa.shape)

    s2_broadband = np.array([np.trapezoid(dni[:, i], x=wl)
                             for i in range(dni.shape[1])]) / SPECTRL2_BAND_COVERAGE
    # solar_efficiency.py's model, evaluated on the same grid.
    toa = float(pvlib.irradiance.get_extra_radiation(dayofyear))
    cz = np.cos(np.radians(zz.ravel()))
    tau = 0.70 ** isa_pressure_ratio(aa.ravel())
    bl = toa * tau ** (1.0 / np.where(cz > 0, cz, 1.0))
    r = (bl / s2_broadband).reshape(aa.shape)
    return alt_grid, zen_grid, m, r


def interp_mismatch(alt_grid, zen_grid, m_grid, alt_m, zenith_deg):
    """Bilinear lookup into the M grid, clipped to grid bounds."""
    a = np.clip(np.asarray(alt_m, float), alt_grid[0], alt_grid[-1])
    z = np.clip(np.asarray(zenith_deg, float), zen_grid[0], zen_grid[-1])
    ai = np.clip(np.searchsorted(alt_grid, a) - 1, 0, len(alt_grid) - 2)
    zi = np.clip(np.searchsorted(zen_grid, z) - 1, 0, len(zen_grid) - 2)
    fa = (a - alt_grid[ai]) / (alt_grid[ai + 1] - alt_grid[ai])
    fz = (z - zen_grid[zi]) / (zen_grid[zi + 1] - zen_grid[zi])
    return ((1 - fa) * (1 - fz) * m_grid[ai, zi]
            + fa * (1 - fz) * m_grid[ai + 1, zi]
            + (1 - fa) * fz * m_grid[ai, zi + 1]
            + fa * fz * m_grid[ai + 1, zi + 1])


# --------------------------------------------------------------------------
# Validation -- run before any number is believed
# --------------------------------------------------------------------------
def validate(sr_wl, sr, u_ref, refs, dayofyear) -> bool:
    """Two independent anchors:

    1. Vacuum limit: strip the atmosphere (p->0, pw->0, aod->0) and SPECTRL2's
       direct beam must collapse onto the extraterrestrial (AM0) spectrum.
       This is the check that the pressure input is doing what we think.
    2. Known physics: c-Si is rated ~10% relative lower at AM0 than AM1.5G,
       so M aloft should land near 0.88-0.93 and M at sea level near 1.0.
       A wildly different number means the integral, not the atmosphere.
    """
    ok = True
    print("=" * 72)
    print("VALIDATION")
    print("=" * 72)

    wl0 = refs.index.to_numpy(float)
    et = refs["extraterrestrial"].to_numpy(float)
    lo, hi = 300.0, 4000.0

    out = spectrum.spectrl2(
        apparent_zenith=0.0, aoi=0.0, surface_tilt=0.0, ground_albedo=0.0,
        surface_pressure=1.0, relative_airmass=1.0,
        precipitable_water=1e-4, ozone=1e-6, aerosol_turbidity_500nm=1e-6,
        dayofyear=dayofyear,
    )
    wl_s = out["wavelength"]
    dni_vac = np.asarray(out["dni"], dtype=float).ravel()
    m_et = (wl0 >= lo) & (wl0 <= hi)
    et_total = np.trapezoid(et[m_et], x=wl0[m_et])
    vac_total = np.trapezoid(dni_vac, x=wl_s)
    # AM0 total is scaled by the earth-sun distance on this dayofyear.
    et_scaled = et_total * float(pvlib.irradiance.get_extra_radiation(dayofyear) / 1366.1)
    err = 100.0 * (vac_total - et_scaled) / et_scaled
    print(f"  [1] vacuum limit: SPECTRL2 DNI {vac_total:7.1f} W/m2 vs AM0 "
          f"{et_scaled:7.1f} W/m2  -> {err:+.1f}%")
    if abs(err) > 5.0:
        print("      FAIL: vacuum limit does not recover AM0")
        ok = False
    else:
        print("      pass (<5%)")

    u_am0 = usable_fraction(wl_s, dni_vac, sr_wl, sr)
    print(f"  [2] M at AM0 limit          = {u_am0 / u_ref:.4f}   "
          f"(expect ~0.88-0.93; c-Si AM0 derate)")
    if not (0.84 <= u_am0 / u_ref <= 0.96):
        print("      FAIL: outside the physically expected AM0 band")
        ok = False
    else:
        print("      pass")

    wl_sl, dni_sl = spectrl2_direct(48.2, 0.0, dayofyear)
    m_sl = usable_fraction(wl_sl, dni_sl[:, 0], sr_wl, sr) / u_ref
    print(f"      M at sea level, AM1.5 geometry (z=48.2 deg) = {m_sl:.4f}   "
          f"(expect ~1.00)")
    if not (0.94 <= m_sl <= 1.06):
        print("      FAIL: sea-level M should reproduce the reference")
        ok = False
    else:
        print("      pass")
    print()
    return ok


# --------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", required=True, help="solar_efficiency.py output CSV")
    p.add_argument("--eqe", default="max7_eqe.csv")
    p.add_argument("--cell-vmpp-stc", type=float, default=CELL_VMPP_STC_V)
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args()

    csv_path = Path(args.csv)
    df = pd.read_csv(csv_path, parse_dates=["time"], low_memory=False)
    dayofyear = int(pd.Timestamp(df["time"].iloc[0]).dayofyear)

    sr_wl, sr = load_spectral_response(Path(args.eqe))
    print(f"Max7 EQE: {len(sr_wl)} points, {sr_wl[0]:.0f}-{sr_wl[-1]:.0f} nm; "
          f"peak SR {sr.max():.3f} A/W at {sr_wl[sr.argmax()]:.0f} nm")

    u_ref, refs = reference_usable_fraction(sr_wl, sr, 300.0, 4000.0)
    print(f"AM1.5G usable fraction U_ref = {u_ref:.4f} A/W "
          f"(i.e. {100*u_ref/ (u_ref):.0f}% baseline by construction)\n")

    if not validate(sr_wl, sr, u_ref, refs, dayofyear):
        print("Validation failed -- refusing to report mismatch numbers.")
        return 1

    print("Building M(altitude, zenith) grid via SPECTRL2 ...")
    alt_grid, zen_grid, m_grid, r_grid = build_mismatch_grid(
        sr_wl, sr, u_ref, dayofyear)
    print(f"  grid {m_grid.shape[0]} altitudes x {m_grid.shape[1]} zeniths; "
          f"M range [{np.nanmin(m_grid):.3f}, {np.nanmax(m_grid):.3f}]")
    print(f"  {'':7s} {'M(z=20)':>9s} {'M(z=40)':>9s} {'M(z=60)':>9s}   "
          f"{'r(z=20)':>9s} {'r(z=40)':>9s} {'r(z=60)':>9s}")
    for a in (0.0, 5000.0, 10000.0, 17600.0):
        zs = np.array([20.0, 40.0, 60.0])
        mm = interp_mismatch(alt_grid, zen_grid, m_grid, np.full(3, a), zs)
        rr = interp_mismatch(alt_grid, zen_grid, r_grid, np.full(3, a), zs)
        print(f"  {a/1000:5.1f} km {mm[0]:9.3f} {mm[1]:9.3f} {mm[2]:9.3f}   "
              f"{rr[0]:9.3f} {rr[1]:9.3f} {rr[2]:9.3f}")
    print("  M = spectral mismatch vs AM1.5G;  r = BeerLambert/SPECTRL2 broadband")
    print("  (r < 1 => solar_efficiency.py under-predicts irradiance there)")
    print()

    # ---- per-sample legs -------------------------------------------------
    zen = 90.0 - df["sun_elevation_deg"].to_numpy(float)
    alt = df["alt_msl_m"].to_numpy(float)
    df["m_spectral"] = interp_mismatch(alt_grid, zen_grid, m_grid, alt, zen)
    df["r_broadband"] = interp_mismatch(alt_grid, zen_grid, r_grid, alt, zen)
    df.loc[zen > MAX_ZENITH_DEG, ["m_spectral", "r_broadband"]] = np.nan

    # Temperature leg from the Vmpp-derived cell temperature (the thermometer).
    tcell = pd.to_numeric(df["tc_string1_est_c"], errors="coerce")
    df["f_temp"] = 1.0 + (POWER_TEMP_COEFF_PCT_PER_C / 100.0) * (tcell - STC_TEMP_C)

    # ---- empirical spectral probe from measured string-1 current ---------
    # Impp ~= Impp_stc * (G_eff/1000) * (1 + alpha*(T-25)); solving for the
    # ratio of effective to modeled broadband irradiance isolates M (up to a
    # single unknown current-scale constant, hence the normalisation below).
    impp_stc = (CELL_EFFICIENCY * (CELL_AREA_CM2 / 10000.0) * 1000.0) / args.cell_vmpp_stc
    i1 = pd.to_numeric(df["pv_current_a_1"], errors="coerce")
    poa1 = pd.to_numeric(df["poa_string1_w_m2"], errors="coerce")
    alpha_corr = 1.0 + (ISC_TEMP_COEFF_PCT_PER_C / 100.0) * (tcell - STC_TEMP_C)
    df["m_empirical"] = (i1 / alpha_corr) / (impp_stc * (poa1 / 1000.0))

    # Only samples where current is genuinely an irradiance measurement:
    # a fresh (non-held) cell-temp estimate means the off-MPP voltage guard
    # and low-signal/flat-current masks all passed for that sample.
    held = df["tc_string1_held"].astype("boolean").fillna(True)
    clean = (tcell.notna() & ~held & (poa1 > 200.0) & (i1 > 0.5)
             & df["m_spectral"].notna() & np.isfinite(df["m_empirical"]))
    print(f"Empirical probe: {int(clean.sum())} of {len(df)} samples usable "
          f"({100*clean.mean():.1f}%)  Impp_stc={impp_stc:.2f} A")

    # ---- decomposition by altitude band ---------------------------------
    P = pd.to_numeric(df["pv_power_actual_w"], errors="coerce")
    eff = pd.to_numeric(df["pre_mppt_efficiency_pct"], errors="coerce")
    avail = P / (eff / 100.0)                      # the model's own denominator

    bands = [(0, 3000, "< 3 km"), (3000, 10000, "3-10 km"),
             (10000, 99000, "> 10 km")]

    def ewm(mask, col):
        """Energy-weighted mean of col, weighted by modeled available power so
        it composes consistently with the energy-weighted ratio."""
        w = avail[mask]
        v = df.loc[mask, col]
        return float(np.nansum(v * w) / np.nansum(w.where(v.notna())))

    print()
    print("=" * 100)
    print("GAP DECOMPOSITION  (energy-weighted within each altitude band)")
    print("=" * 100)
    print(f"{'band':9s} {'n':>7s} {'meas/model':>11s} {'M_spec':>8s} {'f_temp':>8s} "
          f"{'r_bband':>8s} {'net':>8s} {'corrected':>10s} {'residual':>9s}")
    a = df["alt_msl_m"].to_numpy(float)
    base = P.notna() & avail.notna() & (avail > 0) & np.isfinite(avail)
    for lo, hi, lbl in list(bands) + [(0, 99000, "ALL")]:
        m = base & (a >= lo) & (a < hi)
        if m.sum() < 100:
            continue
        if lbl == "ALL":
            print("-" * 100)
        ratio = P[m].sum() / avail[m].sum()
        ms, ft, rb = ewm(m, "m_spectral"), ewm(m, "f_temp"), ewm(m, "r_broadband")
        # P_meas / true_available, where
        #   true_available = (POA_BL / r) * A * eta * T * M * f_temp
        net = ms * ft / rb
        corrected = ratio / net
        print(f"{lbl:9s} {int(m.sum()):7d} {100*ratio:10.1f}% {ms:8.3f} {ft:8.3f} "
              f"{rb:8.3f} {net:8.3f} {100*corrected:9.1f}% {100*(corrected-1):+8.1f}%")

    print()
    print("  meas/model = measured / modeled available, as reported today")
    print("  M_spec     = SPECTRL2 spectral mismatch vs AM1.5G  (<1 => model over-predicts power)")
    print("  f_temp     = 1 + gamma*(T_cell-25), gamma = "
          f"{POWER_TEMP_COEFF_PCT_PER_C} %/degC  (>1 => cold cells, model under-predicts)")
    print("  r_bband    = BeerLambert/SPECTRL2 broadband DNI  (<1 => model under-predicts irradiance)")
    print("  net        = M*f_temp/r ; corrected = meas/model divided by net.")
    print("  corrected == 100% would mean these three effects fully explain the gap.")

    # ---- empirical vs modeled M trend -----------------------------------
    print()
    print("=" * 92)
    print("EMPIRICAL CROSS-CHECK  (string-1 current as an irradiance probe)")
    print("=" * 92)
    print("Absolute scale is contaminated by the Impp_stc anchor, so columns are")
    print("normalised to their own < 3 km baseline; the TREND is the falsifiable")
    print("claim. 'emp raw' divides Impp by the Beer-Lambert POA; 'emp fixed'")
    print("divides by the SPECTRL2-corrected POA (i.e. x r), removing the")
    print("broadband model's own altitude bias from the probe.")
    df["m_empirical_fixed"] = df["m_empirical"] * df["r_broadband"]
    sub = df[clean]
    asub = sub["alt_msl_m"].to_numpy(float)
    base_m = (asub < 3000)
    if base_m.sum() > 50:
        b_raw = sub.loc[base_m, "m_empirical"].median()
        b_fix = sub.loc[base_m, "m_empirical_fixed"].median()
        b_mod = sub.loc[base_m, "m_spectral"].median()
        print(f"\n{'band':9s} {'n':>7s} {'emp raw':>9s} {'emp fixed':>10s} "
              f"{'model M':>9s} {'delta raw':>10s} {'delta fixed':>12s}")
        for lo, hi, lbl in bands:
            mm = (asub >= lo) & (asub < hi)
            if mm.sum() < 50:
                print(f"{lbl:9s} {int(mm.sum()):7d}   (too few samples)")
                continue
            e_ = sub.loc[mm, "m_empirical"].median() / b_raw
            f_ = sub.loc[mm, "m_empirical_fixed"].median() / b_fix
            o_ = sub.loc[mm, "m_spectral"].median() / b_mod
            print(f"{lbl:9s} {int(mm.sum()):7d} {e_:9.3f} {f_:10.3f} {o_:9.3f} "
                  f"{e_-o_:+10.3f} {f_-o_:+12.3f}")
        print(f"\nraw (un-normalised) empirical M at the <3 km baseline: {b_raw:.3f}")
        print("  -- offset from 1.0 is Impp_stc/telemetry calibration, not spectrum")
        print("     (the estimator's own self-calibration already flagged ~5%).")
    else:
        print("  too few clean low-altitude samples to normalise")

    # ---- is any SYSTEMATIC loss left after all three corrections? --------
    # The energy-weighted residual above mixes two very different things: a
    # calibration/systematic error (present at every sample) and time-varying
    # loss (cloud, attitude, off-MPP tracking -- present only some of the
    # time). The distribution separates them: if the CEILING of the corrected
    # ratio sits at ~100%, there is no systematic bias left and the
    # energy-weighted shortfall is entirely time-varying loss.
    print()
    print("=" * 92)
    print("IS ANY SYSTEMATIC LOSS LEFT?  distribution of measured / corrected model")
    print("=" * 92)
    ratio_corr = (P / avail) / (df["m_spectral"] * df["f_temp"] / df["r_broadband"])
    print(f"{'band':9s} {'n':>7s} {'p25':>7s} {'p50':>7s} {'p75':>7s} "
          f"{'p90':>7s} {'p95':>7s} {'p99':>7s}")
    for lo, hi, lbl in bands:
        m = ((a >= lo) & (a < hi) & np.isfinite(ratio_corr) & (P > 5)
             & tcell.notna())
        if m.sum() < 100:
            continue
        q = np.percentile(ratio_corr[m], [25, 50, 75, 90, 95, 99]) * 100
        print(f"{lbl:9s} {int(m.sum()):7d} " + " ".join(f"{v:6.1f}%" for v in q))
    print()
    print("  A p90-p95 near 100% means the corrected model is an accurate")
    print("  CEILING -- no systematic bias remains -- and the gap seen in the")
    print("  energy-weighted numbers is time-varying loss (cloud/attitude/MPPT).")

    if not args.no_plot:
        make_plot(df, csv_path, alt_grid, zen_grid, m_grid, r_grid, clean)
    return 0


def make_plot(df, csv_path, alt_grid, zen_grid, m_grid, r_grid, clean):
    import warnings

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import plot_style
    plot_style.apply()   # seaborn whitegrid at cmr10 -- see plot_style.py

    # Empty altitude bins are expected (the aircraft does not dwell at every
    # kilometre); an all-NaN bin is a gap in the line, not an error.
    warnings.filterwarnings("ignore", message="All-NaN slice encountered")

    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    stem = csv_path.name.replace("_solar_efficiency.csv", "")
    # cmr10 has no underscore glyph -- see plot_style.py -- so the log stem is
    # rendered with hyphens rather than silently mangled into apostrophes.
    fig.suptitle("Spectral mismatch + temperature vs the modeled/measured gap"
                 f"\n{stem.replace('_', '-')}", fontsize=12)

    # (a) M vs altitude, several zeniths
    ax = axes[0, 0]
    for z in (0.0, 30.0, 50.0, 70.0):
        ax.plot(alt_grid / 1000.0,
                interp_mismatch(alt_grid, zen_grid, m_grid,
                                alt_grid, np.full_like(alt_grid, z)),
                label=f"zenith {z:.0f}$^\\circ$")
    for z, ls in ((20.0, "--"), (60.0, "-.")):
        ax.plot(alt_grid / 1000.0,
                interp_mismatch(alt_grid, zen_grid, r_grid,
                                alt_grid, np.full_like(alt_grid, z)),
                ls, color="tab:red", lw=1.0,
                label=f"r broadband, z={z:.0f}$^\\circ$")
    ax.axhline(1.0, color="k", ls=":", lw=0.8)
    ax.set_xlabel("altitude [km]"); ax.set_ylabel("factor")
    ax.set_title("(a) M: unusable-IR loss grows with altitude\n"
                 "r: the gray model under-predicts DNI down low")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # (b) the two correction factors along the flight
    ax = axes[0, 1]
    ax.plot(df["time"], df["m_spectral"], lw=0.7, label="M spectral")
    ax.plot(df["time"], df["f_temp"], lw=0.7, label="f temperature")
    ax.plot(df["time"], df["r_broadband"], lw=0.7, label="r broadband")
    ax.plot(df["time"], df["m_spectral"] * df["f_temp"] / df["r_broadband"],
            lw=1.0, color="k", label="net = M*f/r")
    ax.axhline(1.0, color="k", ls=":", lw=0.8)
    ax.set_ylabel("correction factor")
    ax.set_ylim(0.75, 1.45)   # clipped: r dives to ~0.45 at the grazing sun of
                              # the landing, which is real but off-scale here
    ax.set_title("(b) spectrum and temperature partly cancel;\n"
                 "the broadband bias r is what actually moves")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.tick_params(axis="x", rotation=30, labelsize=7)

    # (c) measured/modeled vs altitude, with corrections overlaid
    ax = axes[1, 0]
    a = df["alt_msl_m"] / 1000.0
    P = pd.to_numeric(df["pv_power_actual_w"], errors="coerce")
    eff = pd.to_numeric(df["pre_mppt_efficiency_pct"], errors="coerce")
    ratio = eff / 100.0
    good = ratio.notna() & (ratio > 0) & (ratio < 1.5) & (P > 5)
    ax.scatter(a[good], ratio[good], s=1, alpha=0.12, color="tab:blue",
               label="measured / modeled")
    bins = np.arange(0, 19, 1.0)
    idx = np.digitize(a[good], bins)
    med = [np.nanmedian(ratio[good][idx == i]) if (idx == i).sum() > 20 else np.nan
           for i in range(1, len(bins))]
    ctr = 0.5 * (bins[:-1] + bins[1:])
    ax.plot(ctr, med, "o-", color="tab:blue", label="median")
    mc = (df["m_spectral"] * df["f_temp"] / df["r_broadband"])[good]
    med_c = [np.nanmedian((ratio[good] / mc)[idx == i]) if (idx == i).sum() > 20 else np.nan
             for i in range(1, len(bins))]
    ax.plot(ctr, med_c, "s-", color="tab:red",
            label="median, all three removed")
    ax.axhline(1.0, color="k", ls=":", lw=0.8)
    ax.set_xlabel("altitude [km]"); ax.set_ylabel("ratio")
    ax.set_ylim(0, 1.4)
    ax.set_title("(c) does removing all three flatten the altitude trend?")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # (d) empirical vs modeled M
    ax = axes[1, 1]
    sub = df[clean]
    if len(sub) > 50:
        ax.scatter(sub["alt_msl_m"] / 1000.0, sub["m_empirical"], s=1, alpha=0.12,
                   color="tab:green", label="empirical (from Impp)")
        aa = sub["alt_msl_m"].to_numpy(float) / 1000.0
        idx = np.digitize(aa, bins)
        me = [np.nanmedian(sub["m_empirical"].to_numpy()[idx == i])
              if (idx == i).sum() > 20 else np.nan for i in range(1, len(bins))]
        ax.plot(ctr, me, "o-", color="tab:green", label="empirical median (raw)")
        mf = [np.nanmedian(sub["m_empirical_fixed"].to_numpy()[idx == i])
              if (idx == i).sum() > 20 else np.nan for i in range(1, len(bins))]
        ax.plot(ctr, mf, "^-", color="tab:purple",
                label="empirical median (broadband-fixed)")
        ax.plot(alt_grid / 1000.0,
                interp_mismatch(alt_grid, zen_grid, m_grid, alt_grid,
                                np.full_like(alt_grid, 40.0)),
                color="tab:red", label="SPECTRL2 M (zenith 40$^\\circ$)")
    ax.axhline(1.0, color="k", ls=":", lw=0.8)
    ax.set_xlabel("altitude [km]"); ax.set_ylabel("M")
    ax.set_ylim(0, 1.6)
    ax.set_title("(d) current-derived M vs modeled M\n"
                 "(scale offset = Impp,stc anchor)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    fig.tight_layout()
    out = csv_path.parent / f"{stem}_spectral_mismatch.png"
    fig.savefig(out, dpi=130)
    print(f"\nSaved plot -> {out}")


if __name__ == "__main__":
    sys.exit(main())
