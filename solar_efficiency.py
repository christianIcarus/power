#!/usr/bin/env python3
"""
solar_efficiency.py

Estimate clear-sky solar irradiance along a flight (from a PX4 ULog) using the
aircraft's own GPS position / altitude / time-of-day, then compare it against
the *measured* pre-MPPT solar array output to estimate DC array efficiency
relative to the cell nameplate rating.

Data used from the log (see e.g. 2026_08_10_SN30solarRT.ulg):
  - vehicle_gps_position   -> latitude_deg, longitude_deg, altitude_msl_m,
                              time_utc_usec (absolute UTC time, not boot time)
  - /zeus/mppt_0, /zeus/mppt_1
                           -> pv_voltage_v, pv_current_a, pv_power_w
                              These are the PANEL-SIDE (pre-MPPT) electrical
                              readings -- what the solar array is actually
                              delivering into the MPPT input, before whatever
                              conversion losses the MPPT itself introduces.
  - /zeus/temperature      -> tc_fuselage_outside_temp_c (Tout proxy; see
                              caveat below). Used only when --apply-temp-derate
                              is passed.

Cell reference data (Maxeon Gen 7 datasheet, 546209 Rev C):
  - Cell area       : ~155 cm^2 per cell
  - Cell efficiency : ~25.4% (Pe/Oe bin boundary -> "typical" cell)
  - Power temp coeff: -0.27 %/degC, relative to STC_TEMP_C (25 degC). Applied
    as an optional adjustment on the theoretical power when --apply-temp-derate
    is passed (off by default -- see "Temperature derating" below). Despite
    the name, this is NOT always a loss: since the coefficient is negative,
    Tout above 25 degC reduces theoretical power (a loss) but Tout BELOW
    25 degC increases it (a gain). The code, plot, and summary all handle
    both signs -- don't assume derating only ever shrinks power.

Theoretical/available array power at each instant is layered in three tiers,
each one a further loss on top of the last (see "Encapsulation transmission"
and "Temperature derating" below for the middle and last terms):

    P_bare(t)        = GHI_clearsky(t) [W/m^2] * total_cell_area [m^2] * cell_efficiency
    P_no_temp_derate(t) = P_bare(t) * etfe_transmission * poe_transmission
    P_theoretical(t) = P_no_temp_derate(t) [* temp_derate_factor(t)]

P_bare is the cell nameplate ceiling with nothing between the sun and the
cells -- not physically real for this array, but a useful "no losses at
all" reference. P_no_temp_derate adds the encapsulation stack's light loss
(ETFE cover + POE encapsulant, both always applied, not opt-in) and is what
"efficiency" is actually measured against when --apply-temp-derate is off.
P_theoretical is the final number used everywhere once temp derating is
added on top.

GHI_clearsky is modeled with a Beer-Lambert atmospheric attenuation model,
driven by the *aircraft's* lat/lon/altitude/time at each sample -- not a
fixed ground station -- so it tracks the flight's actual position and
altitude gain. See SEA_LEVEL_TRANSMITTANCE below for why this replaced
pvlib's Ineichen/Perez clear-sky model: Ineichen's altitude correction is an
unbounded linear term fit from ground weather stations, and it produced GHI
above the physical top-of-atmosphere ceiling once extrapolated to this
aircraft's stratospheric cruise altitude (~55-58 kft on the flight this was
first caught on) -- see git history for the full derivation.

IMPORTANT SIMPLIFICATION -- clear sky, not actual sky: GHI_clearsky is a
zero-cloud idealization; it has no mechanism to represent real cloud cover.
So "pre-MPPT efficiency" here bundles true cell/wiring performance together
with any actual clouds present during the flight, which this model cannot
see. On the one flight analyzed so far, the measured-power trace shows a
repeating floor/spike pattern that was initially (incorrectly) suspected to
be aircraft attitude changing angle-of-incidence to the sun -- but measured
power correlates only weakly with roll/pitch (|r| ~ 0.03-0.13 across the
whole flight, checked against /zeus/flight), which rules attitude out as the
dominant driver. Patchy real cloud cover (invisible to a clear-sky model) is
the far more likely explanation. Treat "efficiency" as a nameplate-relative,
clear-sky-relative estimate, not a lab-grade cell-efficiency measurement --
a rigorous version would need real cloud-cover data (satellite irradiance
product or a nearby ground station) rather than pure clear-sky modeling.
Aircraft attitude (/zeus/flight has pitch_deg/roll_deg/yaw_deg) remains a
secondary, empirically minor caveat on top of that.

Encapsulation transmission (always applied, NOT opt-in -- ETFE_TRANSMISSION
and POE_TRANSMISSION, stacked): the light path from sun to cell isn't just
open air -- it passes through an ETFE outer film AND a POE encapsulant
layer bonding the cells underneath it, and each layer's loss compounds:
    encapsulation_transmission = etfe_transmission * poe_transmission
ETFE_TRANSMISSION is a SPECTRALLY-WEIGHTED average, not a flat datasheet
number: a flat "% transmission" would overstate the loss that matters
electrically, since the cell doesn't respond equally to every wavelength
and the sun doesn't deliver equal power at every wavelength either. So
it's transmission(lambda) weighted by [cell EQE(lambda) x AM1.5G solar
spectral irradiance(lambda)], summed over wavelength -- a wavelength the
cell can't use, or where the sun delivers almost no energy there, barely
moves the number even if the film transmits it poorly (or well). Source
data: a manufacturer ETFE light-transmission chart (%T, 200-870 nm) and a
Maxeon spectral-response chart (EQE % + the ASTM G173-03 "global tilt"
AM1.5G reference spectrum it was measured against, both 300-1200 nm). Both
were hand-digitized off chart images (no source data file), and the ETFE
curve's 870-1100 nm tail (past its chart's right edge, but still inside
the cell's response range) was extrapolated flat at its last plotted
value -- override with --etfe-transmission if better data turns up.
POE_TRANSMISSION (0.92) is a single flat figure (user-supplied, 2026-08-25)
-- no spectral chart was given for it, so unlike ETFE it isn't wavelength-
weighted against the cell's response. Override with --poe-transmission if
a spectral curve for it ever shows up and is worth digitizing the same way.

Temperature derating (Tout proxy, opt-in via --apply-temp-derate):
No field literally named "Tout"/"OAT" exists in this log. Considered and
rejected: /zeus/aeroprobe.temp_external_c (dead sentinel 999.0 for the whole
flight despite the rest of that topic being alive), vehicle_air_data.
ambient_temperature and sensor_baro.temperature (clamped ~58-68 degC --
avionics-bay/die self-heating, not outside air), and /zeus/mppt_*.
temperature_c (MPPT converter board temp, not ambient). The field actually
used, /zeus/temperature's tc_fuselage_outside_temp_c, is an exterior-facing
fuselage-skin thermocouple -- it runs hot from solar heating of the skin in
addition to true ambient air temperature, so it is a PROXY, not a shaded
free-air OAT measurement. That said, it is arguably a closer analog to
actual cell operating temperature than a shaded sensor would be, since the
array is exposed to the same solar gain as the skin. Never label it as
generic "OAT" in output -- it is reported as "fuselage skin temp".

Usage:
    python solar_efficiency.py --ulog "C:\\path\\to\\log.ulg"
    python solar_efficiency.py --ulog log.ulg --no-plot
    python solar_efficiency.py --ulog log.ulg --cell-count 72 --cell-efficiency 0.254
    python solar_efficiency.py --ulog log.ulg --apply-temp-derate
    python solar_efficiency.py --ulog log.ulg --etfe-transmission 0.93
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pvlib
from pyulog import ULog
from timezonefinder import TimezoneFinder

# --------------------------------------------------------------------------
# Reference constants (Maxeon Gen 7 datasheet 546209 Rev C)
# --------------------------------------------------------------------------
DEFAULT_ULOG = r"C:\Users\ChristianHammerly\Downloads\2026_08_10_SN30solarRT.ulg"
DEFAULT_CELL_COUNT = 72
DEFAULT_CELL_AREA_CM2 = 155.0        # per cell
DEFAULT_CELL_EFFICIENCY = 0.254      # 25.4%, typical Pe/Oe bin boundary
POWER_TEMP_COEFF_PCT_PER_C = -0.27   # informational, not applied by default
STC_TEMP_C = 25.0
STC_IRRADIANCE_W_M2 = 1000.0

# Beer-Lambert sea-level zenith transmittance: broadband atmospheric
# transmittance looking straight up through the WHOLE atmosphere, at sea
# level, sun directly overhead. Value and model both borrowed from
# Icarus-Matrix/vehicle-simulation's endurance calculator (apollo.yaml,
# rebuild-endurance-calculator branch @ 2026-08-25), which validates this
# exact form against a reference irradiance workbook. Chosen over pvlib's
# Ineichen/Perez clear-sky model specifically because this form is bounded
# by construction: attenuation is transmittance^pressure_ratio(altitude), and
# pressure_ratio only ever approaches 0 as altitude increases, so
# transmittance only ever approaches 1 (no attenuation) -- GHI can get
# arbitrarily close to but can never exceed toa_irradiance * cos(zenith) at
# any altitude. Ineichen's altitude term has no equivalent ceiling.
SEA_LEVEL_TRANSMITTANCE = 0.70

# ETFE array-cover light transmission, spectrally-weighted -- see the
# "Encapsulation transmission" docstring section above for the why. The table
# below is the actual hand-digitized data behind the number, kept here (not
# just the final scalar) so the derivation can be audited/redone if better
# source data shows up. Columns: wavelength [nm], AM1.5G global-tilt solar
# spectral irradiance [W/m^2/nm] (ASTM G173-03), Maxeon cell EQE [fraction],
# ETFE transmission [fraction] (flat-extrapolated past 850 nm -- see above).
_ETFE_SPECTRAL_DATA = [
    # wl_nm, am15g_w_m2_nm, maxeon_eqe, etfe_transmission
    (300, 0.05, 0.65, 0.895), (350, 0.50, 0.80, 0.900), (400, 1.10, 0.90, 0.905),
    (450, 1.70, 0.96, 0.915), (500, 1.85, 0.98, 0.920), (550, 1.75, 0.99, 0.925),
    (600, 1.60, 0.99, 0.930), (650, 1.50, 0.99, 0.935), (700, 1.40, 0.98, 0.940),
    (750, 1.20, 0.97, 0.940), (800, 1.10, 0.95, 0.945), (850, 0.97, 0.92, 0.945),
    (900, 0.85, 0.85, 0.945), (950, 0.70, 0.70, 0.945), (1000, 0.85, 0.45, 0.945),
    (1050, 0.75, 0.20, 0.945), (1100, 0.70, 0.05, 0.945),
]


def _spectrally_weighted_etfe_transmission() -> float:
    """T_eff = sum(T(lambda) * EQE(lambda) * AM1.5G(lambda)) / sum(EQE(lambda) * AM1.5G(lambda))."""
    weight_sum = sum(am15g * eqe for _, am15g, eqe, _ in _ETFE_SPECTRAL_DATA)
    weighted_transmission = sum(am15g * eqe * t for _, am15g, eqe, t in _ETFE_SPECTRAL_DATA)
    return weighted_transmission / weight_sum


DEFAULT_ETFE_TRANSMISSION = round(_spectrally_weighted_etfe_transmission(), 4)  # ~0.93

# POE encapsulant light transmission -- the layer bonding the cells that sits
# UNDER the ETFE cover, so its loss stacks with ETFE_TRANSMISSION rather than
# replacing it (see "Encapsulation transmission" docstring section). Flat
# figure (user-supplied, 2026-08-25), not spectrally-weighted like ETFE
# above -- no spectral transmission chart was provided for it.
DEFAULT_POE_TRANSMISSION = 0.92

GPS_TOPIC = "vehicle_gps_position"
MPPT_TOPICS = ("/zeus/mppt_0", "/zeus/mppt_1")
TEMP_TOPIC = "/zeus/temperature"
TEMP_FIELD = "tc_fuselage_outside_temp_c"          # Tout proxy -- see docstring caveat
TOUT_CLAMP_RANGE_C = (-40.0, 85.0)                 # defensive: generous aviation/electronics envelope


def detect_launch_timezone(lat: float, lon: float) -> str:
    """IANA timezone name for the launch site, from its first GPS fix."""
    tz = TimezoneFinder().timezone_at(lat=lat, lng=lon)
    if tz is None:
        raise RuntimeError(
            f"Could not resolve a timezone for launch site ({lat:.5f}, {lon:.5f}); "
            "pass --tz explicitly."
        )
    return tz


# --------------------------------------------------------------------------
# ULog loading
# --------------------------------------------------------------------------
def _get_dataset(ulog: ULog, name: str):
    matches = [d for d in ulog.data_list if d.name == name]
    return matches[0] if matches else None


def load_gps(ulog: ULog) -> pd.DataFrame:
    """GPS fixes with absolute UTC time, filtered to valid 3D fixes."""
    gps = _get_dataset(ulog, GPS_TOPIC)
    if gps is None:
        raise RuntimeError(f"Topic '{GPS_TOPIC}' not found in log")

    df = pd.DataFrame(
        {
            "lat": gps.data["latitude_deg"],
            "lon": gps.data["longitude_deg"],
            "alt_msl_m": gps.data["altitude_msl_m"],
            "fix_type": gps.data["fix_type"],
            "time_utc_us": gps.data["time_utc_usec"],
        }
    )
    df = df[(df["fix_type"] >= 3) & (df["time_utc_us"] > 0)]
    df["time"] = pd.to_datetime(df["time_utc_us"], unit="us", utc=True)
    df = df.drop(columns=["fix_type", "time_utc_us"])
    df = df.drop_duplicates("time").set_index("time").sort_index()
    return df


def load_mppt(ulog: ULog, sync_tolerance_s: float) -> pd.DataFrame:
    """Combine both MPPT channels' pre-MPPT PV readings onto one timeline.

    /zeus/mppt_* messages carry their own absolute-UTC field (timestamp_us),
    separate from the boot-relative `timestamp` field used for log ordering.
    """
    frames = []
    for i, topic in enumerate(MPPT_TOPICS):
        m = _get_dataset(ulog, topic)
        if m is None:
            print(f"  (note: {topic} not present in this log, skipping)")
            continue
        df = pd.DataFrame(
            {
                "time_utc_us": m.data["timestamp_us"],
                f"pv_voltage_v_{i}": m.data["pv_voltage_v"],
                f"pv_current_a_{i}": m.data["pv_current_a"],
                f"pv_power_w_{i}": m.data["pv_power_w"],
                f"connected_{i}": m.data["connected"],
            }
        )
        df["time"] = pd.to_datetime(df["time_utc_us"], unit="us", utc=True)
        df = df.drop(columns="time_utc_us").drop_duplicates("time").set_index("time").sort_index()
        frames.append(df)

    if not frames:
        raise RuntimeError("No /zeus/mppt_* channels found in log")

    merged = frames[0]
    tol = pd.Timedelta(seconds=sync_tolerance_s)
    for extra in frames[1:]:
        merged = pd.merge_asof(
            merged.reset_index(), extra.reset_index(),
            on="time", direction="nearest", tolerance=tol,
        ).set_index("time")
    return merged


def load_temperature(ulog: ULog) -> pd.DataFrame | None:
    """Outside-air-temperature (Tout) proxy, used only for --apply-temp-derate.

    /zeus/temperature carries its own absolute-UTC timestamp_us field, the
    same convention as /zeus/mppt_*. tc_fuselage_outside_temp_c is an
    exterior-facing skin thermocouple -- NOT a shaded free-air OAT sensor
    (see module docstring). Other candidates considered and rejected:
    /zeus/aeroprobe.temp_external_c (dead 999.0 sentinel for the whole
    flight), vehicle_air_data.ambient_temperature / sensor_baro.temperature
    (avionics self-heating, non-physical ~58-68 degC), /zeus/mppt_*.
    temperature_c (MPPT converter board temp, not ambient).
    """
    t = _get_dataset(ulog, TEMP_TOPIC)
    if t is None:
        print(f"  (note: {TEMP_TOPIC} not present in this log - temp derating disabled)")
        return None
    df = pd.DataFrame(
        {
            "time_utc_us": t.data["timestamp_us"],
            "tout_c": t.data[TEMP_FIELD],
        }
    )
    df["time"] = pd.to_datetime(df["time_utc_us"], unit="us", utc=True)
    return df.drop(columns="time_utc_us").drop_duplicates("time").set_index("time").sort_index()


# --------------------------------------------------------------------------
# Solar geometry / clear-sky irradiance
# --------------------------------------------------------------------------
def isa_pressure_ratio(alt_m: np.ndarray) -> np.ndarray:
    """ISA static pressure as a fraction of sea level, piecewise to 32 km.

    p/p0 = (rho/rho0) * (T/T0) via the ideal gas law. Vectorized port of
    Icarus-Matrix/vehicle-simulation's air_density()/air_temperature_K()
    (endurance.py); 32 km covers any altitude this aircraft flies (this
    flight peaked at ~17.6 km).
    """
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


def compute_clearsky_irradiance(times: pd.DatetimeIndex, lat: np.ndarray,
                                 lon: np.ndarray, alt_m: np.ndarray) -> pd.DataFrame:
    """Per-sample clear-sky GHI using the aircraft's own position/time/altitude
    (Beer-Lambert attenuation -- see SEA_LEVEL_TRANSMITTANCE for why, not
    pvlib's Ineichen/Perez model).

    dni_w_m2/dhi_w_m2 are a simplified all-beam/no-diffuse split for the CSV
    export only -- this model doesn't compute diffuse sky radiation
    separately (it was a ~1% contribution even under the old Ineichen model
    at this flight's altitudes). ghi_w_m2 is the only one of the three that
    feeds the theoretical-power calc, summary, and plot.
    """
    solpos = pvlib.solarposition.get_solarposition(times, lat, lon, altitude=alt_m)
    cos_zenith = np.cos(np.radians(solpos["apparent_zenith"].values))
    sun_up = cos_zenith > 0.0

    dni_extra = np.asarray(pvlib.irradiance.get_extra_radiation(times), dtype=float)
    tau = SEA_LEVEL_TRANSMITTANCE ** isa_pressure_ratio(alt_m)

    # Plane-parallel airmass (1/cos z). Only meaningful where the sun is up;
    # elsewhere cos_zenith<=0 would blow this up for no reason since ghi/dni
    # are forced to 0 below regardless.
    safe_cos_zenith = np.where(sun_up, cos_zenith, 1.0)
    direct_normal_w_m2 = dni_extra * tau ** (1.0 / safe_cos_zenith)

    ghi = np.where(sun_up, direct_normal_w_m2 * cos_zenith, 0.0)
    dni = np.where(sun_up, direct_normal_w_m2, 0.0)
    dhi = np.zeros_like(ghi)  # not modeled -- see docstring above

    out = pd.DataFrame(
        {
            "sun_elevation_deg": solpos["apparent_elevation"].values,
            "sun_azimuth_deg": solpos["azimuth"].values,
            "ghi_w_m2": ghi,
            "dni_w_m2": dni,
            "dhi_w_m2": dhi,
        },
        index=times,
    )
    return out


# --------------------------------------------------------------------------
# Main analysis
# --------------------------------------------------------------------------
def analyze(args: argparse.Namespace) -> pd.DataFrame:
    ulog_path = Path(args.ulog)
    if not ulog_path.exists():
        raise FileNotFoundError(f"ULog file not found: {ulog_path}")

    print(f"Loading {ulog_path.name} ...")
    topics = [GPS_TOPIC, *MPPT_TOPICS, TEMP_TOPIC]
    ulog = ULog(str(ulog_path), message_name_filter_list=topics)

    gps_df = load_gps(ulog)
    mppt_df = load_mppt(ulog, args.mppt_sync_tolerance_s)
    print(f"  GPS fixes: {len(gps_df)}   MPPT samples: {len(mppt_df)}")

    # Attach nearest GPS fix (position/altitude) to each MPPT sample.
    tol = pd.Timedelta(seconds=args.gps_tolerance_s)
    merged = pd.merge_asof(
        mppt_df.reset_index(), gps_df.reset_index(),
        on="time", direction="nearest", tolerance=tol,
    ).set_index("time")
    merged = merged.dropna(subset=["lat", "lon", "alt_msl_m"])
    if merged.empty:
        raise RuntimeError("No MPPT samples could be matched to a GPS fix - "
                            "try increasing --gps-tolerance-s")

    # Attach the Tout proxy (fuselage skin temp) alongside GPS, on the same
    # MPPT-anchored timeline. Only used downstream if --apply-temp-derate.
    temp_df = load_temperature(ulog)
    if temp_df is not None:
        tol_temp = pd.Timedelta(seconds=args.temp_tolerance_s)
        merged = pd.merge_asof(
            merged.reset_index(), temp_df.reset_index(),
            on="time", direction="nearest", tolerance=tol_temp,
        ).set_index("time")
    else:
        merged["tout_c"] = np.nan

    # Sanity-filter obviously-bad raw readings before they can reach the
    # printed summary, CSV, or plot -- not just the derate math. Seen in the
    # wild: a hard sentinel of exactly -273.15 degC (0 Kelvin), presumably an
    # uninitialized/fault value from the sensor or its logging path. Anything
    # outside the generous TOUT_CLAMP_RANGE_C envelope is treated as missing
    # (not as a real extreme reading) and forward-filled below like any other
    # gap. Runs unconditionally, regardless of --apply-temp-derate, since a
    # sentinel like this corrupts what's *displayed*, not just the derating.
    raw_tout = merged["tout_c"]
    is_bad_tout = raw_tout.notna() & (
        (raw_tout < TOUT_CLAMP_RANGE_C[0]) | (raw_tout > TOUT_CLAMP_RANGE_C[1])
    )
    n_bad_tout = int(is_bad_tout.sum())
    if n_bad_tout:
        print(f"  WARNING: {n_bad_tout} Tout readings ({100 * n_bad_tout / len(merged):.1f}%) fell "
              f"outside the plausible {TOUT_CLAMP_RANGE_C} degC envelope (e.g. a dead/fault "
              f"sentinel) - treated as missing and forward-filled, not plotted/derated as-read.")
        merged.loc[is_bad_tout, "tout_c"] = np.nan

    # Tout is slow-moving (whole-flight range is typically only a few degC),
    # so forward-filling brief match gaps introduces negligible error. Rows
    # with no match at all (before the first sample, or topic missing) stay
    # NaN here and get a no-op derate factor (1.0) below.
    n_missing_before_ffill = merged["tout_c"].isna().sum()
    merged["tout_c"] = merged["tout_c"].ffill()
    n_missing = merged["tout_c"].isna().sum()
    if n_missing:
        print(f"  WARNING: {n_missing} samples ({100 * n_missing / len(merged):.1f}%) had no "
              f"Tout data (before first temperature sample, or topic missing) - "
              f"temp derate factor defaults to 1.0 (no derating) for these rows.")
    elif n_missing_before_ffill:
        print(f"  (note: {n_missing_before_ffill} Tout samples forward-filled across small gaps)")

    # Total measured pre-MPPT array power (sum across whatever channels exist).
    power_cols = [c for c in merged.columns if c.startswith("pv_power_w_")]
    merged["pv_power_actual_w"] = merged[power_cols].sum(axis=1, skipna=True)

    # Flag channels that contributed essentially nothing (likely unused/disconnected).
    total_actual = merged["pv_power_actual_w"].clip(lower=0).sum()
    for col in power_cols:
        channel_total = merged[col].clip(lower=0).sum()
        share = channel_total / total_actual if total_actual > 0 else 0.0
        if share < 0.01:
            print(f"  WARNING: {col} contributed only {share * 100:.2f}% of total "
                  f"measured power - likely disconnected/unused for this flight.")

    # Clear-sky irradiance at the aircraft's own position/time.
    print("Computing solar position (pvlib) + clear-sky irradiance (Beer-Lambert) ...")
    irr = compute_clearsky_irradiance(
        merged.index, merged["lat"].values, merged["lon"].values, merged["alt_msl_m"].values
    )
    merged = merged.join(irr)

    # Optional temperature derating: factor = 1 + coeff%/degC * (Tout - STC_TEMP_C).
    # coeff is negative, so Tout ABOVE STC_TEMP_C (25 degC) is a LOSS
    # (factor < 1) but Tout BELOW STC_TEMP_C is a GAIN (factor > 1) -- this
    # is not just a loss term, and downstream code/plots must not assume it
    # always reduces power. Sentinel/out-of-envelope raw readings are already
    # scrubbed to NaN-then-ffilled above (so they can't reach here at all);
    # the clip below is just a final defensive no-op in case that upstream
    # sanitation is ever skipped, and the resulting factor is floored at 0
    # (power can shrink toward zero in a pathological case, but never go
    # negative) as a last backstop.
    if args.apply_temp_derate:
        tout_clamped = merged["tout_c"].clip(*TOUT_CLAMP_RANGE_C)
        merged["temp_derate_factor"] = (
            1.0 + (POWER_TEMP_COEFF_PCT_PER_C / 100.0) * (tout_clamped - STC_TEMP_C)
        )
        merged["temp_derate_factor"] = merged["temp_derate_factor"].fillna(1.0).clip(lower=0.0)
    else:
        merged["temp_derate_factor"] = 1.0

    # Theoretical array output, three tiers (see module docstring for the
    # full breakdown): bare-cell ceiling -> encapsulation loss (ETFE cover x
    # POE encapsulant, both always applied, stacked -- see "Encapsulation
    # transmission") -> optional temp derate. Bare is kept around
    # unconditionally as a reference curve; the encapsulated "nominal" is
    # only kept as its own column when temp derating is on, so the plot can
    # show all three and make the derate's effect visible on top of the
    # (already-applied) encapsulation loss -- when derating is off, nominal
    # IS the final theoretical number.
    total_area_m2 = args.cell_count * args.cell_area_cm2 / 1e4
    bare_w = merged["ghi_w_m2"] * total_area_m2 * args.cell_efficiency
    merged["pv_power_theoretical_bare_w"] = bare_w
    encapsulation_transmission = args.etfe_transmission * args.poe_transmission
    nominal_w = bare_w * encapsulation_transmission
    if args.apply_temp_derate:
        merged["pv_power_theoretical_nominal_w"] = nominal_w
        merged["pv_power_theoretical_w"] = nominal_w * merged["temp_derate_factor"]
    else:
        merged["pv_power_theoretical_w"] = nominal_w

    valid = (merged["sun_elevation_deg"] > 0) & (merged["pv_power_theoretical_w"] > 1.0)
    merged["pre_mppt_efficiency_pct"] = np.nan
    merged.loc[valid, "pre_mppt_efficiency_pct"] = (
        100.0 * merged.loc[valid, "pv_power_actual_w"] / merged.loc[valid, "pv_power_theoretical_w"]
    )

    merged.attrs["total_area_m2"] = total_area_m2
    merged.attrs["valid_mask"] = valid
    return merged


def print_summary(df: pd.DataFrame, args: argparse.Namespace, tz: str) -> None:
    valid = df.attrs["valid_mask"]
    daylight = df[valid]

    dt_hours = df.index.to_series().diff().dt.total_seconds().fillna(0.0) / 3600.0
    dt_hours = dt_hours.clip(upper=10.0 / 3600.0)  # ignore >10s gaps in energy integration
    energy_actual_wh = (df["pv_power_actual_w"].clip(lower=0) * dt_hours).sum()
    energy_theoretical_wh = (df["pv_power_theoretical_w"].fillna(0) * dt_hours).sum()

    rated_stc_w = args.cell_count * args.cell_area_cm2 / 1e4 * args.cell_efficiency * STC_IRRADIANCE_W_M2

    local_start = df.index[0].tz_convert(tz)
    local_end = df.index[-1].tz_convert(tz)

    print("\n" + "=" * 70)
    print("SOLAR ARRAY SUMMARY  (pre-MPPT, panel-side measurement)")
    print("=" * 70)
    print(f"Launch-site timezone     : {tz}")
    print(f"Flight window (local)    : {local_start.strftime('%Y-%m-%d %H:%M:%S %Z')} "
          f"-> {local_end.strftime('%H:%M:%S %Z')}")
    print(f"Flight window (UTC)      : {df.index[0].strftime('%Y-%m-%d %H:%M:%S %Z')} "
          f"-> {df.index[-1].strftime('%H:%M:%S %Z')}")
    print(f"Flight duration          : {(df.index[-1] - df.index[0])}")
    print(f"Samples analyzed         : {len(df)}  ({len(daylight)} with sun above horizon)")
    print(f"Cells / area / eff.      : {args.cell_count} x {args.cell_area_cm2:.0f} cm^2 "
          f"@ {args.cell_efficiency * 100:.1f}%  -> {df.attrs['total_area_m2']:.3f} m^2 total")
    print(f"Encapsulation transmission: ETFE {args.etfe_transmission * 100:.1f}% x POE "
          f"{args.poe_transmission * 100:.1f}% = {args.etfe_transmission * args.poe_transmission * 100:.1f}%  "
          f"(always applied -- see docstring)")
    print(f"STC-rated array power    : {rated_stc_w:.1f} W  (at 1000 W/m^2, 25 degC, bare cells)")
    print(f"Peak measured PV power   : {df['pv_power_actual_w'].max():.1f} W")
    print(f"Peak modeled GHI         : {df['ghi_w_m2'].max():.1f} W/m^2")
    print(f"Peak sun elevation       : {df['sun_elevation_deg'].max():.1f} deg")
    if "tout_c" in df.columns and df["tout_c"].notna().any():
        print(f"Fuselage skin temp (Tout proxy): mean {df['tout_c'].mean():.1f} degC  "
              f"range [{df['tout_c'].min():.1f}, {df['tout_c'].max():.1f}] degC")
        if args.apply_temp_derate:
            # Positive change_pct = LOSS (Tout above STC), negative = GAIN
            # (Tout below STC) -- coeff is negative, so don't assume a loss.
            change_pct = 100.0 * (1.0 - df["temp_derate_factor"].mean())
            direction = "loss" if change_pct >= 0 else "gain"
            print(f"Temp derate effect       : mean {abs(change_pct):.1f}% {direction}  "
                  f"(STC {STC_TEMP_C:.0f} degC, coeff {POWER_TEMP_COEFF_PCT_PER_C}%/degC)")
    print(f"Energy delivered (meas.) : {energy_actual_wh:.1f} Wh")
    print(f"Energy available (model.): {energy_theoretical_wh:.1f} Wh")
    if len(daylight):
        eff = daylight["pre_mppt_efficiency_pct"]
        print(f"Pre-MPPT efficiency      : mean {eff.mean():.1f}%  median {eff.median():.1f}%  "
              f"p10-p90 [{eff.quantile(0.1):.1f}%, {eff.quantile(0.9):.1f}%]")
    if energy_theoretical_wh > 0:
        print(f"Energy-weighted efficiency: {100 * energy_actual_wh / energy_theoretical_wh:.1f}%")
    print("=" * 70)
    print("NOTE: efficiency is measured against modeled CLEAR-SKY GHI, which has")
    print("no mechanism to represent real cloud cover -- any actual clouds during")
    print("the flight will show up here as 'lower efficiency', not as reduced")
    print("irradiance. Aircraft attitude is a secondary, empirically minor caveat")
    print("(measured power correlates only weakly with roll/pitch, |r| ~ 0.03-0.13).")


def make_plot(df: pd.DataFrame, out_path: Path, tz: str) -> None:
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    local_index = df.index.tz_convert(tz)
    df = df.set_axis(local_index)

    has_tout = "tout_c" in df.columns and df["tout_c"].notna().any()
    fig, axes = plt.subplots(3, 1, figsize=(12, 10.5), sharex=True)

    def legend_outside(ax, *extra_axes):
        """Combine handles from ax + any twinx axes into one legend box
        placed clear of the plotted data, to the right of the panel."""
        handles, labels = ax.get_legend_handles_labels()
        for extra in extra_axes:
            h, l = extra.get_legend_handles_labels()
            handles += h
            labels += l
        ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.06, 1.0), borderaxespad=0.0)

    ax = axes[0]
    ax.plot(df.index, df["ghi_w_m2"], color="tab:orange", label="Modeled Clear-Sky GHI @ Lat/Lon/Alt")
    ax2 = ax.twinx()
    ax2.plot(df.index, df["sun_elevation_deg"], color="tab:gray", alpha=0.6, label="Sun Elevation")
    ax.set_ylabel("GHI: Global Horizontal Irradiance (W/m^2)")
    ax2.set_ylabel("Sun elevation (deg)")
    ax.set_title("Modeled Clear-Sky Irradiance Along the Flight Track", fontweight="bold")
    legend_outside(ax, ax2)

    has_derate = "pv_power_theoretical_nominal_w" in df.columns

    ax = axes[1]
    ax.plot(df.index, df["pv_power_actual_w"], color="tab:blue", label="Measured Pre-MPPT Power")
    # Bare-cell ceiling (no ETFE, no temp derate) is always drawn as a thin
    # reference line -- not physically real for this array, but it's the
    # only way to see the ETFE cover's loss on the plot when temp derating
    # is off (in that case "nominal"/"theoretical" below are the same curve).
    ax.plot(df.index, df["pv_power_theoretical_bare_w"], color="tab:gray", linestyle=":",
            alpha=0.6, label="Theoretical, No Temp Derate, Bare Cells")
    if has_derate:
        ax.plot(df.index, df["pv_power_theoretical_nominal_w"], color="tab:green", linestyle="--",
                alpha=0.6, label="Theoretical, No Temp Derate, Encapsulated")
        ax.plot(df.index, df["pv_power_theoretical_w"], color="tab:olive", linestyle="-.",
                label="Theoretical, Temp-Derated")
        # Tout above STC (25 degC) is a LOSS (derated < nominal); Tout below
        # STC is a GAIN (derated > nominal) since the coefficient is
        # negative -- shade/label each region distinctly rather than
        # assuming it's always a loss. Only label a region if it actually
        # occurs, so the legend doesn't show an entry with no matching area.
        derated = df["pv_power_theoretical_w"]
        nominal = df["pv_power_theoretical_nominal_w"]
        is_loss = derated < nominal
        is_gain = derated > nominal
        ax.fill_between(df.index, derated, nominal, where=is_loss, interpolate=True,
                         color="tab:red", alpha=0.15,
                         label="Temp Derate Loss" if is_loss.any() else None)
        ax.fill_between(df.index, derated, nominal, where=is_gain, interpolate=True,
                         color="tab:cyan", alpha=0.15,
                         label="Temp Derate Gain" if is_gain.any() else None)
    else:
        # pv_power_theoretical_w already has the ETFE loss baked in (see
        # analyze()), so it IS the "No Temp Derate, Encapsulated" tier here --
        # there's no separate nominal column to plot on top of the bare-cell line.
        ax.plot(df.index, df["pv_power_theoretical_w"], color="tab:green", linestyle="--",
                label="Theoretical, No Temp Derate, Encapsulated")
    ax.set_ylabel("Power (W)")
    ax.set_title("Measured vs. Theoretical Array Power", fontweight="bold")
    legend_outside(ax)

    ax = axes[2]
    ax.plot(df.index, df["alt_msl_m"], color="tab:purple", label="Altitude (MSL)")
    ax.set_ylabel("Altitude MSL (m)")
    if has_tout:
        ax2 = ax.twinx()
        ax2.plot(df.index, df["tout_c"], color="tab:brown", label="Fuselage Skin Temp (Tout Proxy)")
        ax2.axhline(STC_TEMP_C, color="black", linewidth=0.8, linestyle=":", label=f"STC ({STC_TEMP_C:.0f} degC)")
        ax2.set_ylabel("Temp (degC)")
        ax.set_title("Environmentals", fontweight="bold")
        legend_outside(ax, ax2)
    else:
        ax.set_title("Environmentals", fontweight="bold")
        legend_outside(ax)

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=df.index.tz))
    axes[-1].set_xlabel(f"Local time, {tz} ({df.index[0].date()})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot -> {out_path}")


def open_in_vscode(path: Path) -> None:
    """Open the generated plot in VS Code, reusing the existing window.

    `-r`/`--reuse-window` targets the already-open VS Code window instead of
    spawning a new one; re-running the script therefore just refreshes the
    same image tab rather than piling up new windows each time. Best-effort:
    if the `code` CLI isn't on PATH (e.g. not installed, or "Shell Command:
    Install 'code' command" was never run) or the call otherwise fails, this
    prints a note instead of failing the whole analysis run.
    """
    # Resolve to the actual code(.CMD) path rather than passing the bare
    # "code" string to subprocess: on Windows the launcher is a .CMD shim,
    # and CreateProcess (which subprocess uses without shell=True) won't
    # apply PATHEXT resolution to find it the way a shell would -- so the
    # unresolved name raises "WinError 2: The system cannot find the file
    # specified" even though shutil.which() locates it just fine.
    code_cmd = shutil.which("code")
    if code_cmd is None:
        print("  (note: 'code' CLI not found on PATH - skipping VS Code auto-open. "
              "Run \"Shell Command: Install 'code' command in PATH\" from VS Code's "
              "command palette to enable this.)")
        return
    try:
        subprocess.run([code_cmd, "-r", str(path)], check=True)
    except (subprocess.CalledProcessError, OSError) as exc:
        print(f"  (note: could not auto-open {path.name} in VS Code: {exc})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ulog", default=DEFAULT_ULOG, help="Path to the PX4 .ulg flight log")
    parser.add_argument("--cell-count", type=int, default=DEFAULT_CELL_COUNT)
    parser.add_argument("--cell-area-cm2", type=float, default=DEFAULT_CELL_AREA_CM2)
    parser.add_argument("--cell-efficiency", type=float, default=DEFAULT_CELL_EFFICIENCY,
                         help="Fractional cell efficiency, e.g. 0.254 for 25.4%%")
    parser.add_argument("--etfe-transmission", type=float, default=DEFAULT_ETFE_TRANSMISSION,
                         help="Fractional light transmission through the ETFE array cover, "
                              "spectrally-weighted by cell EQE x AM1.5G spectrum (see "
                              "docstring). Always applied, not opt-in like --apply-temp-derate. "
                              f"Default {DEFAULT_ETFE_TRANSMISSION:.2f}.")
    parser.add_argument("--poe-transmission", type=float, default=DEFAULT_POE_TRANSMISSION,
                         help="Fractional light transmission through the POE encapsulant "
                              "(stacks with --etfe-transmission -- see docstring). Flat figure, "
                              f"not spectrally-weighted. Always applied. Default {DEFAULT_POE_TRANSMISSION:.2f}.")
    parser.add_argument("--gps-tolerance-s", type=float, default=2.0,
                         help="Max time gap allowed when matching a GPS fix to an MPPT sample")
    parser.add_argument("--mppt-sync-tolerance-s", type=float, default=0.5,
                         help="Max time gap allowed when merging the two MPPT channels together")
    parser.add_argument("--temp-tolerance-s", type=float, default=1.0,
                         help="Max time gap allowed when matching a Tout (fuselage skin temp) "
                              "sample to an MPPT sample")
    parser.add_argument("--apply-temp-derate", action="store_true",
                         help="Derate theoretical power using the datasheet's power temp "
                              "coefficient and the Tout proxy (fuselage skin temp). Off by "
                              "default -- see docstring for why this is opt-in.")
    parser.add_argument("--output-dir", default=None,
                         help="Where to write the CSV/plot (default: alongside the ulog file)")
    parser.add_argument("--no-plot", action="store_true", help="Skip generating the PNG plot")
    parser.add_argument("--no-open", action="store_true",
                         help="Don't auto-open the generated plot in VS Code after saving it "
                              "(on by default; requires the 'code' CLI on PATH)")
    parser.add_argument("--tz", default=None,
                         help="Timezone for the plot's time axis (IANA name, e.g. America/Los_Angeles). "
                              "Default: auto-detected from the launch-site GPS fix. "
                              "Log data is otherwise kept/exported in UTC.")
    args = parser.parse_args()

    df = analyze(args)

    tz = args.tz or detect_launch_timezone(df["lat"].iloc[0], df["lon"].iloc[0])
    print_summary(df, args, tz)

    ulog_path = Path(args.ulog)
    out_dir = Path(args.output_dir) if args.output_dir else ulog_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = ulog_path.stem

    csv_path = out_dir / f"{stem}_solar_efficiency.csv"
    export_cols = [
        "lat", "lon", "alt_msl_m", "sun_elevation_deg", "sun_azimuth_deg",
        "ghi_w_m2", "dni_w_m2", "dhi_w_m2",
        *(["tout_c"] if "tout_c" in df.columns else []),
        *[c for c in df.columns if c.startswith("pv_")],
        *(["temp_derate_factor"] if "temp_derate_factor" in df.columns else []),
        "pre_mppt_efficiency_pct",
    ]
    df[export_cols].to_csv(csv_path)
    print(f"Saved data -> {csv_path}")

    if not args.no_plot:
        plot_path = out_dir / f"{stem}_solar_efficiency.png"
        make_plot(df, plot_path, tz)
        if not args.no_open:
            open_in_vscode(plot_path)


if __name__ == "__main__":
    main()
