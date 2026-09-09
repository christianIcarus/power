#!/usr/bin/env python3
"""
array_power.py

Estimate solar array output power vs time from FLIGHT PATH DATA ALONE.

The only log input is the aircraft's position track -- latitude, longitude,
altitude, and absolute UTC time. Everything else (sun position, irradiance,
panel pointing, ambient temperature, cell temperature, and finally power) is
modeled from that track. No MPPT electrical telemetry, no attitude sensor, no
temperature sensor is read or required.

    vehicle_gps_position -> latitude_deg, longitude_deg, altitude_msl_m,
                            time_utc_usec (absolute UTC, not boot time)

This is deliberately a FORWARD model: it predicts what the array should
produce, and never compares against what it actually produced. If you want
measured-vs-modeled efficiency, that is a different script.

Pipeline
--------
    track (lat/lon/alt/t)
      -> heading            (great-circle bearing between consecutive fixes)
      -> sun position       (pvlib, at the aircraft's own position/altitude)
      -> clear-sky DNI      (Beer-Lambert, altitude-aware)
      -> POA irradiance     (DNI projected onto each string's tilted normal,
                             rotated by heading; wings-level assumption)
      -> ambient temp       (ISA standard atmosphere, from altitude)
      -> cell temp          (NOCT model, from ambient + POA)
      -> array power        (area x efficiency x encapsulation x temp derate)

Per string, then summed. The two MPPT strings are mounted at slightly
different angles and hold different cell counts, so each gets its own POA,
its own cell temperature, and its own power contribution:

    P_string(t) = POA_string(t) [W/m^2]
                  * cell_count_string * cell_area [m^2]
                  * cell_efficiency
                  * etfe_transmission * poe_transmission
                  * temp_derate_factor(T_cell_string(t))
    P_array(t)  = P_string0(t) + P_string1(t)

Heading and panel pointing (wings-level assumption)
---------------------------------------------------
Each string's surface normal is tilted ~8-15 deg off vertical by the wing's
own CAD geometry (see PANEL_NORMAL_BODY_STRING_0/_1). A tilted panel's
pointing depends on which way the aircraft is FACING, so heading matters even
though bank and pitch are not modeled here: as the aircraft turns, that fixed
tilt sweeps around the compass relative to the sun.

Heading is derived from the track itself -- the great-circle initial bearing
between consecutive GPS fixes -- so no attitude sensor is needed. It is
therefore COURSE OVER GROUND, not true heading: in a crosswind the aircraft
crabs, and the two differ by the drift angle. Bearings are smoothed by
circular (sin/cos) averaging over --heading-smoothing-s, and near-stationary
samples (below MIN_GROUND_SPEED_M_S, where bearing is numerically meaningless)
inherit the last good value.

Roll and pitch are assumed ZERO. This is the deliberate simplification that
keeps the model to flight-path data only. It is wrong during banked turns --
a station-keeping aircraft in a continuous orbit is banked essentially all the
time, and bank can swing instantaneous POA by tens of percent either way. The
error is largely symmetric over a full orbit (banking toward the sun on one
side, away on the other), so ENERGY over many orbits is far better estimated
than any INSTANTANEOUS power value. Treat the power trace as an orbit-average
envelope, not a per-second truth. --assume-horizontal drops the tilt model
entirely and uses flat-plate GHI instead.

IMPORTANT SIMPLIFICATION -- clear sky, not actual sky
------------------------------------------------------
The irradiance model is a zero-cloud idealization. It has no mechanism to
represent real cloud cover, haze, or aerosol loading, and will happily report
full sun through an overcast. Estimated power is therefore an upper envelope
under clear conditions, not a forecast of a particular day's actual output.
A rigorous version would need a satellite irradiance product or ground-station
data.

GHI/DNI use a Beer-Lambert atmospheric attenuation model driven by the
aircraft's own lat/lon/altitude/time at each sample. See
SEA_LEVEL_TRANSMITTANCE for why this replaced pvlib's Ineichen/Perez
clear-sky model: Ineichen's altitude correction is an unbounded linear term
fit from ground weather stations, and it produced GHI above the physical
top-of-atmosphere ceiling once extrapolated to this aircraft's stratospheric
cruise altitude (~55-58 kft).

IMPORTANT SIMPLIFICATION -- ISA ambient temperature
----------------------------------------------------
Ambient air temperature comes from the ISA standard atmosphere, which is a
fixed GLOBAL-AVERAGE profile: it is a function of altitude and nothing else.
It does not know latitude, season, or weather. Real atmospheric temperature
at a given altitude varies substantially with all three -- the tropopause sits
near 17 km at the equator but only 9-10 km near the poles, and real deviations
from ISA at cruise altitude are commonly +/-10-20 degC. Because this aircraft
cruises right around ISA's own layer boundaries (11 km and 20 km), it can be
physically in the stratosphere while ISA still models it as tropopause, or
vice versa.

The consequence for this script is modest but real: cell temperature error
propagates through POWER_TEMP_COEFF_PCT_PER_C (-0.27 %/degC), so a 15 degC
ambient error is about 4% of power. That is much smaller than the clear-sky
assumption's error budget, which is why ISA is considered good enough here.
Pass --ambient-offset-c to shift the whole ISA profile if you have reason to
believe the flight ran warm or cold relative to standard.

Cell temperature (NOCT model)
------------------------------
    T_cell = T_ambient + (POA / 800 W/m^2) * (NOCT - 20 degC)

NOCT (Nominal Operating Cell Temperature) is defined at 800 W/m^2, 20 degC
ambient, and 1 m/s wind over an open-rack module. NONE of those conditions
describe this array: it is bonded to a wing skin moving at flight speed, so
convective cooling is far stronger than the 1 m/s the standard assumes, and
the effective NOCT is correspondingly LOWER than a datasheet figure. Thin air
at altitude cuts the other way (less mass flow per unit speed). The default
DEFAULT_NOCT_C is the conventional 45 degC placeholder, NOT a measured or
derived value for this airframe -- override it with --noct-c once real cell
temperature data exists to anchor it. See git history on the `main` branch for
a Vmpp-derived cell-temperature estimator that could provide that anchor.

Cell reference data (Maxeon Gen 7 datasheet, 546209 Rev C)
-----------------------------------------------------------
  - Cell area       : ~155 cm^2 per cell
  - Cell efficiency : ~25.4% (Pe/Oe bin boundary -> "typical" cell)
  - Power temp coeff: -0.27 %/degC relative to STC_TEMP_C (25 degC). Since the
    coefficient is negative, a cell BELOW 25 degC produces MORE than nameplate,
    not less -- and this aircraft's cells spend most of the flight far below
    25 degC. The derate factor here is routinely GREATER than 1. Do not assume
    it only ever shrinks power.

Track input
-----------
Either a flown PX4 log (--ulog) or a planned/hypothetical track CSV
(--track-csv, columns: time, lat, lon, alt_msl_m, and optionally heading_deg).
The physics is identical for both; only where the track comes from differs.

Usage:
    python array_power.py --ulog "C:\\path\\to\\log.ulg"
    python array_power.py --track-csv tracks/rogers_city_5kft_sep29.csv
    python array_power.py --ulog log.ulg --noct-c 35
    python array_power.py --ulog log.ulg --ambient-offset-c -8
    python array_power.py --ulog log.ulg --assume-horizontal
    python array_power.py --ulog log.ulg --no-plot
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
# Array reference constants (Maxeon Gen 7 datasheet 546209 Rev C)
# --------------------------------------------------------------------------
DEFAULT_ULOG = r"C:\Users\ChristianHammerly\Downloads\2026_08_10_SN30solarRT.ulg"
DEFAULT_CELL_COUNT = 72
DEFAULT_STRING1_CELL_COUNT = 58       # of the 72 total -- string 1 is NOT half the array (user, 2026-08-25);
                                       # string 0 gets the remaining 14
DEFAULT_CELL_AREA_CM2 = 155.0         # per cell
DEFAULT_CELL_EFFICIENCY = 0.254       # 25.4%, typical Pe/Oe bin boundary
POWER_TEMP_COEFF_PCT_PER_C = -0.27    # negative: cold cells OVERperform nameplate
STC_TEMP_C = 25.0
STC_IRRADIANCE_W_M2 = 1000.0

# NOCT model constants -- see "Cell temperature" in the module docstring for
# why 45 degC is a placeholder rather than a figure derived for this airframe.
DEFAULT_NOCT_C = 45.0
NOCT_IRRADIANCE_W_M2 = 800.0          # NOCT test condition, fixed by definition
NOCT_AMBIENT_C = 20.0                 # NOCT test condition, fixed by definition

# Beer-Lambert sea-level zenith transmittance: broadband atmospheric
# transmittance looking straight up through the WHOLE atmosphere, at sea
# level, sun directly overhead. Value and model both borrowed from
# Icarus-Matrix/vehicle-simulation's endurance calculator (apollo.yaml,
# rebuild-endurance-calculator branch @ 2026-08-25), which validates this
# exact form against a reference irradiance workbook. Chosen over pvlib's
# Ineichen/Perez clear-sky model specifically because this form is bounded by
# construction: attenuation is transmittance^pressure_ratio(altitude), and
# pressure_ratio only ever approaches 0 as altitude increases, so
# transmittance only ever approaches 1 (no attenuation) -- GHI can get
# arbitrarily close to but never exceed toa_irradiance * cos(zenith) at any
# altitude. Ineichen's altitude term has no equivalent ceiling.
SEA_LEVEL_TRANSMITTANCE = 0.70

# ETFE array-cover light transmission, spectrally-weighted -- a flat "%
# transmission" would overstate the loss that matters electrically, since the
# cell doesn't respond equally to every wavelength and the sun doesn't deliver
# equal power at every wavelength either. So it's transmission(lambda)
# weighted by [cell EQE(lambda) x AM1.5G solar spectral irradiance(lambda)],
# summed over wavelength. Source data: a manufacturer ETFE light-transmission
# chart (%T, 200-870 nm) and a Maxeon spectral-response chart (EQE % + the
# ASTM G173-03 "global tilt" AM1.5G reference spectrum it was measured
# against, both 300-1200 nm). Both were hand-digitized off chart images (no
# source data file), and the ETFE curve's 870-1100 nm tail (past its chart's
# right edge, but still inside the cell's response range) was extrapolated
# flat at its last plotted value -- override with --etfe-transmission if
# better data turns up. The table is kept here rather than just the final
# scalar so the derivation can be audited or redone.
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
# replacing it. Flat figure (user-supplied, 2026-08-25), not spectrally-
# weighted like ETFE above -- no spectral transmission chart was provided.
DEFAULT_POE_TRANSMISSION = 0.92

GPS_TOPIC = "vehicle_gps_position"

# Heading derivation (see "Heading and panel pointing" in the docstring).
DEFAULT_HEADING_SMOOTHING_S = 10.0
MIN_GROUND_SPEED_M_S = 3.0            # below this, consecutive-fix bearing is GPS noise, not travel
EARTH_RADIUS_M = 6_371_000.0

# Solar panel surface normals, in the aircraft's PX4 FRD body frame (+X
# forward/nose, +Y right/starboard, +Z down).
#
# The two MPPT strings are mounted at slightly different angles, so each gets
# its own normal, given directly as a CAD "Direction vector" (Y=0, unstated
# both times, in both cases) in a frame where +X points toward the TAIL and +Z
# points UP -- the opposite sign convention from PX4 FRD on both axes. Since
# both frames are right-handed and share the same three physical axes
# (longitudinal, spanwise, vertical), flipping X and Z but not Y is the only
# PROPER rotation (determinant +1, a 180-degree turn about Y) that reconciles
# them -- flipping all three, or only one, would mirror rather than rotate,
# which two right-handed frames on the same rigid body can never be related by:
def _panel_normal_body(x_cad: float, z_cad: float) -> np.ndarray:
    v = np.array([-x_cad, 0.0, -z_cad])
    return v / np.linalg.norm(v)


PANEL_NORMAL_BODY_STRING_0 = _panel_normal_body(x_cad=0.263, z_cad=0.965)  # magnitude ~1.0002
PANEL_NORMAL_BODY_STRING_1 = _panel_normal_body(x_cad=0.144, z_cad=0.99)   # magnitude ~1.0004


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
# ULog loading -- GPS track only
# --------------------------------------------------------------------------
def load_gps(ulog: ULog) -> pd.DataFrame:
    """The whole log input for this script: valid 3D fixes with absolute UTC time."""
    matches = [d for d in ulog.data_list if d.name == GPS_TOPIC]
    if not matches:
        raise RuntimeError(f"Topic '{GPS_TOPIC}' not found in log")
    gps = matches[0]

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
    if df.empty:
        raise RuntimeError("No valid 3D GPS fixes with absolute UTC time in this log")
    df["time"] = pd.to_datetime(df["time_utc_us"], unit="us", utc=True)
    df = df.drop(columns=["fix_type", "time_utc_us"])
    return df.drop_duplicates("time").set_index("time").sort_index()


# --------------------------------------------------------------------------
# Track-derived heading
# --------------------------------------------------------------------------
def track_heading_deg(df: pd.DataFrame, smoothing_s: float) -> pd.Series:
    """Course over ground [deg clockwise from north], from consecutive fixes.

    Great-circle INITIAL bearing from each fix to the next:

        theta = atan2(sin(dlon) cos(lat2),
                      cos(lat1) sin(lat2) - sin(lat1) cos(lat2) cos(dlon))

    Not the simpler flat-earth atan2(dE, dN): at this aircraft's latitudes the
    two agree closely over a one-second step, but the great-circle form costs
    nothing extra and does not degrade at high latitude.

    Smoothing is CIRCULAR -- a plain rolling mean of degrees would average
    359 and 1 to 180, pointing the panel exactly backwards every time the
    aircraft crosses north. Averaging the unit vectors (sin/cos) and taking
    atan2 of the result has no such wrap discontinuity.

    Bearing is undefined when the aircraft has not moved, and GPS jitter at
    rest produces uniformly-distributed garbage bearings. Samples below
    MIN_GROUND_SPEED_M_S are dropped before smoothing and inherit the last
    good heading afterward (ffill, then bfill for a lead-in of taxi/stationary
    samples before the first real motion).
    """
    lat = np.radians(df["lat"].values)
    lon = np.radians(df["lon"].values)
    dt_s = df.index.to_series().diff().shift(-1).dt.total_seconds().values

    # Forward differences: each fix's bearing is toward the NEXT fix, so the
    # last sample has none and reuses its predecessor's via the ffill below.
    lat1, lat2 = lat[:-1], lat[1:]
    dlon = lon[1:] - lon[:-1]
    bearing = np.full(len(df), np.nan)
    bearing[:-1] = np.degrees(np.arctan2(
        np.sin(dlon) * np.cos(lat2),
        np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon),
    )) % 360.0

    # Great-circle distance over the same step, for the stationary gate.
    hav = (np.sin((lat2 - lat1) / 2.0) ** 2
           + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2)
    step_m = np.full(len(df), np.nan)
    step_m[:-1] = 2.0 * EARTH_RADIUS_M * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0)))
    with np.errstate(invalid="ignore", divide="ignore"):
        speed_m_s = step_m / dt_s
    bearing[~(speed_m_s >= MIN_GROUND_SPEED_M_S)] = np.nan

    raw = pd.Series(bearing, index=df.index)
    window = f"{smoothing_s:.0f}s"
    sin_mean = np.sin(np.radians(raw)).rolling(window, min_periods=1).mean()
    cos_mean = np.cos(np.radians(raw)).rolling(window, min_periods=1).mean()
    heading = np.degrees(np.arctan2(sin_mean, cos_mean)) % 360.0
    return heading.ffill().bfill()


# --------------------------------------------------------------------------
# Atmosphere: ISA temperature and pressure
# --------------------------------------------------------------------------
def isa_temperature_k(alt_m: np.ndarray) -> np.ndarray:
    """ISA standard-atmosphere temperature [K], piecewise to 32 km.

    Three layers: troposphere (linear -6.5 K/km to 11 km), tropopause
    (isothermal 216.65 K to 20 km), lower stratosphere (linear +1 K/km above).
    32 km covers any altitude this aircraft flies.

    GLOBAL AVERAGE ONLY -- a function of altitude and nothing else. See
    "IMPORTANT SIMPLIFICATION -- ISA ambient temperature" in the module
    docstring for what this ignores (latitude, season, weather) and how much
    it matters downstream.
    """
    alt_m = np.asarray(alt_m, dtype=float)
    return np.where(
        alt_m <= 11000.0,
        288.15 - 0.0065 * alt_m,
        np.where(alt_m <= 20000.0, 216.65, 216.65 + 0.001 * (alt_m - 20000.0)),
    )


def isa_pressure_ratio(alt_m: np.ndarray) -> np.ndarray:
    """ISA static pressure as a fraction of sea level, piecewise to 32 km.

    p/p0 = (rho/rho0) * (T/T0) via the ideal gas law. Vectorized port of
    Icarus-Matrix/vehicle-simulation's air_density()/air_temperature_K()
    (endurance.py). Drives the Beer-Lambert optical depth in
    compute_clearsky_irradiance().
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
    return (density / 1.225) * (isa_temperature_k(alt_m) / 288.15)


# --------------------------------------------------------------------------
# Solar geometry / clear-sky irradiance
# --------------------------------------------------------------------------
def compute_clearsky_irradiance(times: pd.DatetimeIndex, lat: np.ndarray,
                                 lon: np.ndarray, alt_m: np.ndarray) -> pd.DataFrame:
    """Per-sample clear-sky irradiance at the aircraft's own position/time/altitude.

    Beer-Lambert attenuation (see SEA_LEVEL_TRANSMITTANCE for why, not pvlib's
    Ineichen/Perez model). Diffuse sky radiation is NOT modeled separately --
    this is an all-beam/no-diffuse split, which was a ~1% contribution even
    under the old Ineichen model at these altitudes. ghi_w_m2 (flat plate
    facing straight up) is kept as a reference curve and as the
    --assume-horizontal fallback; dni_w_m2 is what the POA projection uses.
    """
    solpos = pvlib.solarposition.get_solarposition(times, lat, lon, altitude=alt_m)
    cos_zenith = np.cos(np.radians(solpos["apparent_zenith"].values))
    sun_up = cos_zenith > 0.0

    dni_extra = np.asarray(pvlib.irradiance.get_extra_radiation(times), dtype=float)
    tau = SEA_LEVEL_TRANSMITTANCE ** isa_pressure_ratio(alt_m)

    # Plane-parallel airmass (1/cos z). Only meaningful where the sun is up;
    # elsewhere cos_zenith <= 0 would blow this up for no reason, since ghi/dni
    # are forced to 0 below regardless.
    safe_cos_zenith = np.where(sun_up, cos_zenith, 1.0)
    direct_normal_w_m2 = dni_extra * tau ** (1.0 / safe_cos_zenith)

    return pd.DataFrame(
        {
            "sun_elevation_deg": solpos["apparent_elevation"].values,
            "sun_azimuth_deg": solpos["azimuth"].values,
            "ghi_w_m2": np.where(sun_up, direct_normal_w_m2 * cos_zenith, 0.0),
            "dni_w_m2": np.where(sun_up, direct_normal_w_m2, 0.0),
        },
        index=times,
    )


def panel_normal_ned(roll_deg: np.ndarray, pitch_deg: np.ndarray, yaw_deg: np.ndarray,
                      normal_body: np.ndarray) -> np.ndarray:
    """Rotate a body-frame panel normal into the NED earth frame, per sample.

    Standard aerospace body-to-NED rotation, 3-2-1 (yaw-pitch-roll) Euler
    sequence: R = Rz(yaw) @ Ry(pitch) @ Rx(roll), applied to the constant
    body-frame vector in closed form rather than building an (N, 3, 3) stack.

    Returns an (N, 3) array of unit vectors in NED (North, East, Down).
    """
    r, p, y = np.radians(roll_deg), np.radians(pitch_deg), np.radians(yaw_deg)
    cr, sr = np.cos(r), np.sin(r)
    cp, sp = np.cos(p), np.sin(p)
    cy, sy = np.cos(y), np.sin(y)
    nx, ny, nz = normal_body

    ned_n = (cy * cp) * nx + (cy * sp * sr - sy * cr) * ny + (cy * sp * cr + sy * sr) * nz
    ned_e = (sy * cp) * nx + (sy * sp * sr + cy * cr) * ny + (sy * sp * cr - cy * sr) * nz
    ned_d = (-sp) * nx + (cp * sr) * ny + (cp * cr) * nz
    return np.stack([ned_n, ned_e, ned_d], axis=-1)


def sun_direction_ned(elevation_deg: np.ndarray, azimuth_deg: np.ndarray) -> np.ndarray:
    """Unit vector FROM the aircraft TOWARD the sun, in the NED earth frame.

    pvlib's azimuth (clockwise from north: 0=N, 90=E, 180=S, 270=W) maps
    directly onto NED's N/E axes. NED's "down" is positive, so "toward the
    sun" (above the horizon) has a NEGATIVE down-component.
    """
    el, az = np.radians(elevation_deg), np.radians(azimuth_deg)
    return np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), -np.sin(el)], axis=-1)


def cos_incidence_angle(heading_deg: np.ndarray, elevation_deg: np.ndarray,
                         azimuth_deg: np.ndarray, normal_body: np.ndarray) -> np.ndarray:
    """cos(AOI) between a wings-level tilted panel and the sun.

    Roll and pitch are pinned to zero -- see "Heading and panel pointing" in
    the module docstring for why, and for what that costs during banked turns.
    Clipped at 0 where the sun is behind the panel (no direct beam reaches
    it), never negative.
    """
    zeros = np.zeros_like(heading_deg)
    normal = panel_normal_ned(zeros, zeros, heading_deg, normal_body)
    sun = sun_direction_ned(elevation_deg, azimuth_deg)
    return np.clip(np.einsum("ij,ij->i", normal, sun), 0.0, 1.0)


# --------------------------------------------------------------------------
# Array power model
# --------------------------------------------------------------------------
def track_from_ulog(args: argparse.Namespace) -> pd.DataFrame:
    """Flown track: GPS fixes from a log, with heading derived from the path."""
    ulog_path = Path(args.ulog)
    if not ulog_path.exists():
        raise FileNotFoundError(f"ULog file not found: {ulog_path}")

    print(f"Loading {ulog_path.name} ...")
    ulog = ULog(str(ulog_path), message_name_filter_list=[GPS_TOPIC])
    df = load_gps(ulog)
    print(f"  GPS fixes: {len(df)}")
    df["heading_deg"] = track_heading_deg(df, args.heading_smoothing_s)
    return df


TRACK_CSV_REQUIRED_COLS = ("time", "lat", "lon", "alt_msl_m")


def track_from_csv(args: argparse.Namespace) -> pd.DataFrame:
    """Track from a plain CSV -- for planned/hypothetical flights, or any
    trajectory that didn't come out of a PX4 log.

    Required columns: time (any pandas-parseable timestamp; naive values are
    read as UTC), lat, lon, alt_msl_m. Optional: heading_deg -- supply it to
    pin the panel's pointing explicitly (a planned cruise heading), or omit it
    to have the heading derived from the track the same way a flown log's is.

    A track that holds one fixed lat/lon (a station-keeping stand-in) has no
    motion to derive heading from, so heading_deg is effectively REQUIRED
    there -- without it every bearing falls below MIN_GROUND_SPEED_M_S and the
    whole column resolves to a single fallback value.
    """
    csv_path = Path(args.track_csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"Track CSV not found: {csv_path}")

    print(f"Loading {csv_path.name} ...")
    df = pd.read_csv(csv_path)
    missing = [c for c in TRACK_CSV_REQUIRED_COLS if c not in df.columns]
    if missing:
        raise RuntimeError(
            f"Track CSV is missing required column(s): {', '.join(missing)}. "
            f"Required: {', '.join(TRACK_CSV_REQUIRED_COLS)}; optional: heading_deg."
        )

    df["time"] = pd.to_datetime(df["time"], utc=True)
    df = df.drop_duplicates("time").set_index("time").sort_index()
    if df.empty:
        raise RuntimeError(f"Track CSV has no rows: {csv_path}")
    print(f"  Track points: {len(df)}")

    if "heading_deg" in df.columns:
        print("  Using heading_deg from the CSV (not derived from the track)")
    else:
        df["heading_deg"] = track_heading_deg(df, args.heading_smoothing_s)
    return df


def apply_power_model(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    """Track -> power. Expects lat/lon/alt_msl_m/heading_deg on a UTC index;
    see the module docstring for the full chain.
    """
    print("Computing solar position (pvlib) + clear-sky irradiance (Beer-Lambert) ...")
    df = df.join(compute_clearsky_irradiance(
        df.index, df["lat"].values, df["lon"].values, df["alt_msl_m"].values
    ))

    # Per-string POA. String 0 and string 1 sit at different CAD tilt angles,
    # so they see different irradiance at the same instant -- worth keeping
    # separate rather than picking one normal to stand for the whole array.
    strings = {
        0: (args.cell_count - args.string1_cell_count, PANEL_NORMAL_BODY_STRING_0),
        1: (args.string1_cell_count, PANEL_NORMAL_BODY_STRING_1),
    }
    for idx, (_, normal_body) in strings.items():
        if args.assume_horizontal:
            df[f"poa_string{idx}_w_m2"] = df["ghi_w_m2"]
        else:
            df[f"poa_string{idx}_w_m2"] = df["dni_w_m2"] * cos_incidence_angle(
                df["heading_deg"].values, df["sun_elevation_deg"].values,
                df["sun_azimuth_deg"].values, normal_body,
            )

    # Ambient air temperature from altitude alone (ISA). --ambient-offset-c
    # shifts the whole profile for a flight known to run warm or cold.
    df["t_ambient_c"] = isa_temperature_k(df["alt_msl_m"].values) - 273.15 + args.ambient_offset_c

    encapsulation_transmission = args.etfe_transmission * args.poe_transmission
    cell_area_m2 = args.cell_area_cm2 / 1e4

    for idx, (cell_count, _) in strings.items():
        poa = df[f"poa_string{idx}_w_m2"]

        # NOCT: cell temperature rises above ambient in proportion to
        # irradiance. Zero irradiance -> cell sits at ambient.
        df[f"t_cell_string{idx}_c"] = df["t_ambient_c"] + (poa / NOCT_IRRADIANCE_W_M2) * (
            args.noct_c - NOCT_AMBIENT_C
        )

        # Temperature derate. The coefficient is NEGATIVE, so cells below
        # STC_TEMP_C (which is most of this flight) give a factor ABOVE 1 --
        # a gain, not a loss. Floored at 0 so a pathological input can shrink
        # power toward zero but never invert its sign.
        derate = 1.0 + (args.power_temp_coeff / 100.0) * (df[f"t_cell_string{idx}_c"] - STC_TEMP_C)
        df[f"temp_derate_string{idx}"] = derate.clip(lower=0.0)

        df[f"p_string{idx}_w"] = (
            poa * cell_count * cell_area_m2 * args.cell_efficiency
            * encapsulation_transmission * df[f"temp_derate_string{idx}"]
        )

    df["p_array_w"] = df["p_string0_w"] + df["p_string1_w"]

    df.attrs["strings"] = {i: n for i, (n, _) in strings.items()}
    df.attrs["total_area_m2"] = args.cell_count * cell_area_m2
    df.attrs["encapsulation_transmission"] = encapsulation_transmission
    return df


def print_summary(df: pd.DataFrame, args: argparse.Namespace, tz: str) -> None:
    # Energy integration: >10 s gaps are dropped rather than trapezoided
    # across, so a logging dropout can't manufacture energy that was never
    # sampled.
    dt_hours = df.index.to_series().diff().dt.total_seconds().fillna(0.0) / 3600.0
    dt_hours = dt_hours.clip(upper=10.0 / 3600.0)
    energy_wh = (df["p_array_w"] * dt_hours).sum()

    daylight = df[df["sun_elevation_deg"] > 0]
    rated_stc_w = df.attrs["total_area_m2"] * args.cell_efficiency * STC_IRRADIANCE_W_M2

    local_start = df.index[0].tz_convert(tz)
    local_end = df.index[-1].tz_convert(tz)
    counts = df.attrs["strings"]

    print("\n" + "=" * 70)
    print("ESTIMATED SOLAR ARRAY OUTPUT  (modeled from flight path only)")
    print("=" * 70)
    print(f"Launch-site timezone     : {tz}")
    print(f"Flight window (local)    : {local_start.strftime('%Y-%m-%d %H:%M:%S %Z')} "
          f"-> {local_end.strftime('%H:%M:%S %Z')}")
    print(f"Flight duration          : {df.index[-1] - df.index[0]}")
    print(f"Samples                  : {len(df)}  ({len(daylight)} with sun above horizon)")
    print(f"Altitude range           : {df['alt_msl_m'].min():.0f} - {df['alt_msl_m'].max():.0f} m MSL")
    print("-" * 70)
    print(f"Array                    : {args.cell_count} cells x {args.cell_area_cm2:.0f} cm^2 "
          f"@ {args.cell_efficiency * 100:.1f}%  -> {df.attrs['total_area_m2']:.3f} m^2")
    print(f"  string 0 / string 1    : {counts[0]} / {counts[1]} cells")
    print(f"Encapsulation            : ETFE {args.etfe_transmission * 100:.1f}% x POE "
          f"{args.poe_transmission * 100:.1f}% = {df.attrs['encapsulation_transmission'] * 100:.1f}%")
    print(f"STC-rated array power    : {rated_stc_w:.1f} W  (1000 W/m^2, 25 degC, bare cells)")
    print("-" * 70)
    if args.assume_horizontal:
        print("Panel pointing           : disabled (--assume-horizontal) -- POA == flat-plate GHI")
    else:
        # Positive = tilt+heading GAIN vs a flat plate, negative = loss. No
        # fixed sign: it depends entirely on this flight's track vs the sun.
        lit = df["ghi_w_m2"] > 1.0
        poa_mean = df.loc[lit, ["poa_string0_w_m2", "poa_string1_w_m2"]].mean(axis=1)
        change_pct = 100.0 * (poa_mean / df.loc[lit, "ghi_w_m2"] - 1.0)
        print(f"Panel pointing           : mean {change_pct.mean():+.1f}% vs flat-plate GHI  "
              f"(range [{change_pct.min():+.1f}%, {change_pct.max():+.1f}%], "
              f"heading + CAD tilt, wings-level)")
    print(f"Peak modeled GHI (flat)  : {df['ghi_w_m2'].max():.1f} W/m^2")
    print(f"Peak modeled POA         : string 0 {df['poa_string0_w_m2'].max():.1f}  "
          f"string 1 {df['poa_string1_w_m2'].max():.1f} W/m^2")
    print(f"Peak sun elevation       : {df['sun_elevation_deg'].max():.1f} deg")
    print("-" * 70)
    print(f"Ambient temp (ISA)       : {df['t_ambient_c'].min():.1f} to "
          f"{df['t_ambient_c'].max():.1f} degC")
    print(f"Cell temp (NOCT {args.noct_c:.0f} degC) : string 0 "
          f"[{df['t_cell_string0_c'].min():.1f}, {df['t_cell_string0_c'].max():.1f}]  "
          f"string 1 [{df['t_cell_string1_c'].min():.1f}, {df['t_cell_string1_c'].max():.1f}] degC")
    # Coefficient is negative and these cells run cold, so this is normally a
    # GAIN -- report it with an explicit direction rather than calling it a loss.
    mean_derate = df[["temp_derate_string0", "temp_derate_string1"]].mean(axis=1).mean()
    change_pct = 100.0 * (mean_derate - 1.0)
    print(f"Temp derate effect       : mean {abs(change_pct):.1f}% "
          f"{'gain' if change_pct >= 0 else 'loss'}  "
          f"(STC {STC_TEMP_C:.0f} degC, coeff {args.power_temp_coeff}%/degC)")
    print("-" * 70)
    print(f"Peak estimated power     : {df['p_array_w'].max():.1f} W")
    if len(daylight):
        print(f"Mean estimated power     : {daylight['p_array_w'].mean():.1f} W  (sun above horizon)")
    print(f"Estimated energy         : {energy_wh:.1f} Wh")
    print("=" * 70)
    print("MODELED, NOT MEASURED. Clear-sky irradiance cannot represent cloud")
    print("cover, so this is an upper envelope, not a forecast of actual output.")
    print("Bank angle is not modeled (wings-level assumption), so instantaneous")
    print("power is far less trustworthy than orbit-averaged energy. Ambient temp")
    print("is ISA global-average -- no latitude, season, or weather. See the")
    print("module docstring for the full error budget.")


def print_hourly_table(df: pd.DataFrame, tz: str) -> None:
    """Hourly array output + cumulative energy, on the local clock.

    Cumulative energy is TRAPEZOIDAL over the hourly samples, matching how a
    spreadsheet of hourly readings would normally be integrated -- so this
    column is directly comparable against one. It is not the same number as
    the summary's energy figure, which integrates every sample: hourly
    trapezoid under-reads a curve that is concave near solar noon and
    misplaces the sunrise/sunset ends, and the gap between the two is a
    readout of how much resolution the hourly grid throws away.
    """
    local = df.set_axis(df.index.tz_convert(tz))
    hourly = local["p_array_w"].resample("1h").mean()
    if hourly.empty:
        return

    # Trapezoid on the hourly grid: each step contributes the mean of its two
    # endpoints, in Wh since the step is exactly one hour.
    steps = (hourly + hourly.shift()).div(2.0).fillna(0.0)
    cumulative = steps.cumsum()

    print("\n" + "=" * 52)
    print("HOURLY OUTPUT")
    print("=" * 52)
    print(f"{'Local time':<12}{'Hour (dec)':>12}{'Output (W)':>14}{'Cum. (Wh)':>14}")
    print("-" * 52)
    for ts, watts in hourly.items():
        hour_dec = ts.hour + ts.minute / 60.0
        print(f"{ts.strftime('%H:%M'):<12}{hour_dec:>12.1f}{watts:>14.0f}{cumulative[ts]:>14.0f}")
    print("-" * 52)
    print(f"{'Daily total':<12}{'':>12}{'':>14}{cumulative.iloc[-1]:>14.0f}")
    print("=" * 52)


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
def make_plot(df: pd.DataFrame, out_path: Path, tz: str) -> None:
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    df = df.set_axis(df.index.tz_convert(tz))
    fig, axes = plt.subplots(4, 1, figsize=(12, 13), sharex=True)

    ax = axes[0]
    ax.plot(df.index, df["alt_msl_m"], color="tab:gray", linewidth=1.2, label="Altitude (MSL)")
    ax.set_ylabel("Altitude (m)")
    ax.set_title("Flight Path", fontweight="bold")
    ax_sun = ax.twinx()
    ax_sun.plot(df.index, df["sun_elevation_deg"], color="tab:orange", linewidth=1.2,
                label="Sun elevation")
    ax_sun.axhline(0.0, color="tab:orange", linestyle=":", linewidth=0.8)
    ax_sun.set_ylabel("Sun Elevation (deg)", color="tab:orange")
    ax_sun.tick_params(axis="y", labelcolor="tab:orange")

    ax = axes[1]
    ax.plot(df.index, df["ghi_w_m2"], color="tab:gray", linewidth=1.0, linestyle="--",
            label="GHI (flat plate)")
    ax.plot(df.index, df["poa_string0_w_m2"], color="tab:blue", linewidth=1.0, label="POA string 0")
    ax.plot(df.index, df["poa_string1_w_m2"], color="tab:green", linewidth=1.0, label="POA string 1")
    ax.set_ylabel("Irradiance (W/m$^2$)")
    ax.set_title("Clear-Sky Irradiance (modeled -- no cloud cover)", fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)

    ax = axes[2]
    ax.plot(df.index, df["t_ambient_c"], color="tab:purple", linewidth=1.2, label="Ambient (ISA)")
    ax.plot(df.index, df["t_cell_string0_c"], color="tab:blue", linewidth=1.0, label="Cell string 0")
    ax.plot(df.index, df["t_cell_string1_c"], color="tab:green", linewidth=1.0, label="Cell string 1")
    ax.axhline(STC_TEMP_C, color="tab:red", linestyle=":", linewidth=0.9,
               label=f"STC {STC_TEMP_C:.0f} degC")
    ax.set_ylabel("Temperature (degC)")
    ax.set_title("Ambient (ISA) and Cell (NOCT) Temperature", fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)

    ax = axes[3]
    ax.plot(df.index, df["p_string0_w"], color="tab:blue", linewidth=0.9, label="String 0")
    ax.plot(df.index, df["p_string1_w"], color="tab:green", linewidth=0.9, label="String 1")
    ax.plot(df.index, df["p_array_w"], color="tab:red", linewidth=1.4, label="Array total")
    ax.set_ylabel("Power (W)")
    ax.set_title("Estimated Array Output Power", fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=df.index.tz))
    ax.set_xlabel(f"Local time, {tz} ({df.index[0].date()})")

    for a in axes:
        a.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot -> {out_path}")


def open_in_vscode(path: Path) -> None:
    """Open the generated plot in VS Code, reusing the existing window.

    Best-effort: if the `code` CLI isn't on PATH, prints a note rather than
    failing the run. Resolved via shutil.which() because on Windows the
    launcher is a .CMD shim and CreateProcess (which subprocess uses without
    shell=True) won't apply PATHEXT resolution to find it.
    """
    code_cmd = shutil.which("code")
    if code_cmd is None:
        print("  (note: 'code' CLI not found on PATH - skipping VS Code auto-open)")
        return
    try:
        subprocess.run([code_cmd, "-r", str(path)], check=True)
    except (subprocess.CalledProcessError, OSError) as exc:
        print(f"  (note: could not auto-open {path.name} in VS Code: {exc})")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--ulog", help="Path to a flown PX4 .ulg log (default input if "
                                        "neither --ulog nor --track-csv is given)")
    source.add_argument("--track-csv",
                         help="Path to a planned/hypothetical track CSV instead of a log. "
                              f"Required columns: {', '.join(TRACK_CSV_REQUIRED_COLS)}; "
                              "optional heading_deg (derived from the track if absent -- but "
                              "effectively required for a fixed-position track, which has no "
                              "motion to derive it from).")

    parser.add_argument("--cell-count", type=int, default=DEFAULT_CELL_COUNT,
                         help="Total cells across the whole array")
    parser.add_argument("--string1-cell-count", type=int, default=DEFAULT_STRING1_CELL_COUNT,
                         help="Cells wired into string 1 (NOT half of --cell-count -- the two "
                              "strings aren't equal size). String 0 gets the remainder.")
    parser.add_argument("--cell-area-cm2", type=float, default=DEFAULT_CELL_AREA_CM2)
    parser.add_argument("--cell-efficiency", type=float, default=DEFAULT_CELL_EFFICIENCY,
                         help="Fractional cell efficiency, e.g. 0.254 for 25.4%%")
    parser.add_argument("--etfe-transmission", type=float, default=DEFAULT_ETFE_TRANSMISSION,
                         help="Fractional light transmission through the ETFE cover, "
                              "spectrally-weighted by cell EQE x AM1.5G spectrum. "
                              f"Default {DEFAULT_ETFE_TRANSMISSION:.2f}.")
    parser.add_argument("--poe-transmission", type=float, default=DEFAULT_POE_TRANSMISSION,
                         help="Fractional light transmission through the POE encapsulant "
                              f"(stacks with --etfe-transmission). Default {DEFAULT_POE_TRANSMISSION:.2f}.")

    parser.add_argument("--noct-c", type=float, default=DEFAULT_NOCT_C,
                         help="Nominal Operating Cell Temperature [degC] for the cell-temp "
                              "model. The default is the conventional open-rack placeholder, "
                              "NOT a figure derived for this airframe -- flight-speed airflow "
                              f"means the real value is likely LOWER. Default {DEFAULT_NOCT_C:.0f}.")
    parser.add_argument("--power-temp-coeff", type=float, default=POWER_TEMP_COEFF_PCT_PER_C,
                         help="Power temperature coefficient [%%/degC, negative]. Cold cells "
                              f"OVERperform nameplate. Default {POWER_TEMP_COEFF_PCT_PER_C}.")
    parser.add_argument("--ambient-offset-c", type=float, default=0.0,
                         help="Shift the whole ISA ambient-temperature profile by this many "
                              "degC. ISA is a global average with no latitude/season/weather "
                              "term; use this if the flight is known to have run warm or cold.")

    parser.add_argument("--heading-smoothing-s", type=float, default=DEFAULT_HEADING_SMOOTHING_S,
                         help="Circular-mean smoothing window for the track-derived heading "
                              f"[s]. Default {DEFAULT_HEADING_SMOOTHING_S:.0f}.")
    parser.add_argument("--assume-horizontal", action="store_true",
                         help="Ignore heading and panel tilt; assume the array always faces "
                              "straight up (POA == GHI). Panel pointing is ON by default.")

    parser.add_argument("--output-dir", default=None,
                         help="Where to write the CSV/plot (default: alongside the ulog file)")
    parser.add_argument("--no-plot", action="store_true", help="Skip generating the PNG plot")
    parser.add_argument("--no-open", action="store_true",
                         help="Don't auto-open the plot in VS Code after saving it")
    parser.add_argument("--tz", default=None,
                         help="Timezone for the plot's time axis (IANA name). Default: "
                              "auto-detected from the launch-site GPS fix. Data is exported in UTC.")
    args = parser.parse_args()

    if args.track_csv:
        source_path = Path(args.track_csv)
        df = track_from_csv(args)
    else:
        source_path = Path(args.ulog or DEFAULT_ULOG)
        args.ulog = str(source_path)
        df = track_from_ulog(args)

    df = apply_power_model(df, args)
    tz = args.tz or detect_launch_timezone(df["lat"].iloc[0], df["lon"].iloc[0])
    print_summary(df, args, tz)
    print_hourly_table(df, tz)

    out_dir = Path(args.output_dir) if args.output_dir else source_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = source_path.stem

    csv_path = out_dir / f"{stem}_array_power.csv"
    df[[
        "lat", "lon", "alt_msl_m", "heading_deg",
        "sun_elevation_deg", "sun_azimuth_deg", "ghi_w_m2", "dni_w_m2",
        "poa_string0_w_m2", "poa_string1_w_m2",
        "t_ambient_c", "t_cell_string0_c", "t_cell_string1_c",
        "temp_derate_string0", "temp_derate_string1",
        "p_string0_w", "p_string1_w", "p_array_w",
    ]].to_csv(csv_path)
    print(f"Saved data -> {csv_path}")

    if not args.no_plot:
        plot_path = out_dir / f"{stem}_array_power.png"
        make_plot(df, plot_path, tz)
        if not args.no_open:
            open_in_vscode(plot_path)


if __name__ == "__main__":
    main()
