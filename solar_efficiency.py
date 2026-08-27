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
  - /zeus/flight           -> roll_deg, pitch_deg, yaw_deg. Used for the panel
                              incidence-angle correction (see "Panel incidence
                              angle" below); skipped under --assume-horizontal.

Cell reference data (Maxeon Gen 7 datasheet, 546209 Rev C):
  - Cell area       : ~155 cm^2 per cell
  - Cell efficiency : ~25.4% (Pe/Oe bin boundary -> "typical" cell)
  - Power temp coeff: -0.27 %/degC, relative to STC_TEMP_C (25 degC). Applied
    as an optional adjustment on the estimated power when --apply-temp-derate
    is passed (off by default -- see "Temperature derating" below). Despite
    the name, this is NOT always a loss: since the coefficient is negative,
    Tout above 25 degC reduces estimated power (a loss) but Tout BELOW
    25 degC increases it (a gain). The code, plot, and summary all handle
    both signs -- don't assume derating only ever shrinks power.

Estimated/available array power at each instant is layered in three tiers,
each one a further loss on top of the last (see "Panel incidence angle",
"Encapsulation transmission", and "Temperature derating" below for each term):

    P_bare(t)        = POA_irradiance(t) [W/m^2] * total_cell_area [m^2] * cell_efficiency
    P_no_temp_derate(t) = P_bare(t) * etfe_transmission * poe_transmission
    P_estimated(t) = P_no_temp_derate(t) [* temp_derate_factor(t)]

P_bare is the cell nameplate ceiling with nothing between the sun and the
cells -- not physically real for this array, but a useful "no losses at
all" reference. P_no_temp_derate adds the encapsulation stack's light loss
(ETFE cover + POE encapsulant, both always applied, not opt-in) and is what
"efficiency" is actually measured against when --apply-temp-derate is off.
P_estimated is the final number used everywhere once temp derating is
added on top.

POA_irradiance (plane-of-array) is the clear-sky direct beam projected onto
the ACTUAL panel orientation -- not a flat horizontal assumption -- unless
--assume-horizontal is passed, in which case it falls back to GHI_clearsky
(irradiance on a plate facing straight up). See "Panel incidence angle".

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

Panel incidence angle (roll/pitch/yaw, always applied -- opt OUT with
--assume-horizontal, the reverse of --apply-temp-derate's opt-IN):
The array is mounted on a wing, not a fixed ground rack, so it's tilted
by whatever the aircraft's attitude happens to be at each instant --
banking into a turn can point it more squarely at the sun than a flat
plate ever would, or away from it entirely. The two MPPT strings are
mounted at slightly different angles, so each has its own surface normal
in the aircraft's body frame -- PANEL_NORMAL_BODY_STRING_0 and _STRING_1
(see that constant's comment for the full derivation from each string's
CAD "direction vector" and how the shared coordinate convention was
confirmed, not assumed, against the log's own attitude data). Each
sample: panel_normal_ned() rotates a given body-frame normal into the
earth (NED) frame using /zeus/flight's roll/pitch/yaw at that instant;
cos_incidence_angle() then dots it against the sun's direction
(sun_direction_ned()) to get cos(AOI), clipped to 0 when the sun is
behind the panel. POA_irradiance = dni_w_m2 * cos(AOI) -- reusing the
same clear-sky direct-beam value already computed for GHI, just projected
onto the real orientation instead of assuming flat-horizontal. Samples
with no attitude match fall back to cos(zenith) (flat) rather than NaN.
The array-wide poa_w_m2 (efficiency calc, main plot) uses string 1's
normal specifically; poa_string0_w_m2/poa_string1_w_m2 (see
make_poa_strings_plot()) compute each string's POA separately for
comparison -- combining both into a single per-string-weighted array POA
would be a real improvement but isn't done here. --assume-horizontal
restores the old flat-plate behavior exactly (useful for comparison, or
logs without /zeus/flight) by setting every POA variant equal to GHI.

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
    python solar_efficiency.py --ulog log.ulg --assume-horizontal
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
DEFAULT_STRING1_CELL_COUNT = 58       # of the 72 total -- string 1 is NOT half the array (user, 2026-08-25)
DEFAULT_CELL_AREA_CM2 = 155.0        # per cell
DEFAULT_CELL_EFFICIENCY = 0.254      # 25.4%, typical Pe/Oe bin boundary
POWER_TEMP_COEFF_PCT_PER_C = -0.27   # informational, not applied by default
STC_TEMP_C = 25.0
STC_IRRADIANCE_W_M2 = 1000.0

# Genasun MPPT-12SBB datasheet input voltage window.
MPPT_12SBB_VIN_MIN_V = 8.0
MPPT_12SBB_VIN_MAX_V = 51.0

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
TOUT_CLAMP_RANGE_C = (-100.0, 85.0)                # defensive: generous aviation/electronics envelope --
                                                    # floor wide enough to cover genuine stratospheric
                                                    # cruise temps (this aircraft flies ~55-60kft, where
                                                    # ISA predicts well below -40 degC) without catching
                                                    # real readings as the dead-sensor sentinel
FLIGHT_TOPIC = "/zeus/flight"                      # roll_deg/pitch_deg/yaw_deg -- see load_attitude()

# Solar panel surface normals, in the aircraft's PX4 FRD body frame (+X
# forward/nose, +Y right/starboard, +Z down). Confirmed (2026-08-25) to be
# /zeus/flight's own roll/pitch/yaw convention: cross-checked its reported
# angles against Euler angles decomposed from vehicle_attitude's quaternion
# (standard aerospace 3-2-1 sequence) across 122k matched samples on
# 00007-2026-08-24_20-17-42.ulg -- r = 0.94-0.99, near-identical values.
#
# The two MPPT strings are mounted at slightly different angles, so each
# gets its own normal, given directly as a CAD "Direction vector" (Y=0,
# unstated both times, in both cases -- see string 1's original screenshot)
# in a frame where +X points toward the TAIL and +Z points UP -- the
# opposite sign convention from PX4 FRD on both axes. Since both frames are
# right-handed and share the same three physical axes (longitudinal,
# spanwise, vertical), flipping X and Z but not Y is the only PROPER
# rotation (determinant +1, a 180-degree turn about Y) that reconciles them
# -- flipping all three, or only one, would mirror rather than rotate,
# which two right-handed frames on the same rigid body can never be
# related by:
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


def load_attitude(ulog: ULog) -> pd.DataFrame | None:
    """Aircraft roll/pitch/yaw, used for the panel incidence-angle correction.

    /zeus/flight carries its own absolute-UTC timestamp_us field, the same
    convention as /zeus/mppt_* and /zeus/temperature. Its roll_deg/
    pitch_deg/yaw_deg are standard PX4 FRD body-frame Euler angles -- see
    PANEL_NORMAL_BODY_STRING_0/_1's comment for how that was confirmed, not assumed.
    """
    f = _get_dataset(ulog, FLIGHT_TOPIC)
    if f is None:
        print(f"  (note: {FLIGHT_TOPIC} not present in this log - panel incidence-angle "
              f"correction disabled, falling back to a flat horizontal assumption)")
        return None
    df = pd.DataFrame(
        {
            "time_utc_us": f.data["timestamp_us"],
            "roll_deg": f.data["roll_deg"],
            "pitch_deg": f.data["pitch_deg"],
            "yaw_deg": f.data["yaw_deg"],
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

    dni_w_m2/dhi_w_m2 are a simplified all-beam/no-diffuse split -- this
    model doesn't compute diffuse sky radiation separately (it was a ~1%
    contribution even under the old Ineichen model at this flight's
    altitudes). ghi_w_m2 assumes a flat HORIZONTAL surface (irradiance on
    a plate facing straight up); it's kept as a reference curve, but
    poa_w_m2 (derived from dni_w_m2 via panel_normal_ned(), see below) is
    what actually drives the estimated-power calc unless
    --assume-horizontal is passed.
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


def panel_normal_ned(roll_deg: np.ndarray, pitch_deg: np.ndarray, yaw_deg: np.ndarray,
                      normal_body: np.ndarray) -> np.ndarray:
    """Rotate a body-frame panel normal (e.g. PANEL_NORMAL_BODY_STRING_0/_1)
    into the NED earth frame, one rotation per sample.

    Standard aerospace body-to-NED rotation, 3-2-1 (yaw-pitch-roll) Euler
    sequence: R = Rz(yaw) @ Ry(pitch) @ Rx(roll). This is the same sequence
    used to confirm /zeus/flight's convention against vehicle_attitude's
    quaternion (see PANEL_NORMAL_BODY_STRING_0/_1) -- applied here to the
    constant body-frame vector directly (closed-form per-sample), rather
    than building and multiplying an (N, 3, 3) stack of rotation matrices.

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


def cos_incidence_angle(roll_deg: np.ndarray, pitch_deg: np.ndarray, yaw_deg: np.ndarray,
                         elevation_deg: np.ndarray, azimuth_deg: np.ndarray,
                         normal_body: np.ndarray) -> np.ndarray:
    """cos(AOI) between the tilted, rotating panel and the sun -- 0 where the
    sun is behind the panel (no direct beam reaches it), never negative.

    normal_body selects which panel's normal to use (e.g.
    PANEL_NORMAL_BODY_STRING_0 or _STRING_1 -- the two MPPT strings are
    mounted at slightly different angles, see that constant's comment).

    Where attitude is NaN (no /zeus/flight match for that sample), falls
    back to cos(zenith) -- i.e. the flat-horizontal assumption -- rather
    than propagating NaN into the estimated-power calc for that row.
    """
    normal = panel_normal_ned(roll_deg, pitch_deg, yaw_deg, normal_body)
    sun = sun_direction_ned(elevation_deg, azimuth_deg)
    cos_aoi = np.einsum("ij,ij->i", normal, sun)
    cos_zenith_fallback = np.sin(np.radians(elevation_deg))  # cos(zenith) == sin(elevation)
    cos_aoi = np.where(np.isnan(cos_aoi), cos_zenith_fallback, cos_aoi)
    return np.clip(cos_aoi, 0.0, 1.0)


# --------------------------------------------------------------------------
# Flight phase (for plot markers -- holding high/low altitude, descending)
# --------------------------------------------------------------------------
FLIGHT_PHASE_RATE_THRESHOLD_M_S = 0.5   # |vertical rate| at/below this counts as "holding"
FLIGHT_PHASE_SMOOTHING_S = 60.0         # trailing time window for the rate, not a row count

# Rolling window for characterizing % Difference scatter (make_string1_plot).
# Long enough to average out per-sample noise in the efficiency ratio, short
# enough to still track real trend changes (e.g. cloud passage, altitude).
PCT_DIFF_ROLLING_WINDOW_S = 900.0


def classify_flight_phase(altitude_m: pd.Series) -> pd.Series:
    """Per-sample flight phase from the altitude time series alone: one of
    "climbing", "holding_high", "descending", "holding_low".

    Vertical rate is a TIME-based trailing rolling mean (pandas' '60s'
    window, not a fixed row count), since GPS altitude noise makes a raw
    sample-to-sample diff flap sign constantly, and this timeline's sample
    spacing isn't uniform. "High" vs "low" for a holding period is relative
    to THIS flight's own midpoint altitude ((min+max)/2), not an absolute
    threshold -- works whether the flight ceiling is 2 km or 20 km. Only
    "holding_high"/"descending"/"holding_low" are meant to be shown on
    plots (see shade_flight_phases()); "climbing" -- e.g. the initial
    climb-out -- is left as the unmarked default.
    """
    dt_s = altitude_m.index.to_series().diff().dt.total_seconds()
    rate_m_s = altitude_m.diff() / dt_s
    rate_smooth = rate_m_s.rolling(f"{FLIGHT_PHASE_SMOOTHING_S:.0f}s", min_periods=1).mean()

    mid_alt_m = (altitude_m.min() + altitude_m.max()) / 2.0
    is_holding = rate_smooth.abs() <= FLIGHT_PHASE_RATE_THRESHOLD_M_S
    conditions = [
        rate_smooth < -FLIGHT_PHASE_RATE_THRESHOLD_M_S,
        is_holding & (altitude_m >= mid_alt_m),
        is_holding & (altitude_m < mid_alt_m),
    ]
    choices = ["descending", "holding_high", "holding_low"]
    return pd.Series(np.select(conditions, choices, default="climbing"), index=altitude_m.index)


# Phase -> (background color, legend label). Only these three are drawn --
# "climbing" has no entry and is left unmarked (see classify_flight_phase()).
FLIGHT_PHASE_STYLE = {
    "holding_high": ("gold", "Holding High Altitude"),
    "descending": ("lightsteelblue", "Descending"),
    "holding_low": ("darkseagreen", "Holding Low Altitude"),
}


def shade_flight_phases(axes, phase: pd.Series) -> None:
    """Shade holding-high/descending/holding-low segments as low-zorder
    background spans on every axis in `axes`, so any time-series plot can
    show which flight phase each moment falls in. Draws on each axis
    independently (rather than sharing one span object) so each panel's own
    legend -- built separately per panel via legend_outside() -- picks up
    exactly one labeled handle per phase actually present, in that panel's
    own draw order.
    """
    change = phase.ne(phase.shift()).cumsum()
    seen_per_axis = {id(ax): set() for ax in axes}
    for _, seg in phase.groupby(change):
        name = seg.iloc[0]
        if name not in FLIGHT_PHASE_STYLE:
            continue
        color, label = FLIGHT_PHASE_STYLE[name]
        start, end = seg.index[0], seg.index[-1]
        for ax in axes:
            seen = seen_per_axis[id(ax)]
            ax.axvspan(start, end, color=color, alpha=0.12, zorder=0,
                       label=label if name not in seen else None)
            seen.add(name)


def main_hold_window(phase: pd.Series):
    """(zone1_start, zone3_end): the start of the LONGEST holding_high
    segment through the end of the LONGEST holding_low segment that starts
    after it -- i.e. the main climb-hold / descend / hold-low sequence
    ("zones 1-3"), not a brief pre-takeoff hold or a flickery landing-
    approach hold (both also classified holding_low, but much shorter).
    Picking the longest of each type, rather than the first/last, is what
    tells those apart. Returns (start, None) if no holding_low follows the
    main hold, or None entirely if there's no holding_high at all.
    """
    change = phase.ne(phase.shift()).cumsum()
    segments = [(seg.iloc[0], seg.index[0], seg.index[-1]) for _, seg in phase.groupby(change)]
    highs = [s for s in segments if s[0] == "holding_high"]
    if not highs:
        return None
    zone1_start = max(highs, key=lambda s: s[2] - s[1])[1]
    lows_after = [s for s in segments if s[0] == "holding_low" and s[1] > zone1_start]
    zone3_end = max(lows_after, key=lambda s: s[2] - s[1])[2] if lows_after else None
    return zone1_start, zone3_end


def keep_main_hold(df: pd.DataFrame) -> pd.DataFrame:
    """Keep ONLY zones 1-3 -- the main holding_high / descending /
    holding_low sequence, start to end (see main_hold_window()) -- dropping
    everything before zone 1 (climb-out) and after zone 3 (final approach/
    landing, which includes flickery descending/holding_low blips that
    share a phase label with zones 2/3 but aren't part of that contiguous
    run). One contiguous real-time slice, so no seam or re-indexing is
    needed -- the real clock-time index is kept as-is.
    """
    window = main_hold_window(df["flight_phase"])
    if window is None:
        return df
    zone1_start, zone3_end = window
    if zone3_end is None:
        return df[df.index >= zone1_start]
    return df[(df.index >= zone1_start) & (df.index <= zone3_end)]


# --------------------------------------------------------------------------
# Tout model: piecewise-linear target -> non-minimum-phase 2nd-order response
# --------------------------------------------------------------------------
# All three fitted parameters (T_z, tau_1, tau_2) are TIMES, so they're
# seeded and bounded relative to the descent's own duration rather than as
# absolute seconds -- no retuning needed for a longer or shorter descent.
# Several starts are tried because the inverse-response surface has local
# minima (a shallow-dip/slow-lag fit and a deep-dip/fast-lag fit can both be
# locally optimal); the best is kept. Each start is cheap -- the response is
# closed-form, not integrated.
#   a_1 = tau_1 + tau_2   [s]    a_2 = tau_1*tau_2   [s^2]
#   T_z = RHP zero time   [s]    R   = target ramp rate [degC/s]
# a_1/a_2 (the denominator's coefficients) rather than (tau_1, tau_2): the two
# time constants enter the model symmetrically, so fitting them directly has an
# exact exchange symmetry whose Jacobian is rank-deficient along tau_1 == tau_2
# -- which is where the optimum for this data actually sits. The elementary
# symmetric pair quotients that symmetry out, stays smooth through the repeated
# root, AND spans complex poles (a_1^2 < 4*a_2), which the real-(tau_1, tau_2)
# form cannot reach. Bounds are in units of the descent duration / total dT, so
# nothing needs retuning for a different profile.
TOUT_FIT_BOUNDS_REL = {              # (lo, hi) x descent_duration**power
    "a_1": (5e-3, 5.0),              # x D
    "a_2": (2e-5, 25.0),             # x D^2
    "t_z": (2e-3, 3.0),              # x D
    "ramp": (0.3, 20.0),             # x (dT / D)
}
TOUT_FIT_STARTS_REL = (               # (a_1, a_2, T_z, ramp) in those same units
    (0.20, 1.4e-2, 0.16, 1.0),
    (0.20, 1.0e-2, 0.20, 4.6),
    (0.05, 1.0e-3, 0.05, 2.0),
    (0.50, 6.0e-2, 0.50, 1.5),
    (0.10, 2.5e-3, 0.10, 8.0),
    (1.00, 2.5e-1, 0.30, 1.0),
    (0.30, 2.0e-2, 0.60, 3.0),
    (0.08, 1.6e-3, 0.30, 6.0),
    (0.60, 9.0e-2, 0.10, 2.5),
    (0.15, 5.6e-3, 0.40, 12.0),
)


def _tout_nmp_response(t_s: np.ndarray, high_end_s: float, t_high: float, delta_t: float,
                        t_z: float, a_1: float, a_2: float, ramp_rate: float) -> np.ndarray:
    """Exact response of the non-minimum-phase 2nd-order system

        T(s)/U(s) = (1 - T_z*s) / (a_2*s^2 + a_1*s + 1),   T_z > 0

    (a_1 = tau_1 + tau_2, a_2 = tau_1*tau_2 -- so this spans complex poles,
    a_1^2 < 4*a_2, not just two real lags) to U = the piecewise-linear Tout
    target: flat at t_high, a linear ramp at ramp_rate, then flat at
    t_high + delta_t once the ramp has covered delta_t. Evaluated at t_s
    (seconds from window start); high_end_s is where the ramp starts.

    The ramp's DURATION is T_r = delta_t/ramp_rate, a fitted quantity, NOT
    pinned to the descent's length. That matters: pinning it makes the dip
    depth and the terminal error the same number (both ~ramp_rate*(T_z+a_1)),
    so a deep dip forces the response to end that far below t_low -- and the
    measured trace has a deep dip AND arrives on the low hold on time. Letting
    T_r float decouples them. See fit_tout_piecewise_model() for the physical
    caveat this buys.

    Why the RHP (positive) zero: the measured skin temp dips BELOW its own
    cruise level early in the descent, then rises. Driven by this U -- which
    never goes below t_high -- no minimum-phase lag can do that. For a
    first-order lag, T' = (U - T)/tau is >= 0 whenever T <= t_high <= U, so
    T can never leave downward; for a 2nd-order system with no zero and
    T'(0) = 0, T'''(0) = +r/(tau_1*tau_2) > 0, same conclusion. The zero has
    to be in the transfer function. Cross-multiplying the TF gives the ODE

        tau_1*tau_2*T'' + (tau_1 + tau_2)*T' + T = U - T_z*U'

    and it is the -T_z*U' term -- switching on the instant the ramp starts --
    that drives the response the "wrong way" first.

    NOTE on interpretation: the dip physically originates as a FORCING
    excursion (the convection-weighted equilibrium temperature collapsing
    toward the ~-56 degC ambient as air density builds during descent, while
    at cruise thin air and absorbed sunlight hold the skin ~18 degC above
    ambient). Because the piecewise-linear target carries no such excursion,
    that physics gets folded into the response instead, and an RHP zero is
    what it looks like once folded -- which is why "two competing paths with
    opposite signs", K1/(tau_1*s+1) - K2/(tau_2*s+1), is algebraically the
    same model (match coefficients: K1 - K2 = 1, K1*tau_2 - K2*tau_1 = -T_z).
    Consequence for reuse: T_z is PROFILE-SPECIFIC. It stands in for physics
    that scales with descent rate and with the cruise skin-minus-ambient
    offset, NOT with the target's ramp slope, so refit it if either changes
    materially. Generalizing across profiles means moving the dip to the
    input side (build U from altitude + ISA + a solar term) rather than
    leaning harder on T_z.

    Solved in closed form rather than numerically integrated: for U a single
    ramp, U - T_z*U' is linear, so the exact solution is two decaying
    exponentials plus a linear particular solution. The saturating
    ("ramp-and-hold") target is then just SUPERPOSITION of two such ramp
    responses offset by T_r -- one rising ramp minus the same ramp restarted
    at T_r, which cancels the slope and leaves a constant. That is exact at
    any sample spacing (this log's timestamps are irregular), fully
    vectorized, and needs no per-region boundary algebra, which matters
    because the optimizer evaluates it repeatedly over ~32k samples.

    Continuity is automatic everywhere, in value AND slope: the response
    leaves the high hold tangentially (y(0) = y'(0) = 0, matching the hold's
    dT/dt = 0), and as t grows past T_r the superposition tends to
    ramp_rate*T_r = delta_t with zero slope, i.e. it settles onto the low
    hold with no junction step of either kind. (The previous formulation's
    4.15 degC step at that junction is structurally gone, not merely fitted
    away.)

    Numerical safety: poles solve a_2*s^2 + a_1*s + 1 = 0 and both have
    Re(s) < 0 whenever a_1, a_2 > 0 (Routh), so every exponential decays --
    nothing can overflow however far the optimizer wanders. The complex sqrt
    handles real and complex pole pairs in one expression. Only a repeated
    root is degenerate (the A/B solve divides by s1 - s2), so the roots are
    nudged apart by a relative epsilon; unlike the (tau_1, tau_2) form, the
    optimizer has no reason to sit exactly there.

    NOTE the -T_z*U' term's effect is on CURVATURE, not slope: under ramp
    forcing U' steps 0 -> ramp_rate at onset, which puts a step of
    -T_z*ramp_rate on the ODE's right-hand side, and with relative degree 1
    a jump in U^(m) shows up in y^(m+1). So y''(0+) = -T_z*ramp_rate/a_2 < 0
    -- the response leaves the hold tangentially and THEN bends the wrong
    way, which is precisely what preserves dT/dt = 0 while still dipping.
    """
    disc = np.sqrt(np.complex128(a_1 * a_1 - 4.0 * a_2))
    s_1 = (-a_1 + disc) / (2.0 * a_2)
    s_2 = (-a_1 - disc) / (2.0 * a_2)
    if abs(s_1 - s_2) < 1e-9 * abs(s_1):
        s_2 = s_1 * (1.0 + 1e-9)

    lag = a_1 + t_z          # G'(0) = -(a_1 + T_z): the zero adds to apparent lag
    span = ramp_rate * lag
    # y_p = ramp_rate*(t - lag); homogeneous coeffs from y(0) = y'(0) = 0.
    coef_a = (-ramp_rate - s_2 * span) / (s_1 - s_2)
    coef_b = span - coef_a

    def ramp_response(t):
        return (np.real(coef_a * np.exp(s_1 * t) + coef_b * np.exp(s_2 * t))
                + ramp_rate * (t - lag))

    ramp_s = delta_t / ramp_rate                    # T_r: fitted ramp duration
    elapsed = np.maximum(t_s - high_end_s, 0.0)     # 0 while still on the high hold
    y = ramp_response(elapsed)
    past = elapsed > ramp_s
    y[past] -= ramp_response(elapsed[past] - ramp_s)
    return t_high + y


def fit_tout_piecewise_model(df: pd.DataFrame):
    """Model Tout (fuselage skin temp) across the "zones 1-3" window as a
    piecewise-linear target driving a non-minimum-phase 2nd-order response.

    Structure, per the constraints:
      - The INPUT U is the piecewise-linear target: holding_high and
        holding_low are LINEAR with d2T/dt2 = 0 AND dT/dt = 0, i.e. constant,
        joined by a linear ramp across the descent. The least-squares
        constant fit to a segment is just that segment's mean, so each hold
        level IS the measured mean over it -- no optimizer needed for those.
      - The RESPONSE to that input is the 2nd-order system with an RHP zero
        (see _tout_nmp_response()), whose three time parameters T_z, tau_1,
        tau_2 are fit by least squares against the measured descent Tout.

    Fit residuals are taken over the DESCENT only -- that's the segment being
    modeled, and the two hold levels are already pinned to their own means.
    The low hold is then scored separately ("rmse_low_hold_c") as an
    out-of-sample check that the tail settles where it should, rather than
    being folded into the objective.

    Returns a dict with the model Series plus fit metadata, or None if the
    window lacks the columns/segments needed to anchor the fit.
    """
    from scipy.optimize import least_squares

    if "flight_phase" not in df.columns or "tout_c" not in df.columns:
        return None
    phase, tout = df["flight_phase"], df["tout_c"]

    # Same longest-segment identification as main_hold_window() -- a brief
    # pre-takeoff or flickery landing-approach blip can share a phase label
    # with the real holds, so length is what tells them apart.
    change = phase.ne(phase.shift()).cumsum()
    segments = [(seg.iloc[0], seg.index[0], seg.index[-1]) for _, seg in phase.groupby(change)]
    highs = [s for s in segments if s[0] == "holding_high"]
    if not highs:
        return None
    zone1 = max(highs, key=lambda s: s[2] - s[1])
    lows_after = [s for s in segments if s[0] == "holding_low" and s[1] > zone1[2]]
    if not lows_after:
        return None
    zone3 = max(lows_after, key=lambda s: s[2] - s[1])

    high_end, low_start = zone1[2], zone3[1]
    in_zone1 = (df.index >= zone1[1]) & (df.index <= high_end)
    in_zone3 = (df.index >= low_start) & (df.index <= zone3[2])
    in_descent = (df.index > high_end) & (df.index < low_start)

    t_high = tout[in_zone1].mean()   # least-squares constant fit == the mean
    t_low = tout[in_zone3].mean()
    if not np.isfinite(t_high) or not np.isfinite(t_low) or in_descent.sum() < 10:
        return None

    # One shared time axis in seconds, so every region is evaluated by the
    # same closed-form call rather than each region re-deriving its own origin.
    t_s = (df.index - df.index[0]).total_seconds().values
    high_end_s = float((high_end - df.index[0]).total_seconds())
    low_start_s = float((low_start - df.index[0]).total_seconds())
    duration_s = low_start_s - high_end_s
    if duration_s <= 0 or in_descent.sum() < 10:
        return None

    # in_zone1/in_zone3/in_descent are plain numpy bool arrays (DatetimeIndex
    # comparisons, not Series) -- index them directly, no .values.
    measured = tout.values
    fit_mask = in_descent & np.isfinite(measured)
    if fit_mask.sum() < 10:
        return None

    delta_t = float(t_low - t_high)
    # Natural scales for each parameter, so all four are fit as O(1) log
    # multipliers of them -- keeps the Jacobian's columns comparable.
    scales = np.array([duration_s, duration_s ** 2, duration_s, abs(delta_t) / duration_s])

    def evaluate(params, mask=None):
        a_1, a_2, t_z, ramp_rate = np.exp(params) * scales
        target = t_s if mask is None else t_s[mask]
        return _tout_nmp_response(target, high_end_s, t_high, delta_t,
                                   t_z, a_1, a_2, ramp_rate)

    def residual(params):
        return evaluate(params, fit_mask) - measured[fit_mask]

    keys = ("a_1", "a_2", "t_z", "ramp")
    lo = np.log([TOUT_FIT_BOUNDS_REL[k][0] for k in keys])
    hi = np.log([TOUT_FIT_BOUNDS_REL[k][1] for k in keys])
    best = None
    for start in TOUT_FIT_STARTS_REL:
        x0 = np.clip(np.log(np.asarray(start, dtype=float)), lo, hi)
        try:
            trial = least_squares(residual, x0=x0, bounds=(lo, hi),
                                   x_scale="jac", xtol=1e-12, ftol=1e-12)
        except (ValueError, FloatingPointError):
            continue
        if best is None or trial.cost < best.cost:
            best = trial
    if best is None:
        return None

    a_1, a_2, t_z, ramp_rate = (float(v) for v in np.exp(best.x) * scales)
    # Report the poles the fit actually landed on: real pair (two lags) when
    # a_1^2 >= 4*a_2, otherwise a complex pair, which the (a_1, a_2) chart can
    # reach but a real (tau_1, tau_2) parameterization structurally cannot.
    disc = a_1 * a_1 - 4.0 * a_2
    if disc >= 0:
        root = np.sqrt(disc)
        poles = f"real tau={2 * a_2 / (a_1 - root):.0f}s/{2 * a_2 / (a_1 + root):.0f}s"
    else:
        poles = f"complex zeta={a_1 / (2 * np.sqrt(a_2)):.3f} w_n={1 / np.sqrt(a_2) * 1e3:.3f}mrad/s"

    model = pd.Series(evaluate(best.x), index=df.index)

    resid_descent = residual(best.x)
    low_mask = in_zone3 & np.isfinite(measured)
    rmse_low = (float(np.sqrt(np.mean((model.values[low_mask] - measured[low_mask]) ** 2)))
                if low_mask.sum() else float("nan"))

    # Dip metrics: the whole point of the RHP zero is reproducing an excursion
    # BELOW the cruise level, so report modeled vs. measured depth (positive =
    # degC below the high-hold level) to show whether it landed. Two measured
    # figures, because they differ a lot and only one is a fair comparison:
    # the raw minimum is the bottom of a ~+/-4 degC periodic oscillation
    # (visible sawtooth, same period as the POA irradiance ripple), whereas
    # the TREND minimum is the actual excursion a smooth model can represent.
    # Quoting the raw figure alone makes the fit look far worse than it is.
    desc_idx = np.flatnonzero(in_descent)
    trend = tout.rolling("15min", center=True, min_periods=1).mean().values
    modeled_dip = float(t_high - model.values[desc_idx].min())
    measured_dip = float(t_high - np.nanmin(measured[desc_idx]))
    measured_dip_trend = float(t_high - np.nanmin(trend[desc_idx]))
    dip_at_s = float(t_s[desc_idx][int(np.argmin(model.values[desc_idx]))] - high_end_s)
    # Residual against the trend, alongside the oscillation's own RMS about
    # that trend -- together these show how much of rmse_c is irreducible
    # measurement ripple rather than model error.
    rmse_trend = float(np.sqrt(np.nanmean((model.values[desc_idx] - trend[desc_idx]) ** 2)))
    osc_rms = float(np.sqrt(np.nanmean((measured[desc_idx] - trend[desc_idx]) ** 2)))

    return {
        "model": model,
        "t_z": t_z,
        "a_1": a_1,
        "a_2": a_2,
        "poles": poles,
        "ramp_s": delta_t / ramp_rate,
        "t_high_c": float(t_high),
        "t_low_c": float(t_low),
        "descent_s": duration_s,
        "rmse_c": float(np.sqrt(np.mean(resid_descent ** 2))),
        "rmse_low_hold_c": rmse_low,
        "n_descent": int(fit_mask.sum()),
        "modeled_dip_c": modeled_dip,
        "measured_dip_c": measured_dip,
        "measured_dip_trend_c": measured_dip_trend,
        "dip_at_s": dip_at_s,
        "rmse_vs_trend_c": rmse_trend,
        "oscillation_rms_c": osc_rms,
    }


# --------------------------------------------------------------------------
# Main analysis
# --------------------------------------------------------------------------
def analyze(args: argparse.Namespace) -> pd.DataFrame:
    ulog_path = Path(args.ulog)
    if not ulog_path.exists():
        raise FileNotFoundError(f"ULog file not found: {ulog_path}")

    print(f"Loading {ulog_path.name} ...")
    topics = [GPS_TOPIC, *MPPT_TOPICS, TEMP_TOPIC, FLIGHT_TOPIC]
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

    # Attach aircraft attitude, on the same MPPT-anchored timeline, for the
    # panel incidence-angle correction (see cos_incidence_angle()). Skipped
    # entirely under --assume-horizontal, matching the old flat-plate
    # behavior without needing /zeus/flight in the log at all.
    if not args.assume_horizontal:
        attitude_df = load_attitude(ulog)
        if attitude_df is not None:
            tol_attitude = pd.Timedelta(seconds=args.attitude_tolerance_s)
            merged = pd.merge_asof(
                merged.reset_index(), attitude_df.reset_index(),
                on="time", direction="nearest", tolerance=tol_attitude,
            ).set_index("time")
            n_missing_attitude = merged["roll_deg"].isna().sum()
            if n_missing_attitude:
                print(f"  WARNING: {n_missing_attitude} samples ({100 * n_missing_attitude / len(merged):.1f}%) "
                      f"had no attitude match - falling back to a flat horizontal assumption for those rows.")
        else:
            merged["roll_deg"] = np.nan
            merged["pitch_deg"] = np.nan
            merged["yaw_deg"] = np.nan

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

    # Plane-of-array irradiance: the panel isn't flat/horizontal (see
    # PANEL_NORMAL_BODY_STRING_0/_1), so project the direct beam onto the
    # actual, rotating panel normal instead of assuming it always faces
    # straight up. ghi_w_m2 (flat-plate) is kept as-is for comparison;
    # poa_w_m2 is what actually drives pv_power_estimated_bare_w below,
    # unless --assume-horizontal was passed (in which case poa_w_m2 ==
    # ghi_w_m2). poa_w_m2 uses string 1's normal specifically -- it's the
    # array-wide POA the rest of the pipeline (efficiency, main plot) has
    # always used; poa_string0_w_m2/poa_string1_w_m2 below are for the
    # dedicated per-string comparison plot, not for the summary/efficiency
    # numbers. Combining the two strings' actual differing orientations into
    # one array-wide POA would be a real improvement but is out of scope
    # here -- see PANEL_NORMAL_BODY_STRING_0/_1 comment.
    if args.assume_horizontal:
        merged["poa_w_m2"] = merged["ghi_w_m2"]
        merged["poa_string0_w_m2"] = merged["ghi_w_m2"]
        merged["poa_string1_w_m2"] = merged["ghi_w_m2"]
    else:
        roll, pitch, yaw = merged["roll_deg"].values, merged["pitch_deg"].values, merged["yaw_deg"].values
        elevation, azimuth = merged["sun_elevation_deg"].values, merged["sun_azimuth_deg"].values
        cos_aoi_string0 = cos_incidence_angle(roll, pitch, yaw, elevation, azimuth,
                                               normal_body=PANEL_NORMAL_BODY_STRING_0)
        cos_aoi_string1 = cos_incidence_angle(roll, pitch, yaw, elevation, azimuth,
                                               normal_body=PANEL_NORMAL_BODY_STRING_1)
        merged["poa_string0_w_m2"] = merged["dni_w_m2"] * cos_aoi_string0
        merged["poa_string1_w_m2"] = merged["dni_w_m2"] * cos_aoi_string1
        merged["poa_w_m2"] = merged["poa_string1_w_m2"]

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

    # Estimated array output, three tiers (see module docstring for the
    # full breakdown): bare-cell ceiling -> encapsulation loss (ETFE cover x
    # POE encapsulant, both always applied, stacked -- see "Encapsulation
    # transmission") -> optional temp derate. Bare is kept around
    # unconditionally as a reference curve; the encapsulated "nominal" is
    # only kept as its own column when temp derating is on, so the plot can
    # show all three and make the derate's effect visible on top of the
    # (already-applied) encapsulation loss -- when derating is off, nominal
    # IS the final estimated number.
    total_area_m2 = args.cell_count * args.cell_area_cm2 / 1e4
    bare_w = merged["poa_w_m2"] * total_area_m2 * args.cell_efficiency
    merged["pv_power_estimated_bare_w"] = bare_w
    encapsulation_transmission = args.etfe_transmission * args.poe_transmission
    nominal_w = bare_w * encapsulation_transmission
    if args.apply_temp_derate:
        merged["pv_power_estimated_nominal_w"] = nominal_w
        merged["pv_power_estimated_w"] = nominal_w * merged["temp_derate_factor"]
    else:
        merged["pv_power_estimated_w"] = nominal_w

    valid = (merged["sun_elevation_deg"] > 0) & (merged["pv_power_estimated_w"] > 1.0)
    merged["pre_mppt_efficiency_pct"] = np.nan
    merged.loc[valid, "pre_mppt_efficiency_pct"] = (
        100.0 * merged.loc[valid, "pv_power_actual_w"] / merged.loc[valid, "pv_power_estimated_w"]
    )

    # Same three-tier estimate, restricted to MPPT string 1 alone: its own
    # POA (PANEL_NORMAL_BODY_STRING_1) and its own cell count
    # (--string1-cell-count -- NOT half of --cell-count; the two strings
    # aren't equal size). Encapsulation and temp-derate factors are reused
    # as-is since the ETFE/POE cover and Tout apply array-wide, not per
    # string. Only computed if this log actually has a string-1 channel;
    # feeds make_string1_plot() only, not the array-wide summary above.
    if "pv_power_w_1" in merged.columns:
        area_string1_m2 = args.string1_cell_count * args.cell_area_cm2 / 1e4
        bare_string1_w = merged["poa_string1_w_m2"] * area_string1_m2 * args.cell_efficiency
        merged["pv_power_estimated_bare_string1_w"] = bare_string1_w
        nominal_string1_w = bare_string1_w * encapsulation_transmission
        if args.apply_temp_derate:
            merged["pv_power_estimated_nominal_string1_w"] = nominal_string1_w
            merged["pv_power_estimated_string1_w"] = nominal_string1_w * merged["temp_derate_factor"]
        else:
            merged["pv_power_estimated_string1_w"] = nominal_string1_w

        valid_string1 = (merged["sun_elevation_deg"] > 0) & (merged["pv_power_estimated_string1_w"] > 1.0)
        merged["pre_mppt_efficiency_string1_pct"] = np.nan
        merged.loc[valid_string1, "pre_mppt_efficiency_string1_pct"] = (
            100.0 * merged.loc[valid_string1, "pv_power_w_1"]
            / merged.loc[valid_string1, "pv_power_estimated_string1_w"]
        )
        merged.attrs["string1_area_m2"] = area_string1_m2

    # Flight phase, for plot markers -- see classify_flight_phase()/
    # shade_flight_phases(). Altitude-only, so this could run right after
    # the GPS merge; kept here instead so it's next to everything else that
    # only matters for the plots/CSV, not the efficiency calc itself.
    merged["flight_phase"] = classify_flight_phase(merged["alt_msl_m"])

    merged.attrs["total_area_m2"] = total_area_m2
    merged.attrs["valid_mask"] = valid
    return merged


def print_summary(df: pd.DataFrame, args: argparse.Namespace, tz: str) -> None:
    valid = df.attrs["valid_mask"]
    daylight = df[valid]

    dt_hours = df.index.to_series().diff().dt.total_seconds().fillna(0.0) / 3600.0
    dt_hours = dt_hours.clip(upper=10.0 / 3600.0)  # ignore >10s gaps in energy integration
    energy_actual_wh = (df["pv_power_actual_w"].clip(lower=0) * dt_hours).sum()
    energy_estimated_wh = (df["pv_power_estimated_w"].fillna(0) * dt_hours).sum()

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
    print(f"Peak modeled GHI (flat)  : {df['ghi_w_m2'].max():.1f} W/m^2")
    print(f"Peak modeled POA (tilted): {df['poa_w_m2'].max():.1f} W/m^2")
    if args.assume_horizontal:
        print(f"Panel incidence angle    : disabled (--assume-horizontal) -- POA == GHI")
    else:
        # Positive change_pct = tilt GAIN vs flat (banked toward the sun),
        # negative = LOSS (banked away) -- unlike temp derate this has no
        # fixed sign bias, it depends entirely on this flight's maneuvering.
        tilt_valid = df["ghi_w_m2"] > 1.0
        change_pct = 100.0 * (df.loc[tilt_valid, "poa_w_m2"] / df.loc[tilt_valid, "ghi_w_m2"] - 1.0)
        print(f"Panel incidence angle    : mean {change_pct.mean():+.1f}% vs flat-plate GHI  "
              f"(range [{change_pct.min():+.1f}%, {change_pct.max():+.1f}%], roll/pitch/yaw applied)")
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
    print(f"Energy available (model.): {energy_estimated_wh:.1f} Wh")
    if len(daylight):
        eff = daylight["pre_mppt_efficiency_pct"]
        print(f"Pre-MPPT efficiency      : mean {eff.mean():.1f}%  median {eff.median():.1f}%  "
              f"p10-p90 [{eff.quantile(0.1):.1f}%, {eff.quantile(0.9):.1f}%]")
    if energy_estimated_wh > 0:
        print(f"Energy-weighted efficiency: {100 * energy_actual_wh / energy_estimated_wh:.1f}%")
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
    fig, axes = plt.subplots(4, 1, figsize=(12, 13.5), sharex=True)

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

    ax = axes[1]
    ax.plot(df.index, df["ghi_w_m2"], color="tab:orange", alpha=0.6,
            label="Modeled Clear-Sky GHI @ Lat/Lon/Alt (Flat, Horizontal)")
    ax.plot(df.index, df["poa_w_m2"], color="tab:red",
            label="Modeled POA Irradiance (Actual Panel Tilt)")
    ax2 = ax.twinx()
    ax2.plot(df.index, df["sun_elevation_deg"], color="tab:gray", alpha=0.6, label="Sun Elevation")
    ax.set_ylabel("Irradiance (W/m^2)")
    ax2.set_ylabel("Sun elevation (deg)")
    ax.set_title("Modeled Clear-Sky Irradiance Along the Flight Track", fontweight="bold")
    legend_outside(ax, ax2)

    has_derate = "pv_power_estimated_nominal_w" in df.columns

    ax = axes[2]
    ax.plot(df.index, df["pv_power_actual_w"], color="tab:blue", label="Measured Pre-MPPT Power")
    # Bare-cell ceiling, the pre-temp-derate "No Temp Derate" curves, and the
    # temp derate loss/gain shading are all computed and still exported to
    # the CSV (pv_power_estimated_bare_w / _nominal_w), but hidden from the
    # plot itself (2026-08-25) to declutter it now that there are several
    # estimated tiers -- only the final Estimated, Temp-Derated curve (which
    # already has bare/ETFE/POE/tilt all baked in) is drawn against
    # Measured. See print_summary()'s "Temp derate effect" line for the
    # loss/gain magnitude instead of reading it off the plot.
    if has_derate:
        ax.plot(df.index, df["pv_power_estimated_w"], color="tab:olive", linestyle="-.",
                label="Estimated, Temp-Derated, Encapsulated")
    else:
        # pv_power_estimated_w already has the ETFE loss baked in (see
        # analyze()), so it IS the "No Temp Derate, Encapsulated" tier here --
        # there's no separate nominal column to plot on top of the bare-cell line.
        ax.plot(df.index, df["pv_power_estimated_w"], color="tab:green", linestyle="--",
                label="Estimated, No Temp Derate, Encapsulated")
    ax.set_ylabel("Power (W)")
    ax.set_title("Measured vs. Estimated Array Power", fontweight="bold")
    legend_outside(ax)

    ax = axes[3]
    ax.plot(df.index, df["pre_mppt_efficiency_pct"], color="tab:brown",
            label="Measured / Estimated Power")
    ax.axhline(100.0, color="black", linewidth=0.8, linestyle=":", label="100% (measured = estimated)")
    ax.set_ylabel("Efficiency (%)")
    ax.set_title("Measured / Estimated Power", fontweight="bold")
    legend_outside(ax)

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=df.index.tz))
    axes[-1].set_xlabel(f"Local time, {tz} ({df.index[0].date()})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot -> {out_path}")


def rolling_band(series: pd.Series, window_s: float):
    """Rolling mean +/- a rolling standard deviation.

    Returns (center, half_width) so callers can plot center as a line and
    (center - half_width, center + half_width) as a shaded band.
    """
    window = f"{window_s:.0f}s"
    center = series.rolling(window, min_periods=1).mean()
    half_width = series.rolling(window, min_periods=1).std()
    return center, half_width


def plot_pct_diff_panel(ax, df: pd.DataFrame, show_band: bool = True,
                         pct_diff: pd.Series = None, title: str = None) -> None:
    """Draws a String 1 % Difference (measured vs. estimated) panel --
    flight-phase shading, raw trace, 0% reference, and (if show_band)
    a rolling mean/std band -- onto `ax`. Factored out of make_string1_plot()
    so make_string1_plot() (as its bottom panel), make_string1_pct_diff_plot()
    (as its own standalone figure), and make_normal_sweep_plot() (as its top
    and bottom panels) render identically instead of drifting apart.

    show_band=False gives the plain raw trace with no rolling overlay --
    used for the standalone copy, since the decorated version now lives
    alongside the angle-sweep plot instead.

    pct_diff overrides the default df["pre_mppt_efficiency_string1_pct"] - 100
    -- e.g. make_normal_sweep_plot()'s bottom panel passes in
    compute_pct_diff_at_angle()'s result to show the % Difference under the
    empirically optimal angle instead of the surveyed one. title likewise
    overrides the default title, so that panel can say which angle it used.

    Does not add a legend -- callers place that differently (inside vs.
    outside the axes, extra handles or not), so each calls its own
    ax.legend(...)/legend_outside(...) afterward.
    """
    shade_flight_phases([ax], df["flight_phase"])
    if pct_diff is None:
        # % difference from estimated (measured/estimated - 1) x 100, not
        # the raw ratio -- 0% = measured matches estimated, +/- reads
        # directly as over/under, rather than needing to mentally subtract
        # 100 each time.
        pct_diff = df["pre_mppt_efficiency_string1_pct"] - 100.0
    if show_band:
        window_min = PCT_DIFF_ROLLING_WINDOW_S / 60.0
        center, half_width = rolling_band(pct_diff, PCT_DIFF_ROLLING_WINDOW_S)
        ax.fill_between(df.index, center - half_width, center + half_width,
                         color="tab:brown", alpha=0.15, zorder=1,
                         label=f"{window_min:.0f}-min Rolling Std")
        ax.plot(df.index, pct_diff, color="tab:brown", alpha=0.5, linewidth=0.8,
                label="String 1 % Difference (Measured vs. Estimated)")
        ax.plot(df.index, center, color="black", linewidth=1.2,
                label=f"{window_min:.0f}-min Rolling Mean")
    else:
        ax.plot(df.index, pct_diff, color="tab:brown",
                label="String 1 % Difference (Measured vs. Estimated)")
    ax.axhline(0.0, color="black", linewidth=0.8, linestyle=":", label="0% (measured = estimated)")
    ax.set_ylabel("% Difference")
    ax.set_title(title or "String 1: % Difference, Measured vs. Estimated Power", fontweight="bold")


def make_string1_plot(df: pd.DataFrame, out_path: Path, tz: str) -> None:
    """The same environmentals/irradiance/power/efficiency workflow as
    make_plot(), applied to MPPT string 1 alone instead of the whole array:
    its own POA (poa_string1_w_m2, PANEL_NORMAL_BODY_STRING_1), its own
    measured channel (pv_power_w_1), and its own cell count
    (--string1-cell-count). Environmentals (altitude/Tout) isn't actually
    string-specific -- it's array-wide, same as in make_plot() -- but is
    repeated here so this plot stands alone. Only the final estimated tier
    is drawn against measured, matching make_plot()'s decluttered convention.

    Only covers "zones 1-3" -- the main holding_high / descending /
    holding_low sequence (see main_hold_window()) -- dropping the climb-out
    before it and the final approach/landing after it. One contiguous
    real-time slice (see keep_main_hold()), so the x-axis is plain real
    clock time throughout, same as make_plot().
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    local_index = df.index.tz_convert(tz)
    df = df.set_axis(local_index)
    df = keep_main_hold(df)

    has_tout = "tout_c" in df.columns and df["tout_c"].notna().any()
    fig, axes = plt.subplots(4, 1, figsize=(12, 13.5), sharex=True)
    # Flight-phase shading only goes on the Efficiency panel (axes[3]) here,
    # not Environmentals/Irradiance/Power -- see the call just before that
    # panel's plotting below.

    def legend_outside(ax, *extra_axes):
        handles, labels = ax.get_legend_handles_labels()
        for extra in extra_axes:
            h, l = extra.get_legend_handles_labels()
            handles += h
            labels += l
        ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.06, 1.0), borderaxespad=0.0)

    ax = axes[0]
    ax.plot(df.index, df["alt_msl_m"], color="tab:purple", label="Altitude (MSL)")
    ax.set_ylabel("Altitude MSL (m)")
    if has_tout:
        ax2 = ax.twinx()
        ax2.plot(df.index, df["tout_c"], color="tab:brown", label="Fuselage Skin Temp (Tout Proxy)")
        # Piecewise-linear target -> non-minimum-phase (RHP-zero) 2nd-order
        # response, fit to this window's own Tout -- see
        # fit_tout_piecewise_model() / _tout_nmp_response().
        tout_fit = fit_tout_piecewise_model(df)
        if tout_fit is not None:
            ax2.plot(df.index, tout_fit["model"], color="tab:blue", linewidth=1.8,
                     label=f"Tout 2nd-Order Model, RHP Zero "
                           f"(T_z={tout_fit['t_z']:.0f}s, RMSE {tout_fit['rmse_c']:.2f} degC)")
            print(f"  Tout model: PWL target {tout_fit['t_high_c']:.1f} -> {tout_fit['t_low_c']:.1f} degC "
                  f"into a 2nd-order response with an RHP zero; descent fit over "
                  f"{tout_fit['descent_s'] / 60.0:.0f} min, {tout_fit['n_descent']} samples")
            print(f"    T_z={tout_fit['t_z']:.0f}s  poles: {tout_fit['poles']}  "
                  f"target ramp T_r={tout_fit['ramp_s'] / 60.0:.1f} min")
            print(f"    descent RMSE {tout_fit['rmse_c']:.2f} degC vs raw / "
                  f"{tout_fit['rmse_vs_trend_c']:.2f} degC vs 15-min trend "
                  f"(measurement ripple about that trend is {tout_fit['oscillation_rms_c']:.2f} degC, "
                  f"irreducible); low-hold RMSE {tout_fit['rmse_low_hold_c']:.2f} degC "
                  f"(not fit, out-of-sample)")
            print(f"    inverse-response dip below cruise: modeled "
                  f"{tout_fit['modeled_dip_c']:.1f} degC at +{tout_fit['dip_at_s'] / 60.0:.0f} min, "
                  f"measured {tout_fit['measured_dip_trend_c']:.1f} degC on the trend "
                  f"({tout_fit['measured_dip_c']:.1f} degC at the raw oscillation trough)")
        ax2.axhline(STC_TEMP_C, color="black", linewidth=0.8, linestyle=":", label=f"STC ({STC_TEMP_C:.0f} degC)")
        ax2.set_ylabel("Temp (degC)")
        ax.set_title("Environmentals", fontweight="bold")
        legend_outside(ax, ax2)
    else:
        ax.set_title("Environmentals", fontweight="bold")
        legend_outside(ax)

    ax = axes[1]
    ax.plot(df.index, df["ghi_w_m2"], color="tab:orange", alpha=0.6,
            label="Modeled Clear-Sky GHI @ Lat/Lon/Alt (Flat, Horizontal)")
    ax.plot(df.index, df["poa_string1_w_m2"], color="tab:red", label="String 1 POA Irradiance")
    ax2 = ax.twinx()
    ax2.plot(df.index, df["sun_elevation_deg"], color="tab:gray", alpha=0.6, label="Sun Elevation")
    ax.set_ylabel("Irradiance (W/m^2)")
    ax2.set_ylabel("Sun elevation (deg)")
    ax.set_title("String 1: Modeled Clear-Sky Irradiance", fontweight="bold")
    legend_outside(ax, ax2)

    ax = axes[2]
    ax.plot(df.index, df["pv_power_w_1"], color="tab:blue", label="String 1 Measured Power")
    ax.plot(df.index, df["pv_power_estimated_string1_w"], color="tab:olive", linestyle="-.",
            label="String 1 Estimated Power")
    ax.set_ylabel("Power (W)")
    ax.set_title("String 1: Measured vs. Estimated Power", fontweight="bold")
    legend_outside(ax)

    ax = axes[3]
    plot_pct_diff_panel(ax, df, show_band=False)
    legend_outside(ax)

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=df.index.tz))
    axes[-1].set_xlabel(f"Local time, {tz} ({df.index[0].date()})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot -> {out_path}")


def make_string1_pct_diff_plot(df: pd.DataFrame, out_path: Path, tz: str) -> None:
    """Standalone copy of make_string1_plot()'s bottom panel -- % Difference
    (measured vs. estimated), String 1 -- as its own single-panel figure.
    Same data and "zones 1-3" window (keep_main_hold()) as the full 4-panel
    plot, but plain (no rolling mean/std band, see plot_pct_diff_panel()) --
    the decorated version now lives in make_normal_sweep_plot()'s top panel
    instead, so this one stays simple.
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    local_index = df.index.tz_convert(tz)
    df = df.set_axis(local_index)
    df = keep_main_hold(df)

    fig, ax = plt.subplots(figsize=(12, 4.5))
    plot_pct_diff_panel(ax, df, show_band=False)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=df.index.tz))
    ax.set_xlabel(f"Local time, {tz} ({df.index[0].date()})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot -> {out_path}")


def make_strings_plot(df: pd.DataFrame, out_path: Path, tz: str) -> None:
    """Per-MPPT-string voltage and current, one line per string with a
    consistent color across both panels. Separate from make_plot()'s summed
    pv_power_actual_w, which hides channel-specific behavior -- e.g. one
    string sagging in voltage, or dropping out, while the total looks fine.

    String IDs are discovered from the dataframe's own columns (pv_voltage_
    v_<id>) rather than hardcoded, so this naturally adapts to however many
    /zeus/mppt_* channels were actually present in this log (see load_mppt()
    -- a missing channel is skipped there, not backfilled with a placeholder).
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    local_index = df.index.tz_convert(tz)
    df = df.set_axis(local_index)

    string_ids = sorted(
        c[len("pv_voltage_v_"):] for c in df.columns if c.startswith("pv_voltage_v_")
    )
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, axes = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True)

    def legend_outside(ax):
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

    ax = axes[0]
    for i, sid in enumerate(string_ids):
        ax.plot(df.index, df[f"pv_voltage_v_{sid}"], color=colors[i % len(colors)],
                label=f"String {sid} Voltage")
    ax.axhline(MPPT_12SBB_VIN_MIN_V, color="black", linewidth=0.8, linestyle="--",
               label=f"MPPT-12SBB Input Range ({MPPT_12SBB_VIN_MIN_V:.0f}-{MPPT_12SBB_VIN_MAX_V:.0f}V)")
    ax.axhline(MPPT_12SBB_VIN_MAX_V, color="black", linewidth=0.8, linestyle="--")
    ax.set_ylabel("Voltage (V)")
    ax.set_title("Per-String PV Voltage", fontweight="bold")
    legend_outside(ax)

    ax = axes[1]
    for i, sid in enumerate(string_ids):
        ax.plot(df.index, df[f"pv_current_a_{sid}"], color=colors[i % len(colors)],
                label=f"String {sid} Current")
    ax.set_ylabel("Current (A)")
    ax.set_title("Per-String PV Current", fontweight="bold")
    legend_outside(ax)

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=df.index.tz))
    axes[-1].set_xlabel(f"Local time, {tz} ({df.index[0].date()})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot -> {out_path}")


def make_poa_strings_plot(df: pd.DataFrame, out_path: Path, tz: str) -> None:
    """Modeled POA irradiance per string (poa_string0_w_m2 vs poa_string1_w_m2)
    -- each string has its own surface normal (PANEL_NORMAL_BODY_STRING_0/_1),
    so the two can genuinely diverge whenever the aircraft's attitude favors
    one string's mounting angle over the other's, even though both see the
    same clear-sky DNI. This is the modeled/estimated-side counterpart to
    make_strings_plot()'s measured voltage/current -- separate plots because
    one is W/m^2 (irradiance model) and the other is V/A (electrical
    measurement); nothing here is compared against pv_power_w_<id> directly.
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    local_index = df.index.tz_convert(tz)
    df = df.set_axis(local_index)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(df.index, df["poa_string0_w_m2"], color="tab:blue", label="String 0 POA Irradiance")
    ax.plot(df.index, df["poa_string1_w_m2"], color="tab:orange", label="String 1 POA Irradiance")
    ax.set_ylabel("POA Irradiance (W/m^2)")
    ax.set_title("Modeled Plane-of-Array Irradiance, String 0 vs String 1", fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=df.index.tz))
    ax.set_xlabel(f"Local time, {tz} ({df.index[0].date()})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot -> {out_path}")


def compute_pct_diff_at_angle(window: pd.DataFrame, args, theta_deg: float) -> pd.Series:
    """String 1's % Difference (measured vs. estimated, already offset -100)
    recomputed with the panel normal fixed at theta_deg -- same x_cad(theta) =
    -cos(theta), z_cad(theta) = sin(theta) parameterization as
    sweep_panel_normal_angle() -- instead of the surveyed PANEL_NORMAL_BODY_
    STRING_1. Factored out of that function's inner loop so the identical
    per-angle calculation can also drive a real time-series plot (see
    make_normal_sweep_plot()'s bottom panel) instead of only a reduced
    summary statistic.

    `window` must already be restricted to whatever time range the caller
    wants (e.g. keep_main_hold()'s "zones 1-3") -- this function does no
    windowing of its own, so the same `window` can be reused across many
    calls without recomputing it each time.
    """
    roll = window["roll_deg"].values
    pitch = window["pitch_deg"].values
    yaw = window["yaw_deg"].values
    elevation = window["sun_elevation_deg"].values
    azimuth = window["sun_azimuth_deg"].values
    dni = window["dni_w_m2"].values
    temp_derate_factor = window["temp_derate_factor"].values
    measured = window["pv_power_w_1"]
    sun_up = window["sun_elevation_deg"] > 0

    area_string1_m2 = args.string1_cell_count * args.cell_area_cm2 / 1e4
    encapsulation_transmission = args.etfe_transmission * args.poe_transmission

    theta = np.radians(theta_deg)
    normal_body = _panel_normal_body(x_cad=-np.cos(theta), z_cad=np.sin(theta))
    cos_aoi = cos_incidence_angle(roll, pitch, yaw, elevation, azimuth, normal_body=normal_body)
    estimated_w = dni * cos_aoi * area_string1_m2 * args.cell_efficiency * encapsulation_transmission
    estimated_w = estimated_w * temp_derate_factor
    estimated_w = pd.Series(estimated_w, index=window.index)

    valid = sun_up & (estimated_w > 1.0)
    pct_diff = pd.Series(np.nan, index=window.index)
    pct_diff.loc[valid] = 100.0 * measured.loc[valid] / estimated_w.loc[valid] - 100.0
    return pct_diff


def sweep_panel_normal_angle(df: pd.DataFrame, args, angle_min_deg: float = 0.0,
                              angle_max_deg: float = 180.0, angle_step_deg: float = 1.0) -> pd.DataFrame:
    """Sweep MPPT string 1's panel-normal unit vector through the aircraft's
    X-Z (longitudinal/vertical, CAD-frame) plane and recompute the % Difference
    rolling-std band at each angle -- a sensitivity check on the surveyed
    normal (PANEL_NORMAL_BODY_STRING_1), rather than a fixed assumption.

    Parameterized as x_cad(theta) = -cos(theta), z_cad(theta) = sin(theta),
    which traces the unit circle through:
        theta =   0 deg -> (x=-1, z= 0)  pointing toward the nose, flat
        theta =  90 deg -> (x= 0, z= 1)  straight up
        theta = 180 deg -> (x=+1, z= 0)  pointing toward the tail, flat
    (angle_min_deg/angle_max_deg need not span the full 0-180 -- e.g. a
    narrow, fine-stepped range around a known minimum). Each (x_cad, z_cad)
    is converted to a body-frame normal via _panel_normal_body() -- the same
    CAD->body transform PANEL_NORMAL_BODY_STRING_0/_1 themselves use -- so
    this sweeps the same kind of vector, just varied instead of fixed.

    Restricted to the same "zones 1-3" window as make_string1_plot() (see
    keep_main_hold()), and always applies the tilt correction regardless of
    --assume-horizontal -- the whole point here is comparing different tilts.

    Returns a DataFrame indexed by angle_deg with one column,
    mean_rolling_std: the flight-mean of the % Difference rolling std
    (rolling_band(), same PCT_DIFF_ROLLING_WINDOW_S window) at that angle.
    NaN at angles where the panel would face away from the sun for the
    entire window (no valid estimate to compare against).
    """
    window = keep_main_hold(df)
    angles = np.arange(angle_min_deg, angle_max_deg + angle_step_deg / 2.0, angle_step_deg)
    mean_rolling_std = []
    for theta_deg in angles:
        pct_diff = compute_pct_diff_at_angle(window, args, theta_deg)
        _, half_width = rolling_band(pct_diff, PCT_DIFF_ROLLING_WINDOW_S)
        mean_rolling_std.append(half_width.mean())  # pandas .mean() skips NaN; all-NaN -> NaN, no warning

    return pd.DataFrame({"angle_deg": angles, "mean_rolling_std": mean_rolling_std}).set_index("angle_deg")


def make_normal_sweep_plot(df: pd.DataFrame, sweep: pd.DataFrame, out_path: Path, tz: str, args) -> None:
    """Three-panel figure:
      1. String 1's % Difference (measured vs. estimated) under the ASSUMED
         normal (PANEL_NORMAL_BODY_STRING_1), with its rolling mean/std band
         (plot_pct_diff_panel(show_band=True)).
      2. sweep_panel_normal_angle()'s result -- mean rolling std vs. assumed
         panel-normal angle. A dip identifies the angle that minimizes
         measured-vs-estimated scatter; both the assumed normal and the
         sweep's own minimum are marked (angle + std value in the legend).
      3. The same % Difference as panel 1, but recomputed under the sweep's
         OPTIMAL angle (compute_pct_diff_at_angle()) -- a direct visual
         before/after of what that angle actually buys.
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    local_index = df.index.tz_convert(tz)
    ts = df.set_axis(local_index)
    ts = keep_main_hold(ts)

    # Invert _panel_normal_body()'s (x_cad, z_cad) -> (nx, ny, nz) = normalize(
    # [-x_cad, 0, -z_cad]) to recover the angle this sweep would assign the
    # assumed String 1 normal: cos(theta) = nx, sin(theta) = -nz (the
    # normalization factor cancels out of the ratio, so it doesn't matter
    # that PANEL_NORMAL_BODY_STRING_1's inputs weren't exactly unit length).
    nx, _, nz = PANEL_NORMAL_BODY_STRING_1
    assumed_theta = np.degrees(np.arctan2(-nz, nx))

    valid = sweep["mean_rolling_std"].dropna()
    min_theta = valid.idxmin() if not valid.empty else None
    min_value = valid.min() if not valid.empty else None

    # Value of the curve at the assumed angle itself (nearest grid point --
    # assumed_theta generally doesn't land exactly on a sweep step), so the
    # assumed normal's own scatter can be compared numerically against the
    # minimum, not just located on the x-axis via the vline.
    nearest_pos = sweep.index.get_indexer([assumed_theta], method="nearest")[0]
    assumed_grid_theta = sweep.index[nearest_pos]
    assumed_value = sweep["mean_rolling_std"].iloc[nearest_pos]
    if pd.isna(assumed_value):
        assumed_value = None

    fig, (ax_ts, ax_sweep, ax_opt) = plt.subplots(3, 1, figsize=(11, 14.5))

    plot_pct_diff_panel(ax_ts, ts, show_band=True,
                        title=f"String 1: % Difference, Measured vs. Estimated Power "
                              f"(Assumed Normal, {assumed_theta:.2f} deg)")
    ax_ts.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    ax_ts.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=ts.index.tz))
    ax_ts.set_xlabel(f"Local time, {tz} ({ts.index[0].date()})")

    ax_sweep.plot(sweep.index, sweep["mean_rolling_std"], color="tab:brown", marker=".", markersize=3)
    ax_sweep.axvline(assumed_theta, color="black", linewidth=0.8, linestyle="--",
                      label="Assumed String 1 Normal")
    if assumed_value is not None:
        ax_sweep.plot(assumed_grid_theta, assumed_value, marker="o", color="tab:blue", markersize=7,
                       zorder=5,
                       label=f"Assumed Normal Std: {assumed_value:.2f} at {assumed_theta:.2f} deg")
    if min_theta is not None:
        ax_sweep.plot(min_theta, min_value, marker="o", color="tab:red", markersize=7, zorder=5,
                       label=f"Minimum Std: {min_value:.2f} at {min_theta:.2f} deg")
    ax_sweep.set_xlabel("Assumed Panel-Normal Angle (deg): 0 = -X (nose), 90 = +Z (up), 180 = +X (tail)")
    ax_sweep.set_ylabel(f"Mean {PCT_DIFF_ROLLING_WINDOW_S / 60.0:.0f}-min Rolling Std of % Difference")
    ax_sweep.set_title("% Difference Rolling Deviation vs Panel-Normal Angle Sweep", fontweight="bold")
    ax_sweep.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

    if min_theta is not None:
        optimized_pct_diff = compute_pct_diff_at_angle(ts, args, min_theta)
        plot_pct_diff_panel(ax_opt, ts, show_band=True, pct_diff=optimized_pct_diff,
                             title=f"String 1: % Difference, Measured vs. Estimated Power "
                                   f"(Optimized Normal, {min_theta:.2f} deg)")
        ax_opt.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
        ax_opt.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=ts.index.tz))
        ax_opt.set_xlabel(f"Local time, {tz} ({ts.index[0].date()})")

        # Same y-axis on both % Difference panels (assumed vs. optimized) so
        # they're directly comparable at a glance -- otherwise each
        # autoscales to its own outlier spikes and a real difference between
        # the two can be masked by a mere axis-scale difference.
        y_min = min(ax_ts.get_ylim()[0], ax_opt.get_ylim()[0])
        y_max = max(ax_ts.get_ylim()[1], ax_opt.get_ylim()[1])
        ax_ts.set_ylim(y_min, y_max)
        ax_opt.set_ylim(y_min, y_max)
    else:
        ax_opt.axis("off")

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
    parser.add_argument("--string1-cell-count", type=int, default=DEFAULT_STRING1_CELL_COUNT,
                         help="Cells wired into MPPT string 1 specifically (NOT half of "
                              "--cell-count -- the two strings aren't necessarily equal size). "
                              "Used only for the dedicated string-1 efficiency plot.")
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
    parser.add_argument("--attitude-tolerance-s", type=float, default=0.5,
                         help="Max time gap allowed when matching a roll/pitch/yaw sample "
                              "(/zeus/flight) to an MPPT sample")
    parser.add_argument("--apply-temp-derate", action="store_true",
                         help="Derate estimated power using the datasheet's power temp "
                              "coefficient and the Tout proxy (fuselage skin temp). Off by "
                              "default -- see docstring for why this is opt-in.")
    parser.add_argument("--assume-horizontal", action="store_true",
                         help="Ignore aircraft attitude and assume the panel always faces "
                              "straight up (the old behavior), instead of projecting irradiance "
                              "onto the actual, rotating panel normal via /zeus/flight's "
                              "roll/pitch/yaw. Panel incidence-angle correction is ON by "
                              "default -- see docstring for why this is opt-OUT.")
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
        "lat", "lon", "alt_msl_m", "flight_phase", "sun_elevation_deg", "sun_azimuth_deg",
        "ghi_w_m2", "dni_w_m2", "dhi_w_m2", "poa_w_m2", "poa_string0_w_m2", "poa_string1_w_m2",
        *(["roll_deg", "pitch_deg", "yaw_deg"] if "roll_deg" in df.columns else []),
        *(["tout_c"] if "tout_c" in df.columns else []),
        *[c for c in df.columns if c.startswith("pv_")],
        *(["temp_derate_factor"] if "temp_derate_factor" in df.columns else []),
        "pre_mppt_efficiency_pct",
        *(["pre_mppt_efficiency_string1_pct"] if "pre_mppt_efficiency_string1_pct" in df.columns else []),
    ]
    df[export_cols].to_csv(csv_path)
    print(f"Saved data -> {csv_path}")

    if not args.no_plot:
        plot_path = out_dir / f"{stem}_solar_efficiency.png"
        make_plot(df, plot_path, tz)
        if not args.no_open:
            open_in_vscode(plot_path)

        strings_plot_path = out_dir / f"{stem}_strings.png"
        make_strings_plot(df, strings_plot_path, tz)
        if not args.no_open:
            open_in_vscode(strings_plot_path)

        poa_strings_plot_path = out_dir / f"{stem}_poa_strings.png"
        make_poa_strings_plot(df, poa_strings_plot_path, tz)
        if not args.no_open:
            open_in_vscode(poa_strings_plot_path)

        if "pv_power_w_1" in df.columns:
            string1_plot_path = out_dir / f"{stem}_string1.png"
            make_string1_plot(df, string1_plot_path, tz)
            if not args.no_open:
                open_in_vscode(string1_plot_path)

            pct_diff_plot_path = out_dir / f"{stem}_string1_pct_diff.png"
            make_string1_pct_diff_plot(df, pct_diff_plot_path, tz)
            if not args.no_open:
                open_in_vscode(pct_diff_plot_path)

            print("Sweeping String 1 panel-normal angle (90-105 deg, 0.05 deg steps) ...")
            sweep = sweep_panel_normal_angle(df, args, angle_min_deg=90.0, angle_max_deg=105.0,
                                              angle_step_deg=0.05)
            sweep_plot_path = out_dir / f"{stem}_normal_sweep.png"
            make_normal_sweep_plot(df, sweep, sweep_plot_path, tz, args)
            if not args.no_open:
                open_in_vscode(sweep_plot_path)
        else:
            print("  (note: no pv_power_w_1 in this log - skipping the string-1 plots)")


if __name__ == "__main__":
    main()
