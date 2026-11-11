import os
import math
import numpy as np
import ephem
from astropy import units as u
from astropy.time import Time

from src.txt_reader import *


# PARAMETERS
tles_dir = os.path.join(os.path.dirname(__file__),"src", "tles")
gs_dir = os.path.dirname(os.path.dirname(__file__))

# Simulation parameters
# Used by BOTH Hypatia (orbit propagation / dynamic state) and ns-3 (simulation_end_time_ns,
# dynamic_state_update_interval_ns): main.py writes them into config_ns3.properties.
# One fstate_<t>.txt and one gsl_if_bandwidth_<t>.txt are written per time step
# -> number of files = duration / time_step.
# 3600 s covers the first LEO-GEO window of both Sentinels (S1-1: 00:00-00:23, S1-0: 00:34-01:12 UTC)
# and the hand-over between them (see validation/results/report.md).
class simulation:
    duration = 3600          # s
    time_step = 1000000000   # ns (1 s)
    # The run always starts at conf.epoch (ns-3 starts at the TLE epoch, so no offset is supported).
    # To simulate a later window, move conf.epoch (the TLEs are regenerated with it).

# Satellites to simulate (edit these lists to add/remove satellites)
# Satellite index in Hypatia/ns-3 = order of appearance: LEO groups first, then GEO.
class constellations:
    # LEO: (name, inclination [deg], mean motion [rev/day], mean anomaly of each satellite [deg])
    leo = [
        ("Sentinel-1", 98.18, 14.60, [0, 180]),   # sats 0, 1
        ("Sentinel-2", 98.62, 14.31, [0, 180]),   # sats 2, 3
        ("Sentinel-3", 98.65, 14.26, [0, 140]),   # sats 4, 5  (S3A/S3B are 140 deg apart)
    ]
    # GEO: (name, longitude [deg East])
    geo = [
        ("EDRS-A", 9.0),                          # sat 6
        ("EDRS-C", 31.0),                         # sat 7
    ]

# Geometry
class geometric:
    Re = 6378.137         # km
    ATM = 100.0           # km
    leo_sat = sum(len(c[3]) for c in constellations.leo)   # number of LEO satellites
    geo_sat = len(constellations.geo)                      # number of GEO satellites
    sat_tot = leo_sat + geo_sat
    gs_tot = get_num_ground_stations(os.path.join(gs_dir, "ground_stations", "ground_stations.txt"))

# General configuration
class conf :
    epoch = Time("2026-01-01 00:00:00", scale="utc")

# Communication parameters
# NOTE on the Hypatia/ns-3 model: every node has ONE GSL interface and its rate is that node's
# TRANSMIT rate (satellite -> downlink, ground station -> uplink). ISLs share one rate.
class communication: # datasheet parameters
    optical_oisl_bps = 1.8e9     # LEO -> GEO: EDRS Laser Communication Terminal, 1.8 Gbit/s
                                 #   (Sentinel-1 uses it for its 520 Mbit/s payload stream)
    ka_geo_gs_down_bps = 1.8e9   # GEO -> GS: EDRS Ka-band feeder link, 1.8 Gbit/s
                                 #   (ESA ARTES EDRS: "data rate requirement for the download of 1.8 Gbit/s").
                                 #   The 300 Mbit/s figure is the EDRS-A Ka-band INTER-SATELLITE link, not the feeder.
    direct_leo_gs_bps  = 520e6   # LEO -> GS: Sentinel-1 X-band payload downlink, 2 x 260 Mbit/s
    gs_tc_uplink_bps   = 64e3    # GS -> LEO: Sentinel-1 S-band telecommand uplink (reference only)
    gs_source_bps      = 520e6   # GS transmit rate used in ns-3 (see ground_station_gsl_rate_mbps)

    isl_rate_mbps              = optical_oisl_bps / 1_000_000   # 1800.0 Mbps
    gsl_base_rate_mbps         = 10000.0                        # 10 Gbps base channel scale
    leo_gs_rate_mbps           = direct_leo_gs_bps / 1_000_000  # 520.0 Mbps
    geo_gs_downlink_rate_mbps  = ka_geo_gs_down_bps / 1_000_000 # 1800.0 Mbps
    # Ground-station transmit rate. In Hypatia, traffic (TCP flows, pings) can only start and end
    # at ground stations, so the source ground station stands in for the Sentinel payload:
    # its rate is set to the Sentinel payload rate (520 Mbit/s) so the injected traffic matches what
    # a Sentinel actually delivers. The physical uplink is only 64 kbit/s (S-band TC); using that
    # value would cap every flow at 64 kbit/s.
    ground_station_gsl_rate_mbps = gs_source_bps / 1_000_000    # 520.0 Mbps

# Visibility
class visibility: # datasheet parameters
    theta_geo_max = 90         # º  GEO terminal half-cone around nadir
    theta_leo_max = 90         # º  LEO terminal half-cone around zenith
    max_range = 800000         # km (effectively no limit; LCT tested to 75,000 km)
    min_elevation_gs = 5.0     # º  ground station elevation mask (satellite must be above it)
    max_gsl_length = 45000000  # m  upper bound on GS-satellite slant range (GEO from Europe ~ 38,000-40,000 km)

# ns-3 traffic: pings between ground stations to measure latency / outages through the relay.
# Node ids in ns-3: satellites 0..n_sat-1, then ground stations (n_sat + gid).
# Ground stations: 0 Redu, 1 Weilheim, 2 Matera, 3 Hawaii  ->  ns-3 ids 8, 9, 10, 11 (with 8 satellites).
class ns3:
    enable_pingmesh = True
    pingmesh_interval_ns = 100000000            # one ping every 100 ms per pair
    pingmesh_pairs_gid = [(3, 1), (2, 1)]       # Hawaii -> Weilheim, Matera -> Weilheim (ground-station ids)


# BASIC CONFIGURATION FUNCTIONS


# set/overwrite keys in config_ns3.properties (data rates, simulation time, ...)
def update_ns3_properties(config_path, rate_properties):
    with open(config_path, "r") as f_in:
        lines = f_in.readlines()

    updated_keys = set()
    for line_index, line in enumerate(lines):
        key, separator, _ = line.partition("=")
        key = key.strip()
        if separator and key in rate_properties:
            lines[line_index] = f"{key}={rate_properties[key]}\n"
            updated_keys.add(key)

    for key, value in rate_properties.items():
        if key not in updated_keys:
            lines.append(f"{key}={value}\n")

    with open(config_path, "w") as f_out:
        f_out.writelines(lines)


update_ns3_rate_properties = update_ns3_properties  # old name

# Orbital radius
def radius(sat, epoch_str, date_str):
    
    # Returns the GEOCENTRIC position vector in METRES (equatorial frame of date).
    #
    # NOTE: an ephem.Observer at lat=0, lon=0, elevation=0 is on the Earth's SURFACE
    # (Gulf of Guinea), not at its centre, so sat.range is a topocentric slant range.
    # Geocentric direction -> g_ra / g_dec; geocentric distance -> rebuilt from the
    # sub-satellite point and the height above the WGS-84 ellipsoid (sat.elevation).
    from satgen.distance_tools.distance_tools import geodetic2cartesian

    try:
        sat.compute(str(date_str).replace('-', '/'), epoch=str(epoch_str).replace('-', '/'))
        ra = float(sat.g_ra)
        dec = float(sat.g_dec)
        r = float(np.linalg.norm(geodetic2cartesian(
            math.degrees(sat.sublat), math.degrees(sat.sublong), sat.elevation)))
    except RuntimeError as e:
        print(f"Error en satélite: {sat.name}")
        print(f"Época interna TLE del satélite (sat._epoch): {getattr(sat, '_epoch', 'N/A')}")
        raise e

    # 4. Spherical to Cartesian ECI conversion:
    # x = r * cos(dec) * cos(ra)
    # y = r * cos(dec) * sin(ra)
    # z = r * sin(dec)
    x = r * math.cos(dec) * math.cos(ra)
    y = r * math.cos(dec) * math.sin(ra)
    z = r * math.sin(dec)

    return np.array([x, y, z])

