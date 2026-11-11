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
class simulation:
    duration = 200           # s
    time_step = 1000000000   # ns

# Geometry
class geometric:
    Re = 6378.137         # km
    ATM = 100.0           # km
    _, leo_sat = analyze_tle_file(os.path.join(tles_dir, "sentinel1_tle.txt"))
    leo_sat = int(leo_sat)
    _, geo_sat = analyze_tle_file(os.path.join(tles_dir, "edcrs_tle.txt"))
    geo_sat = int(geo_sat)
    sat_tot = leo_sat + geo_sat    
    gs_tot = get_num_ground_stations(os.path.join(gs_dir, "ground_stations", "ground_stations.txt"))

# General configuration
class conf :
    epoch = Time("2026-01-01 00:00:00", scale="utc")

# Communication parameters
class communication: # datasheet parameters
    optical_oisl_bps = 1.8e9   # LEO -> GEO: Laser Inter-Satellite Link (1.8 Gbps)
    ka_geo_gs_down_bps = 3e8 # GEO -> GS: Ka-band Feeder Downlink (300 Mbps)
    direct_leo_gs_bps  = 520e6 # LEO -> GS: Direct X-band Downlink (520 Mbps)
    
    isl_rate_mbps              = optical_oisl_bps / 1_000_000  # 1800.0 Mbps
    gsl_base_rate_mbps         = 10000.0                       # 10 Gbps base channel scale
    leo_gs_rate_mbps           = direct_leo_gs_bps / 1_000_000 # 520.0 Mbps
    geo_gs_downlink_rate_mbps  = ka_geo_gs_down_bps / 1_000_000# 300.0 Mbps

# Visibility
class visibility: # datasheet parameters
    theta_geo_max = 90         # º 
    theta_leo_max = 90         # º
    max_range = 800000         # km
    max_gsl_length = 45000000  # m 
    #                            Set to this value to allow GS -> GEO connections (TODO: review max Ground-to-Satellite link range) 


# BASIC CONFIGURATION FUNCTIONS


# synchronize the data rates in config_ns3.properties with the values defined in the communication class
def update_ns3_rate_properties(config_path, rate_properties):
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

# Orbital radius
def radius(sat, epoch_str, date_str):
    
    # 1. Place the observer at the Earth's center (lat=0, lon=0, elevation=0)
    observer = ephem.Observer()
    observer.epoch = str(epoch_str).replace('-', '/')
    observer.date = str(date_str).replace('-', '/')
    observer.lat = 0
    observer.lon = 0
    observer.elevation = 0

    # 2. Compute satellite position relative to observer
    sat.compute(observer)

    # 3. Extract Geocentric Right Ascension (RA), Declination (DEC), and Range (distance to Earth center)
    # ra = float(sat.g_ra)       # Right Ascension in radians
    # dec = float(sat.g_dec)     # Declination in radians
    # r = float(sat.range)*10**(-3)       # Distance from observer/center in km

    try:
        ra = float(sat.g_ra)
        dec = float(sat.g_dec)
        r = float(sat.range)
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

