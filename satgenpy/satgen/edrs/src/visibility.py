import numpy as np
from satgen.edrs.config import geometric, visibility
# Analysis based on following paper:
# LEO–GEO visibility. Satellite Communications and Networks. Wiley. https://doi.org/10.1002/sat.70066
#
# UNITS: r_leo, r_geo, r_a, r_b are geocentric position vectors in METRES (as returned by config.radius()).
#        config values are given in km and are converted here.

_RE_M = geometric.Re * 1e3          # km -> m
_H_GRAZE_M = geometric.ATM * 1e3    # km -> m
_MAX_RANGE_M = visibility.max_range * 1e3


def angle(a, b):
    """Angle between two vectors in degrees."""
    cos = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return np.rad2deg(np.arccos(np.clip(cos, -1.0, 1.0)))


def los_min_distance(r1, r2):
    """Minimum distance from Earth's centre to the segment r1 -> r2 (m)."""
    d = r2 - r1                                   # r1 -> r2
    t = -np.dot(r1, d) / np.dot(d, d)             # closest point P(t) = r1 + t d
    t = min(max(t, 0.0), 1.0)                     # clamp to the segment
    return np.linalg.norm(r1 + t * d)


def contact_leo_geo(r_leo, r_geo):
    """True if the LEO-GEO link is geometrically feasible."""
    d = r_geo - r_leo                             # LEO -> GEO

    # 1. Earth + atmosphere blockage
    if los_min_distance(r_leo, r_geo) <= _RE_M + _H_GRAZE_M:
        return False

    # 2. GEO terminal: off-NADIR angle of the LEO seen from GEO
    if angle(-d, -r_geo) >= visibility.theta_geo_max:
        return False

    # 3. LEO terminal: off-ZENITH angle of the GEO seen from LEO
    if angle(d, r_leo) >= visibility.theta_leo_max:
        return False

    # 4. Range
    if np.linalg.norm(d) > _MAX_RANGE_M:
        return False

    return True


def contact_leo_leo(r_a, r_b):
    """True if the LEO-LEO link is geometrically feasible."""
    if los_min_distance(r_a, r_b) <= _RE_M + _H_GRAZE_M:
        return False
    if np.linalg.norm(r_b - r_a) > _MAX_RANGE_M:
        return False
    return True


def contact_gs_sat(ground_station, satellite, epoch_str, date_str):
    """
    Ground station <-> satellite link feasibility.
    Returns (feasible, slant_range_m). Feasible if the satellite is above the ground station's
    elevation mask (visibility.min_elevation_gs) and within visibility.max_gsl_length.
    """
    import ephem
    observer = ephem.Observer()
    observer.epoch = epoch_str
    observer.date = date_str
    observer.lat = str(ground_station["latitude_degrees_str"])
    observer.lon = str(ground_station["longitude_degrees_str"])
    observer.elevation = ground_station["elevation_m_float"]
    satellite.compute(observer)
    elevation_deg = np.degrees(float(satellite.alt))
    slant_range_m = float(satellite.range)
    feasible = (elevation_deg >= visibility.min_elevation_gs) and (slant_range_m <= visibility.max_gsl_length)
    return feasible, slant_range_m
