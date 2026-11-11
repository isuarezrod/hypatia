from satgen.distance_tools.distance_tools import distance_m_between_satellites
import numpy as np
from satgen.edrs.config import geometric, visibility
# Analysis based on followig paper:
# LEO–GEO visibility. Satellite Communications and Networks. Wiley. https://doi.org/10.1002/sat.70066

def angle(a, b):
    cos = np.dot(a, b)/ (np.linalg.norm(a)*np.linalg.norm(b))
    return np.rad2deg(np.arccos(cos))

def contact_leo_geo(r_leo, r_geo):

    # Returns True if there is contact LEO - GEO
    theta_geo_max = visibility.theta_geo_max
    theta_leo_max = visibility.theta_leo_max
    max_range = visibility.max_range

    # 1. Earth block
    d = r_leo - r_geo # distance from LEO to GEO 
    t = -np.dot(r_leo, d)/np.dot(d,d)
    if t < 0:
        d_min = r_leo
    elif t > 1:
        d_min = r_geo
    else:
        d_min = np.cross(r_leo, r_geo)/np.dot(d,d)

    d_min_val = np.linalg.norm(d_min)
    if d_min_val <= geometric.Re + geometric.ATM:
        return False

    # 2. GEO pointing to LEO
    theta_geo = angle(d, r_geo) 
    if theta_geo >= theta_geo_max:
        return False

    # 3. LEO pointing GEO
    theta_leo = angle(-d, -r_leo)
    if theta_leo >= theta_leo_max:
        return False

    # 4. Range
    if np.linalg.norm(d) > max_range:
        return False


    return True

def contact_leo_leo(r_a, r_b):
    max_range = visibility.max_range

     # 1. Earth blocking
    d = r_a - r_b # distance from LEO to GEO 
    t = -np.dot(r_a, d)/np.dot(d,d)
    if t < 0:
        d_min = np.linalg.norm(r_a)
    elif t > 1:
        d_min = np.linalg.norm(r_b)
    else:
        d_min = np.linalg.norm(np.cross(r_a, r_b)/np.dot(d,d))

    if d_min <= geometric.Re + geometric.ATM:
        return False

     # 2. Range
    if np.linalg.norm(d) > max_range:
        return False

    return True