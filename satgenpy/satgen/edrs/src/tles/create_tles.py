import os
from satgen.tles.generate_tles_from_scratch import *
from satgen.edrs.config import *

# All TLE files are written next to this script (src/tles/)
scenario_dir = os.path.dirname(os.path.abspath(__file__))


def geo_mean_anomaly_for_longitude(lon_east_deg, raan_deg=0.0, argp_deg=0.0):
    """
    For a (near-)equatorial circular GEO the sub-satellite longitude at the epoch is
        lon = RAAN + argp + M - GMST(epoch)
    so the mean anomaly that puts the satellite over lon_east_deg is
        M = lon + GMST(epoch) - RAAN - argp
    """
    gmst_deg = conf.epoch.sidereal_time("mean", "greenwich").deg
    return (lon_east_deg + gmst_deg - raan_deg - argp_deg) % 360.0


def create_all_tles():
    """
    Writes one TLE file per group (src/tles/<name>_tle.txt) for every entry in
    config.constellations, then merges them into src/tles/tle.txt.
    Returns the path of tle.txt.
    """
    files = []
    offset = 0   # global satellite index of the first satellite of each group

    # LEO groups (Sentinel-1, -2, -3, ...)
    for name, inclination, mean_motion, mean_anomalies in constellations.leo:
        filename = os.path.join(scenario_dir, name.lower().replace("-", "") + "_tle.txt")
        generate_tles_from_scratch_manual(
            filename_out=filename,
            constellation_name=name,
            num_orbits=1,
            num_sats_per_orbit=len(mean_anomalies),
            phase_diff=False,
            inclination_degree=inclination,
            eccentricity=0.0001,
            arg_of_perigee_degree=0.0,
            mean_motion_rev_per_day=mean_motion,
            satellite_id_offset=offset,
            mean_anomaly_degree_list=mean_anomalies,
        )
        files.append(filename)
        offset += len(mean_anomalies)

    # GEO relays (EDRS): one orbit, satellites placed over their longitude at the epoch
    filename = os.path.join(scenario_dir, "edrs_tle.txt")
    generate_tles_from_scratch_manual(
        filename_out=filename,
        constellation_name="EDRS",
        num_orbits=1,
        num_sats_per_orbit=len(constellations.geo),
        phase_diff=False,
        inclination_degree=0.01,            # avoid the i = 0 singularity in SGP4
        eccentricity=0.0,
        arg_of_perigee_degree=0.0,
        mean_motion_rev_per_day=1.0027379,  # one revolution per sidereal day
        satellite_id_offset=offset,
        # NOT the longitude itself: the mean anomaly is inertial, the longitude is Earth-fixed
        mean_anomaly_degree_list=[geo_mean_anomaly_for_longitude(lon) for _, lon in constellations.geo],
    )
    files.append(filename)

    return merge_tle_files(files)


def merge_tle_files(files):
    """Merges Hypatia TLE files into src/tles/tle.txt."""
    output_file = os.path.join(scenario_dir, "tle.txt")
    bodies = []
    for filename in files:
        with open(filename, "r") as f:
            bodies += f.readlines()[1:]          # skip each file's header line
    total_sats = len(bodies) // 3
    with open(output_file, "w") as fout:
        # Header "<orbits> <satellites per orbit>": Hypatia and ns-3 only need orbits x sats = total
        fout.write("%d 1\n" % total_sats)
        fout.writelines(bodies)
    print(f"TLE file with {total_sats} satellites written to {output_file}")
    return output_file
