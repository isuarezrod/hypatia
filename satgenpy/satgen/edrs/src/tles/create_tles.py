import os
from satgen.tles.generate_tles_from_scratch import *
from satgen.edrs.config import *


def eo_tles(sentinel):
    # scenario_dir = "satgen/edrs/tles"
    scenario_dir = os.path.dirname(os.path.abspath(__file__)) # Saves tle.txt file in folder where create_tles.py is in, in this way it doesnt depend 
                                                              # on the directory, like it does on the line above
    if sentinel == "Sentinel-1":
        generate_tles_from_scratch_with_sgp(
            # Copernicus Sentinel-1
            filename_out = os.path.join(scenario_dir, "sentinel1_tle.txt"),
            constellation_name="Sentinel-1",
            num_orbits=1,
            num_sats_per_orbit=2,
            phase_diff=True,
            inclination_degree=98.18,
            eccentricity=0.0001,
            arg_of_perigee_degree=0.0,
            mean_motion_rev_per_day=14.60
        )
    elif sentinel == "Sentinel-2":
        generate_tles_from_scratch_with_sgp(
            # Copernicus Sentinel-2
            filename_out = os.path.join(scenario_dir, "sentinel2_tle.txt"),
            constellation_name="Sentinel-2",
            num_orbits=1,
            num_sats_per_orbit=2,
            phase_diff=True,
            inclination_degree=98.62,
            eccentricity=0.0001,
            arg_of_perigee_degree=0.0,
            mean_motion_rev_per_day=14.31
        )

    elif sentinel == "Sentinel-3":
        # correction: they are 140º apart
        generate_tles_from_scratch_manual(
            # Copernicus Sentinel-3
            filename_out = os.path.join(scenario_dir, "sentinel3_tle.txt"),
            constellation_name="Sentinel-3",
            num_orbits=1,
            num_sats_per_orbit=2,
            phase_diff=True,
            inclination_degree=98.65,
            eccentricity=0.0001,
            arg_of_perigee_degree=0.0,
            mean_motion_rev_per_day=14.26,
            satellite_id_offset=0,
            mean_anomaly_degree_list=[0, 140]
        )

    else: 
        print('Inserted constellation does not match any of Copernicus Sentinel LEO constellations')

    print("TLE file generated")


def edrs_tles():
    # one orbit, two satellites: EDRS-A and EDRS-C
    scenario_dir = os.path.dirname(os.path.abspath(__file__))

    generate_tles_from_scratch_manual(
        filename_out=os.path.join(scenario_dir, "edcrs_tle.txt"),
        constellation_name="EDCRS",
        num_orbits=1,
        num_sats_per_orbit= 2,
        phase_diff=True,
        inclination_degree= 0.01, # avoid singularity w/ sgp4
        eccentricity=0.0,
        arg_of_perigee_degree=0.0,
        mean_motion_rev_per_day=1.0027379,  # one rotation every sidereal day
        satellite_id_offset= geometric.leo_sat, # Fixed, all Sentinel constellations have 2 satellites
        mean_anomaly_degree_list=[9, 31]
    )

def merge_tle_files(file1, file2):
    scenario_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = os.path.join(scenario_dir, "tle.txt")
    with open(file1, "r") as f1:
        lines1 = f1.readlines()

    with open(file2, "r") as f2:
        lines2 = f2.readlines()

    with open(output_file, "w") as fout:
        # 'lines' are a list, we need to define variable 'header' to split the list1[1] = '1 2'
        header1 = lines1[0].split()
        header2 = lines2[0].split()

        total_sats_1 = int(header1[0]) * int(header1[1])
        total_sats_2 = int(header2[0]) * int(header2[1])
        num_orbits = int(header1[0]) + int(header2[0])
        total_sats = total_sats_1 
        # total_sats_leo constellation = total_sats_geo constellation

        fout.write("%d %d\n" % (num_orbits, total_sats))

        # skip first line of both input files
        fout.writelines(lines1[1:])
        fout.writelines(lines2[1:])
