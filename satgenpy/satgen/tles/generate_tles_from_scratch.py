# The MIT License (MIT)
#
# Copyright (c) 2020 ETH Zurich
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# GENERATE SATELLITE CONSTELLATION TLEs

##############################################################################################################################################
##############################################################################################################################################
# WARNING: Epoch has been modified, as well as mean anomallies (so they don't necessarily need to be equally spaced )
# and satellites offset so it is possible to join LEO and GEO satellites in one tle.txt with coherent IDs
##############################################################################################################################################
##############################################################################################################################################

import math
from sgp4.exporter import export_tle 
from sgp4.api import Satrec, WGS72
from astropy.time import Time
from satgen.edrs.config import conf
from sgp4.api import jday
        # export_tle = converts satellite record into TLE text lines
        # Satrec = satellite record object
        # WGS72 = gravity model
        # jday = computes Julian day

def generate_tles_from_scratch_with_sgp(  # uses sgp4 library
        filename_out,
        constellation_name,
        num_orbits,
        num_sats_per_orbit,
        phase_diff,
        inclination_degree,
        eccentricity,
        arg_of_perigee_degree,
        mean_motion_rev_per_day
):

    with open(filename_out, "w+") as f_out:

        # First line:
        #
        # <number of orbits> <number of satellites per orbit>
        #
        f_out.write("%d %d\n" % (num_orbits, num_sats_per_orbit))

        # Each of the subsequent (number of orbits * number of satellites per orbit) blocks
        # define a satellite as follows:
        #
        # <constellation_name> <global satellite id>
        # <TLE line 1>
        # <TLE line 2>
        satellite_counter = 0
        for orbit in range(0, num_orbits): # each iteration creates an orbital plane

            # Orbit-dependent
            raan_degree = orbit * 360.0 / num_orbits # orbital planes have a RAAN equally spaced to cover the 360º
            orbit_wise_shift = 0
            if orbit % 2 == 1:
                if phase_diff:
                    orbit_wise_shift = 360.0 / (num_sats_per_orbit * 2.0)

            # For each satellite in the orbit
            for n_sat in range(0, num_sats_per_orbit):
                mean_anomaly_degree = orbit_wise_shift + (n_sat * 360 / num_sats_per_orbit) # distributes no.sat/plane equally space around the orbit, and adds an orbital_shift in odd_numbered planes so planos contiguos tienen los satélites intercalados entre planos orbitales

                # Epoch is set in configuration variables
                epoch = conf.epoch
                # Convert Astropy Time object 'epoch' to Python datetime components
                dt = epoch.datetime
                jd, fr = jday(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)

                # creates a SGP4 with all the data needed and then exports it into TLE, which is divided into line 1 and line 2 and then line 1 is debugged
                # Use SGP-4 to generate TLE
                sat_sgp4 = Satrec()

                # Based on: https://pypi.org/project/sgp4/
                sat_sgp4.sgp4init(
                    WGS72,                  # Gravity model [1]
                    'i',                    # Operating mode (a = old AFPSC mode, i = improved mode)
                    satellite_counter + 1,  # satnum:  satellite number
                    (jd + fr) - 2433281.5,  # epoch:   days since 1949 December 31 00:00 UT [2]
                    0.0,                    # bstar:   drag coefficient (kg/m2er)
                    0.0,                    # ndot:    ballistic coefficient (revs/day)
                    0.0,                    # nndot:   second derivative of mean motion (revs/day^3)
                    eccentricity,           # ecco:    eccentricity
                    math.radians(arg_of_perigee_degree),              # argpo:   argument or perigee (radians)
                    math.radians(inclination_degree),                 # inclo:    inclination(radians)
                    math.radians(mean_anomaly_degree),                # mo:       mean anomaly (radians)
                    mean_motion_rev_per_day * 60 / 13750.9870831397,  # no_kazai: mean motion (radians/minute) [3]
                    math.radians(raan_degree)                         # nodeo:    right ascension of
                                                                      #           ascending node (radians)
                )

                # Side notes:
                # [1] WGS72 is also used in the NS-3 model
                # [2] Due to a bug in sgp4init, the TLE below irrespective of the value here gives zeros.
                # [3] Conversion factor from:
                #     https://www.translatorscafe.com/unit-converter/en-US/ (continue on next line)
                #     velocity-angular/1-9/radian/second-revolution/day/
                #

                # Export TLE from the SGP-4 object
                line1, line2 = export_tle(sat_sgp4)

                """
                # Line 1 has some problems: there are unknown characters entered for the international
                # designator, and the Julian date is not respected
                # As such, we set our own bogus international designator 00000ABC
                # and we set our own epoch date as 1 January, 2000
                # Why it's 00001.00000000: https://www.celestrak.com/columns/v04n03/#FAQ04
                tle_line1 = line1[:7] + "U 00000ABC 00001.00000000 " + line1[33:]
                tle_line1 = tle_line1[:68] + str(calculate_tle_line_checksum(tle_line1[:68]))
                tle_line2 = line2
                """
                ##########
                # Modified epoch
                dt = conf.epoch.datetime
                yy = dt.year % 100
                day_of_year = dt.timetuple().tm_yday
                fractional_day = (dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6) / 86400.0

                # 1. Format epoch string (14 chars)
                epoch_str = "%02d%03d.%08d" % (yy, day_of_year, int(round(fractional_day * 1e8)))

                # 2. Build EXACT 68-character Line 1 without checksum
                satellite_id = satellite_counter + 1
                tle_line1_68 = "1 %05dU 00000ABC %s  .00000000  00000-0  00000+0 0    0" % (
                    satellite_id,
                    epoch_str
                )

                # 3. Calculate Modulo-10 checksum digit and append it
                checksum_digit = calculate_tle_line_checksum(tle_line1_68)
                tle_line1 = tle_line1_68 + str(checksum_digit)
                tle_line2 = line2

                ##########
                # Check that the checksum is correct ( each TLE line must be 68 elements)
                if len(tle_line1) != 69 or calculate_tle_line_checksum(tle_line1[:68]) != int(tle_line1[68]):
                    raise ValueError("TLE line 1 checksum failed")
                if len(tle_line2) != 69 or calculate_tle_line_checksum(tle_line2[:68]) != int(tle_line2[68]):
                    raise ValueError("TLE line 2 checksum failed")

                # Write TLE to file (saves Tle file for each satellite)
                f_out.write(constellation_name + " " + str(orbit * num_sats_per_orbit + n_sat) + "\n")
                f_out.write(tle_line1 + "\n")
                f_out.write(tle_line2 + "\n")

                # One more satellite there
                satellite_counter += 1


def generate_tles_from_scratch_manual( # same as before but instead using sgp4 for exporting the TLE, it builds the TLE lines by string formatting
        filename_out,
        constellation_name,
        num_orbits,
        num_sats_per_orbit,
        phase_diff,
        inclination_degree,
        eccentricity,
        arg_of_perigee_degree,
        mean_motion_rev_per_day,
        satellite_id_offset,
        mean_anomaly_degree_list=None,  # Added for cases in which satellites are not equially distributed in the orbit(EDRS)
):

    with open(filename_out, "w+") as f_out:

        # First line:
        #
        # <number of orbits> <number of satellites per orbit>
        #
        f_out.write("%d %d\n" % (num_orbits, num_sats_per_orbit))

        # Each of the subsequent (number of orbits * number of satellites per orbit) blocks
        # define a satellite as follows:
        #
        # <constellation_name> <global satellite id>
        # <TLE line 1>
        # <TLE line 2>
        satellite_counter = 0
        for orbit in range(0, num_orbits):
            # Orbit-dependent
            raan_degree = orbit * 360.0 / num_orbits
            orbit_wise_shift = 0
            if orbit % 2 == 1:
                if phase_diff:
                    orbit_wise_shift = 360.0 / (num_sats_per_orbit * 2.0)

          # For each satellite in the orbit
            for n_sat in range(0, num_sats_per_orbit):
                satellite_id = satellite_counter + 1 + satellite_id_offset 
                # OPTION 1: Equally spaced
                mean_anomaly_degree = orbit_wise_shift + (n_sat * 360 / num_sats_per_orbit)

                ##############################################################################################################################
                ##############################################################################################################################

                # OPTION 2: if there is a list, use it
                if mean_anomaly_degree_list is not None:
                    if len(mean_anomaly_degree_list) != num_orbits * num_sats_per_orbit:
                        raise ValueError(
                            "mean_anomaly_degree_list must have exactly "
                            f"{num_orbits * num_sats_per_orbit} elements"
                        )
                    mean_anomaly_degree = mean_anomaly_degree_list[satellite_counter]

                # Epoch is set based on configuration parameters
                dt = conf.epoch.datetime

                # Two digits year (e.g. 2027 to 27)
                yy = dt.year % 100

                # Day of the year + fraction of day
                day_of_year = dt.timetuple().tm_yday
                fractional_day = (dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6) / 86400.0

                # String of 14 characters (YYDDD.DDDDDDDD)
                epoch_str = "%02d%03d.%08d" % (
                    yy, 
                    day_of_year, 
                    int(round(fractional_day * 1e8))
                )

                # Build TLE1 (exact 68 characters before checksum)
                tle_line1 = "1 %05dU 00000ABC %s  .00000000  00000-0  00000+0 0    0" % (
                    satellite_id,
                    epoch_str
                )
                ##############################################################################################################################
                ##############################################################################################################################

                tle_line2 = "2 %05d %s %s %s %s %s %s    0" % (
                    satellite_id,
                    ("%3.4f" % inclination_degree).rjust(8),
                    ("%3.4f" % raan_degree).rjust(8),
                    ("%0.7f" % eccentricity)[2:],
                    ("%3.4f" % arg_of_perigee_degree).rjust(8),
                    ("%3.4f" % mean_anomaly_degree).rjust(8),
                    ("%2.8f" % mean_motion_rev_per_day).rjust(11),
                )

                # Append checksums
                tle_line1 = tle_line1 + str(calculate_tle_line_checksum(tle_line1))
                tle_line2 = tle_line2 + str(calculate_tle_line_checksum(tle_line2))

                # Write TLE to file
                f_out.write(constellation_name + " " + str(orbit * num_sats_per_orbit + n_sat + satellite_id_offset) + "\n")
                f_out.write(tle_line1 + "\n")
                f_out.write(tle_line2 + "\n")

                # One more satellite there
                satellite_counter += 1

def calculate_tle_line_checksum(tle_line_without_checksum):
    if len(tle_line_without_checksum) != 68:
       raise ValueError("Must have exactly 68 characters")
    s = 0
    for i in range(len(tle_line_without_checksum)):
        if tle_line_without_checksum[i].isnumeric():
            s += int(tle_line_without_checksum[i])
        if tle_line_without_checksum[i] == "-":
            s += 1
    return s % 10
