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


def generate_plus_grid_isls(output_filename_isls, n_orbits, n_sats_per_orbit, isl_shift, idx_offset=0): # creates an ISL topology instead of reading one
    """
    Generate plus grid ISL file.

    :param output_filename_isls     Output filename
    :param n_orbits:                Number of orbits
    :param n_sats_per_orbit:        Number of satellites per orbit
    :param isl_shift:               ISL shift between orbits (e.g., if satellite id in orbit is X,
                                    does it also connect to the satellite at X in the adjacent orbit)
                        isl_shift = 0 = connect straight across to the same position in next orbit
                        isl_shift = 1 = connect to the next-position satellite in next orbit
    :param idx_offset:              Index offset (e.g., if you have multiple shells) (e.g: shell 1 IDs from 0 - 999, shell 2 IDs from 1000 - 1999)
    """

    if n_orbits < 3 or n_sats_per_orbit < 3:
        raise ValueError("Number of x and y must each be at least 3")

    list_isls = []
    for i in range(n_orbits):
        for j in range(n_sats_per_orbit):
            sat = i * n_sats_per_orbit + j # creates satellite index

            # Link to the next in the orbit
            sat_same_orbit = i * n_sats_per_orbit + ((j + 1) % n_sats_per_orbit) # stay in i-th orbit indexs, move to next slot: j+1, and wrap around using modulo: 
            
            # Link to the next orbit                                                                            # so for j=3, if no.sat/orb = 4, sat 3 links with sat 0, like a ring
            sat_adjacent_orbit = ((i + 1) % n_orbits) * n_sats_per_orbit + ((j + isl_shift) % n_sats_per_orbit) # go to orbit i+1, and wrap around w/ modulo so it comes around, and last orbit connects back to orbit 0
                                                                                                                # in that orbit connect to position j + isl_shift
            # Same orbit
            list_isls.append((idx_offset + min(sat, sat_same_orbit), idx_offset + max(sat, sat_same_orbit)))

            # Adjacent orbit
            list_isls.append((idx_offset + min(sat, sat_adjacent_orbit), idx_offset + max(sat, sat_adjacent_orbit)))

    with open(output_filename_isls, 'w+') as f:
        for (a, b) in list_isls:
            f.write(str(a) + " " + str(b) + "\n")

    return list_isls

