import os
from math import sqrt

def distance(r1, r2):
    d = sqrt((r1[0]-r2[0])**2 + 
             (r1[1]-r2[1])**2 + 
             (r1[2]-r2[2])**2)
    return d

def read_sat_tle(filename, sat_name, sat_index):
    with open(filename, "r") as f:
        lines = f.readlines()

    body = lines[1:]  # skip header

    for i in range(0, len(body), 3):
        name_line = body[i].strip()
        tle1 = body[i + 1].strip()
        tle2 = body[i + 2].strip()

        part1 = name_line.split()
        part2 = tle2.split()
        
        name = part1[0]
        idx = int(part2[1])

        if name == sat_name and idx == sat_index:
            return tle1, tle2

    raise ValueError(f"Satellite {sat_name} {sat_index} not found in {filename}")
