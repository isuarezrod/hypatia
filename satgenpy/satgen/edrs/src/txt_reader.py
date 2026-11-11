import os

# Satellite identification
def sat_id(constellation, sat):
    base_dir = os.path.join(os.path.dirname(__file__), "tles") 
    # "Take folder where current file is (config.py) and search for tles folder, put it in there"
    if constellation == "Sentinel-1":
        tle = os.path.join(base_dir, "sentinel1_tle.txt")
    elif constellation == "Sentinel-2":
        tle = os.path.join(base_dir, "sentinel2_tle.txt")
    elif constellation == "Sentinel-3":
        tle = os.path.join(base_dir, "sentinel3_tle.txt")
    elif constellation == "edrs":
        tle = os.path.join(base_dir, "tle.txt")

    with open(tle, "r") as f:
        data = f.readlines()

    row = data[3 * sat - 1]
    satellite_id = row[2:7]

    return satellite_id

# Analyze TLE file
def analyze_tle_file(tle_path):
    leo_ids = []
    with open(tle_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # If first row begins like " 1 2 " we skip it
    start_idx = 0
    if len(lines) > 0 and len(lines[0].split()) == 2 and lines[0].split()[0].isdigit() and lines[0].split()[1].isdigit():
        start_idx = 1

    # Each satellite is made of three rows (name, TLE 1, TLE 2)
    for i in range(start_idx, len(lines), 3):
        if i + 1 < len(lines):
            tle_line_1 = lines[i + 1]
            # verify that TLE 1 begins with a "1 "
            if tle_line_1.startswith('1 '):
                # ID is from columns of strings 3 to 7 (or second word without 'U')
                # example: '1 00001U ...' -> '00001'
                id = tle_line_1.split()[1].rstrip('U').rstrip('C').rstrip('S')
                leo_ids.append(id)
    
    # total satellites
    num_sats = len(leo_ids)
    num_sats_int = int(num_sats)
    
    return leo_ids, num_sats_int

# Number of ground stations 

def get_num_ground_stations(gs_path):
    with open(gs_path, 'r') as f:
        
        lines = [line.strip() for line in f if line.strip()]

    # get last row, first column
    last_line = lines[-1]
    last_gs_id = int(last_line.split(',')[0])
    
    num_gs = last_gs_id + 1
    num_gs_str = str(num_gs)
    
    return num_gs_str