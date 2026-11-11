import shutil
import subprocess
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
satgenpy_dir = os.path.abspath(os.path.join(current_dir, "..", ".."))

if satgenpy_dir not in sys.path:
    sys.path.insert(0, satgenpy_dir)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from sgp4.api import Satrec, jday
from datetime import datetime # sgp4 works with julian days: jd, fr (julian integer, fraction part)

from src.tles.create_tles import *
from satgen.tles import read_tles
from satgen.ground_stations.read_ground_stations import read_ground_stations_basic
from satgen.ground_stations.extend_ground_stations import extend_ground_stations
from config import *
from src.txt_reader import *
from src.verifications import *
from src.simple_isls import *
from satgen.isls import read_isls
from satgen.dynamic_state import generate_dynamic_state

# 1. TLEs
# EO Satellite
# Choose Copernicus constellation and satellite
eo_sat = geometric.leo_sat
print('Generating TLEs for LEO Constellation...')
print('Choose Copernicus constellation to simulate, between the following options: Sentinel-1, Sentinel-2, Sentinel-3')
sentinel = input()
eo_tles(sentinel)

print(f'Choose {sentinel} satellite: 1 or 2')
sat_eo = int(input())
leo_id = sat_id(sentinel, sat_eo)
print(f"Copernicus Satellite ID:{leo_id}")

# EDRS
edrs_tles()

tles_dir = os.path.join(os.path.dirname(__file__),"src", "tles")
if sentinel == "Sentinel-1":
    file1 = os.path.join(tles_dir, "sentinel1_tle.txt")
elif sentinel == "Sentinel-2":
    file1 = os.path.join(tles_dir, "sentinel2_tle.txt")
elif sentinel == "Sentinel-3":
    file1 = os.path.join(tles_dir, "sentinel3_tle.txt")

print('Generating TLEs for GEO Constellation...')
file2 = os.path.join(tles_dir, "edcrs_tle.txt")
merge_tle_files(file1, file2)

print(f'Choose EDRS satellite: 1 for EDRS-A or 2 for EDRS-C')
sat_edrs = int(input()) + eo_sat
geo_id = sat_id("edrs", sat_edrs)
print(f"EDRS Satellite ID:{geo_id}")


# 2. Ground Stations
# TODO : Review Hawaii ground station use and data
base_dir = os.path.dirname(os.path.dirname(__file__)) # os.path lo que hace es subir niveles (file -> folder -> folder) 
gs_file = os.path.join(base_dir, "ground_stations", "ground_stations.txt")
gs_list = read_ground_stations_basic(gs_file)

# Verification of TLEs (through sgp4 = orbit propagator w/perturbations)
    # Steps: TLE -> Satrec -> sgp4
    # Sat IDs in these TLEs are wrong because it is date is taken before merging

filename = os.path.join(os.path.dirname(__file__),"src", "tles", "tle.txt")
leo_tle1, leo_tle2 = read_sat_tle(filename, sentinel, sat_eo)
geo_tle1, geo_tle2 = read_sat_tle(filename, "EDCRS", sat_edrs)

leo = Satrec.twoline2rv(leo_tle1, leo_tle2) # satrec transforms tle strings into a propagable satellite
geo = Satrec.twoline2rv(geo_tle1, geo_tle2)

dt = conf.epoch.utc.datetime   # check the satellites at the simulated epoch
jd, fr = jday(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)

e1, r_leo, v_leo = leo.sgp4(jd, fr)
e2, r_geo, v_geo = geo.sgp4(jd, fr)

if e1 != 0:
    raise RuntimeError(f"LEO propagation error: {e1}")
if e2 != 0:
    raise RuntimeError(f"GEO propagation error: {e2}")

d = distance(r_leo, r_geo)

print(f"LEO: ({r_leo[0]:.2f}, {r_leo[1]:.2f}, {r_leo[2]:.2f})")
print(f"GEO: ({r_geo[0]:.2f}, {r_geo[1]:.2f}, {r_geo[2]:.2f})")
print(f"Distance [km]: {d:.2f}")
        # LEO altitude: 160–2000 km
        # GEO altitude: 35786 km

# 3. Inter-satellite links lists
leo_id_tot, _ = analyze_tle_file(file1)
geo_id_tot, _ = analyze_tle_file(file2)

list_isls_leo_geo = generate_leo_geo_isls(leo_id_tot, geo_id_tot)
list_isls_leo_leo = generate_leo_leo_isls(leo_id_tot)

# 4. Propagation
# 4.1 Define paths and output directory
tle_file_path = os.path.join(os.path.dirname(__file__), "src", "tles", "tle.txt")
output_dynamic_state_dir = os.path.join(os.path.dirname(__file__), "src", "output")

# 4.2. Configure time parameters (single source of truth: config.simulation, also used for ns-3)
epoch = conf.epoch
simulation_end_time_ns = int(simulation.duration * 1_000_000_000)
offset_ns = 0   # ns-3 starts at the TLE epoch, so Hypatia must start there too
print(f"Simulation: {simulation.duration} s, step {simulation.time_step / 1e9:g} s "
      f"-> {simulation_end_time_ns // simulation.time_step} fstate files")

# 4.3. Load input topology elements
tles_info = read_tles(tle_file_path)
satellites = tles_info["satellites"]
ground_stations = gs_list

gsl_base_rate_mbps = communication.gsl_base_rate_mbps
leo_gsl_fraction = communication.leo_gs_rate_mbps / gsl_base_rate_mbps
geo_gsl_fraction = communication.geo_gs_downlink_rate_mbps / gsl_base_rate_mbps
ground_station_gsl_fraction = (
    communication.ground_station_gsl_rate_mbps / gsl_base_rate_mbps
)
list_gsl_interfaces_info = (
    [{"number_of_interfaces": 1, "aggregate_max_bandwidth": leo_gsl_fraction}]
    * geometric.leo_sat
    + [{"number_of_interfaces": 1, "aggregate_max_bandwidth": geo_gsl_fraction}]
    * (len(satellites) - geometric.leo_sat)
    + [{"number_of_interfaces": 1, "aggregate_max_bandwidth": ground_station_gsl_fraction}]
    * len(ground_stations)
)

# 4.4. Range limits
max_gsl_length_m = visibility.max_gsl_length
max_isl_length_m = visibility.max_range * 1000   # km -> m (feasibility itself is decided in visibility.py)

# 4.5. Remove dynamic-state files of previous runs (a shorter run would otherwise leave old files behind)
os.makedirs(output_dynamic_state_dir, exist_ok=True)
for f in os.listdir(output_dynamic_state_dir):
    if f.startswith(("fstate_", "gsl_if_bandwidth_")) and f.endswith(".txt"):
        os.remove(os.path.join(output_dynamic_state_dir, f))

# 5. GENERATE DYNAMIC STATE FILES. Call generate_dynamic_state. 
print("Generating dynamic state files...")
generate_dynamic_state (
    output_dynamic_state_dir= output_dynamic_state_dir,
    epoch= epoch,
    simulation_end_time_ns= simulation_end_time_ns,
    time_step_ns= simulation.time_step,
    offset_ns= offset_ns,
    satellites= satellites,
    ground_stations= ground_stations,
    list_isls_leo_geo = list_isls_leo_geo,
    list_isls_leo_leo = list_isls_leo_leo ,
    list_gsl_interfaces_info= list_gsl_interfaces_info,
    max_gsl_length_m= max_gsl_length_m,
    max_isl_length_m= max_isl_length_m,
    dynamic_state_algorithm= "algorithm_free_one_only_over_isls",
    enable_verbose_logs= True
)

os.makedirs(output_dynamic_state_dir, exist_ok=True)

# Copy topology files into the ns-3 run directory.
shutil.copy2(tle_file_path, os.path.join(output_dynamic_state_dir, "tles.txt"))
extend_ground_stations(
    gs_file,
    os.path.join(output_dynamic_state_dir, "ground_stations.txt")
)

# Combine both ISL lists into ns-3's "a b" format.
with open(os.path.join(output_dynamic_state_dir, "isls.txt"), "w") as f_out:
    for a, b in list_isls_leo_geo + list_isls_leo_leo:
        f_out.write(f"{a} {b}\n")

# Write one interface row per node: satellites first, then ground stations.
with open(os.path.join(output_dynamic_state_dir, "gsl_interfaces_info.txt"), "w") as f_out:
    for node_id, info in enumerate(list_gsl_interfaces_info):
        f_out.write(
            f"{node_id},{info['number_of_interfaces']},"
            f"{info['aggregate_max_bandwidth']}\n"
        )

# 6. RUN NS-3 SIMULATION VIA SUBPROCESS
print("Feeding state files into ns-3 network simulator...")

# Define the absolute path to your ns-3 installation directory
current_dir = os.path.dirname(os.path.abspath(__file__))
hypatia_root = os.path.abspath(os.path.join(current_dir, "..", "..", ".."))
ns3_dir = os.path.join(hypatia_root, "ns3-sat-sim", "simulator")
waf_executable = os.path.join(ns3_dir, "waf")

# Absolute path to output dynamic state directory
dynamic_state_abs_dir = os.path.abspath(output_dynamic_state_dir)

# Path to template config_ns3.properties
# Typically found in satgenpy templates or another example scenario folder
template_config_path = os.path.join(current_dir, "src", "output", "config_ns3.properties.template")
destination_config_path = os.path.join(dynamic_state_abs_dir, "config_ns3.properties")

# Always rebuild config_ns3.properties from the template, so ns-3 runs over exactly the same
# time span and time step as the Hypatia propagation above
shutil.copyfile(template_config_path, destination_config_path)

ns3_properties = {
    "simulation_end_time_ns": simulation_end_time_ns,
    "dynamic_state_update_interval_ns": simulation.time_step,
    "isl_data_rate_megabit_per_s": communication.isl_rate_mbps,
    "gsl_data_rate_megabit_per_s": communication.gsl_base_rate_mbps,
}
if ns3.enable_pingmesh:
    n_sat = len(satellites)
    pairs = ", ".join(f"{n_sat + a}->{n_sat + b}" for a, b in ns3.pingmesh_pairs_gid)
    ns3_properties.update({
        "enable_pingmesh_scheduler": "true",
        "pingmesh_interval_ns": ns3.pingmesh_interval_ns,
        "pingmesh_endpoint_pairs": f"set({pairs})",
    })
update_ns3_properties(destination_config_path, ns3_properties)
print(f"ns-3 config: {simulation.duration} s, update every {simulation.time_step / 1e9:g} s")
    
waf_cmd = [
    waf_executable,
    "--run",
    f"main_satnet --run_dir={dynamic_state_abs_dir}"
]

print(f"Launching ns-3 from directory: {ns3_dir}")

try:
    process = subprocess.Popen(
        waf_cmd,
        cwd=ns3_dir,  # Ensures ns-3 context executes inside the ns-3 folder
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    for line in process.stdout:
        print(line, end="")

    process.wait()
    if process.returncode == 0:
        print("\nSUCCESS: ns-3 simulation finished successfully")
    else:
        print(f"\nERROR: ns-3 failed with return code {process.returncode}")

except Exception as e:
    print(f"\nFailed to launch ns-3: {e}")