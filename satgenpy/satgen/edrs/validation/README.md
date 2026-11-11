# Step 2 – Geometry validation

Checks that the contact windows and ranges produced by the Hypatia path
(`ephem` → `config.radius()` → `src/visibility.py` / `distance_tools`) match an
independent propagation within tolerance.

| source | propagation | geometry | purpose |
|---|---|---|---|
| **hypatia** | ephem (libastro SGP4/SDP4) | `src/visibility.py` | system under test |
| **ref** | python-sgp4 (Vallado) | `validation/geometry.py` (no satgen imports) | independent reference |
| **gmat** | GMAT RK89, JGM-3 8×8 from the SGP4 state at epoch | `validation/geometry.py` | independent dynamics |
| **baseline** | whatever the paper used (ISRO) | – | external windows |

All commands run on the **host**, inside `satgenpy/satgen/edrs` (same Python env you use for `main.py`).

```bash
# 1. Hypatia vs SGP4 (24 h, 10 s). Exit code 0 = PASS
python -m validation.validate_geometry --duration 86400 --step 10

# 2. GMAT: generate the script, run it, feed the report back
python -m validation.make_gmat_script --duration 86400 --step 10
GmatConsole -r validation/gmat/edrs_validation.script        # or open it in the GMAT GUI and Run
python -m validation.validate_geometry --duration 86400 --step 10 --gmat validation/gmat/report_edrs.txt

# 3. External baseline (ISRO windows from the paper)
python -m validation.validate_geometry --baseline validation/isro_baseline.csv
```

`isro_baseline.csv` format (seconds since the epoch in `config.conf.epoch`):

```
pair,start_s,end_s
0-2,1830.0,4380.0
```

The baseline only makes sense if the scenario (orbits, epoch, FOV limits, grazing altitude)
is set to the paper's case. Otherwise compare statistics (number of windows/day, mean duration, range span)
rather than individual windows.

## Default tolerances and why

| test | range | window edge | origin of the error |
|---|---|---|---|
| hypatia vs sgp4 | 1 km | = step | two SGP4 implementations (~0.3–0.5 km) + sampling (edges quantised to the step) |
| gmat vs sgp4 | 10 km | 30 s | analytical mean theory vs numerical; grows with time (~3–4 km/day in LEO with J2 only) |

Outputs in `validation/results/`: `report.md`, `windows.csv`, `window_comparison.csv`, `ranges.csv`, `geometry_validation.png`.
