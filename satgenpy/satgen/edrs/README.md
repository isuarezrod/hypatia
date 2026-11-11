# EDRS – Copernicus relay simulation (Hypatia)

Simulates a Copernicus Sentinel satellite (LEO) talking to the European Data Relay System
(EDRS-A at 9°E, EDRS-C at 31°E, GEO) through laser inter-satellite links (ISLs), and feeds the
result into the ns-3 network simulator of Hypatia.

```
edrs/
├── main.py                 ← run this: builds the scenario, runs Hypatia + ns-3
├── config.py               ← all parameters (time, geometry, data rates, visibility limits)
├── src/
│   ├── tles/create_tles.py ← generates the TLEs of Sentinel and EDRS
│   ├── visibility.py       ← decides, at each instant, if a link is physically possible
│   ├── simple_isls.py      ← list of candidate links (every LEO–GEO pair, every LEO–LEO pair)
│   ├── txt_reader.py, verifications.py  ← small helpers
│   └── output/             ← everything main.py produces (input for ns-3)
├── validation/             ← Step 2: checks the geometry against SGP4 and GMAT
└── viz/                    ← 3D viewer of Earth + satellites + links
```

All commands below are typed in the **Ubuntu (WSL) terminal**, inside
`~/simulation/hypatia/satgenpy/satgen/edrs`, unless it says **Windows** (GMAT, browser).

---

## How `main.py` works

```bash
python main.py
```

It runs six stages, top to bottom.

### 1. TLEs (which satellites, which orbits)

`main.py` asks you three questions in the terminal:

| Prompt | Answer | Meaning |
|---|---|---|
| `Choose Copernicus constellation...` | `Sentinel-1`, `Sentinel-2` or `Sentinel-3` | which LEO constellation (2 satellites each) |
| `Choose Sentinel-X satellite: 1 or 2` | `1` or `2` | which of the two Sentinels you focus on |
| `Choose EDRS satellite: 1 for EDRS-A or 2 for EDRS-C` | `1` or `2` | which relay you focus on |

Behind the scenes:
- `create_tles.eo_tles()` writes `src/tles/sentinelX_tle.txt` (2 LEO satellites, sun-synchronous orbit).
- `create_tles.edrs_tles()` writes `src/tles/edcrs_tle.txt` (2 GEO satellites). The mean anomaly is
  computed from the Earth's rotation angle at the epoch so EDRS-A really sits over 9°E and EDRS-C over 31°E.
- `merge_tle_files()` joins both into `src/tles/tle.txt`: satellites 0–1 are the Sentinels, 2–3 the EDRS.
- All TLEs share the epoch of `config.conf.epoch` (2026-01-01 00:00 UTC).

### 2. Ground stations

Reads `satgen/ground_stations/ground_stations.txt` (Redu, Weilheim, Matera, Hawaii).

### 3. Quick SGP4 sanity check

Propagates the chosen Sentinel and EDRS satellite with SGP4 to one date and prints their positions and
distance. It is only a print-out to see that the TLEs are sensible (LEO ≈ 7,000 km from Earth's centre,
GEO ≈ 42,164 km).

### 4. Candidate inter-satellite links

`simple_isls.py` lists every link that *could* exist:
- LEO–GEO: every Sentinel with every EDRS → `(0,2) (0,3) (1,2) (1,3)`
- LEO–LEO: the two Sentinels with each other → `(0,1)`

These are only candidates. Whether a link is actually usable at a given moment is decided in stage 5.

### 5. Dynamic state: propagation + visibility (the core)

`generate_dynamic_state()` (in `satgen/dynamic_state/`) steps through time from 0 to
`simulation.duration` every `simulation.time_step`. At **each time step**:

1. **Propagation.** Hypatia (ephem library) propagates every satellite from its TLE.
   `config.radius()` returns each satellite's position vector from Earth's centre, in metres.
2. **Visibility.** For every candidate link, `src/visibility.py` decides True/False:

   | Check | LEO–GEO | LEO–LEO |
   |---|---|---|
   | Line of sight does not cross Earth + 100 km of atmosphere | ✔ | ✔ |
   | GEO terminal: LEO within `theta_geo_max` of the GEO's nadir | ✔ | – |
   | LEO terminal: GEO within `theta_leo_max` of the LEO's zenith | ✔ | – |
   | Distance ≤ `max_range` | ✔ | ✔ |

   Only links that pass every check are added to the network graph for that instant.
3. **Routing.** Hypatia computes the shortest paths through the graph and writes the forwarding
   tables `src/output/fstate_<t_ns>.txt` (only the entries that changed since the previous step).

The limits used by `visibility.py` are in `config.py` (`geometric` and `visibility` classes).

> ### ✅ Check the geometry: run the validation (Step 2)
>
> Whenever you change `visibility.py`, the TLEs or the limits in `config.py`, check that the
> contact windows and distances are right by comparing against two independent propagators:
> **SGP4** (analytical, same TLE) and **GMAT** (numerical). Three stages: terminal → GMAT → terminal.
>
> **1. Ubuntu (WSL) terminal**: compare Hypatia against SGP4 and create the GMAT script
> ```bash
> python -m validation.validate_geometry --duration 86400 --step 10
> python -m validation.make_gmat_script --out /mnt/c/Users/triiso/gmat_edrs
> ```
> The first line should end with `Overall: PASS`. The second writes the GMAT script on your C: drive.
>
> **2. Windows: GMAT**
> 1. Open GMAT.
> 2. **File → Open** → `C:\Users\triiso\gmat_edrs\edrs_validation.script`
> 3. Press **Run** (blue triangle, or F5).
> 4. Wait for `Mission run completed` in the message window. GMAT writes
>    `C:\Users\triiso\gmat_edrs\report_edrs.txt`. (Errors about Python/MATLAB plugins at start-up can be ignored.)
>
> **3. Back in the Ubuntu (WSL) terminal**: final check with GMAT included
> ```bash
> python -m validation.validate_geometry --duration 86400 --step 10 \
>     --gmat /mnt/c/Users/triiso/gmat_edrs/report_edrs.txt --tol-range-km-gmat 30
> ```
> Look at the last line: `Overall: PASS`. Results are in `validation/results/`
> (`report.md`, `geometry_validation.png`, CSVs). Open the folder with `explorer.exe validation/results`.
>
> Expected numbers: Hypatia vs SGP4 ≤ 0.5 km; GMAT vs SGP4 ≈ 1–13 km over 24 h (normal: different
> dynamics models); all contact windows identical within one time step.
> Details in [`validation/README.md`](validation/README.md).

### 6. ns-3 files and simulation

`main.py` copies everything ns-3 needs into `src/output/`:
`tles.txt`, `ground_stations.txt`, `isls.txt`, `gsl_interfaces_info.txt`, `config_ns3.properties`
(data rates are synchronised from `config.communication`), and then launches

```
waf --run "main_satnet --run_dir=<edrs>/src/output"
```

from `ns3-sat-sim/simulator`. ns-3 creates every ISL listed in `isls.txt` as a laser link and, at every
time step, routes packets according to `fstate_*.txt`. So **a link only carries traffic while
`visibility.py` says it is feasible**. Results land in `src/output/logs_ns3/`
(e.g. `isl_utilization.csv`).

### Things to keep in mind

- **Simulation time vs. contact windows.** `simulation.duration` is 200 s. The first LEO–GEO window of
  Sentinel 0 opens at ~00:34 UTC, of Sentinel 1 at 00:00 UTC (until ~00:23). Check
  `validation/results/report.md` and choose a duration / offset that falls inside a window,
  otherwise ns-3 will not see any ISL.
- **`main.py` currently crashes before stage 5** (when it builds the GSL interface list) at `communication.ground_station_gsl_rate_mbps`: that
  attribute does not exist in `config.py`. Add it to the `communication` class (the data rate you want
  for ground-station interfaces).
- The SGP4 check of stage 3 uses `datetime(2027, 1, 1)`, one year after the TLE epoch. Use
  `conf.epoch` instead so it checks the simulated day.

---

## Reading `validation/results/geometry_validation.png`

One panel per satellite pair (`0-2` = Sentinel-1 0 ↔ EDRS-A, `0-3` = Sentinel-1 0 ↔ EDRS-C, … ,
`0-1` = the two Sentinels). Horizontal axis: hours since the epoch.

| Element | What it is |
|---|---|
| **Blue line** | Distance between the two satellites computed with SGP4 (reference) |
| **Orange dashed line** | Same distance computed by Hypatia. If you cannot tell it apart from the blue line, they agree |
| **Tall light-blue bands** | Times when the link is feasible according to the reference (SGP4) |
| **Short orange bands** | Times when the link is feasible according to Hypatia (`visibility.py`) |

What to look for: **every blue band has an orange band at the same place, with the same width.**
A blue band without orange means Hypatia misses a contact; orange without blue means Hypatia invents one.

Why the distance waves up and down: the Sentinel orbits the Earth every ~98 min. When it is on the same
side of the Earth as the GEO, the distance drops to ~35,000 km and the link is usually feasible (band).
When it goes behind the Earth, the distance grows to ~49,000 km and the Earth blocks the line of sight
(no band). The pattern slowly changes over the day because the Earth rotates under the orbit.

The `0-1` panel (LEO–LEO) has no bands: the two Sentinels are 180° apart on the same orbit, so the
Earth is always between them. Its y-axis is zoomed in (≈ 14,143–14,148 km): the distance between them
barely changes.

---

## 3D viewer (Earth + satellites + links)

```bash
python -m viz.orbit_viewer                     # 24 h, one point every 30 s
explorer.exe viz/output/orbit_viewer.html      # opens it in your Windows browser
```

- The Earth rotates and the satellites move along their orbits (inertial frame).
  Tick **Earth-fixed view** to see it from the ground: the EDRS satellites stay over Europe/Africa.
- **Green lines** = links that are geometrically feasible at that instant (same criteria as `visibility.py`).
  The panel lists the active links and their length.
- Controls: drag to rotate, scroll to zoom, right-drag to pan; speed slider, time slider, pause.
- Options: `--duration 21600 --step 10` (6 h, finer), `--no-geometry-links`.
- Needs internet the first time it opens (it loads three.js and the Earth texture from a CDN).

**Adding the ns-3 links later.** The viewer reads links in this format (JSON list):

```json
[{"a": 1, "b": 2, "kind": "LEO-GEO", "source": "ns3", "windows": [[0, 1380], [4844, 7306]]}]
```

`a`, `b` = satellite indices as in `tle.txt`; `windows` = `[start_s, end_s]` intervals, in seconds from
the epoch, when ns-3 used that ISL. Write a small script that builds this list from
`src/output/logs_ns3/isl_utilization.csv` (or the `fstate_*.txt` files), save it as e.g. `ns3_links.json`, then:

```bash
python -m viz.orbit_viewer --links-json ns3_links.json
```

Links with `"source": "ns3"` are drawn in **red**, on top of the green geometric ones, so you can see
which feasible links ns-3 actually used.
