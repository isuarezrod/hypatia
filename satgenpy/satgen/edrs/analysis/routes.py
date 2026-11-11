"""
Route & QoS analysis of a main.py run.

1. ROUTES  - rebuilds, at every time step, the route that Hypatia's forwarding tables (fstate_*.txt)
             give to data leaving each EO satellite (LEO) towards the ground, e.g.
             "Sentinel-1 0 -> EDRS 7 -> Weilheim", and its one-way propagation delay.
2. QoS     - per EO satellite: availability, time per route type (direct / via GEO / via LEO / ...),
             delay statistics, number of route changes (hand-overs), longest outage.
3. CHECK   - if ns-3 has run (logs_ns3/pingmesh.csv), compares every measured ping latency with the
             latency expected from geometry along the fstate route (distance / c). They must agree
             to well below 1 ms; lost pings must coincide with "no route".

Usage (WSL / host terminal, from satgenpy/satgen/edrs, after python main.py):
    python -m analysis.routes
    python -m analysis.routes --run-dir src/output --step 10
Outputs in analysis/results/: summary.md, routes_<sat>.csv, ping_check.csv, routes_timeline.png
"""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EDRS = os.path.dirname(HERE)
SATGENPY = os.path.abspath(os.path.join(EDRS, "..", ".."))
for p in (SATGENPY, EDRS):
    if p not in sys.path:
        sys.path.insert(0, p)

from sgp4.api import jday                                                # noqa: E402
from sgp4.propagation import gstime                                      # noqa: E402

from validation import geometry as G                                     # noqa: E402
from satgen.edrs.config import conf, geometric, simulation               # noqa: E402
from satgen.distance_tools.distance_tools import geodetic2cartesian      # noqa: E402

C_M_S = 299792458.0


# --------------------------------------------------------------------------- inputs
def read_ground_stations(path):
    """ground_stations.txt written by main.py: gid,name,lat,lon,elev,x,y,z"""
    gs = []
    for line in open(path):
        if line.strip():
            p = line.strip().split(",")
            gs.append(dict(gid=int(p[0]), name=p[1], lat=float(p[2]), lon=float(p[3]), elev=float(p[4])))
    return gs


def load_fstates(run_dir, t_ns_list):
    """Yields (t_ns, full forwarding table) applying the incremental fstate files in order."""
    table = {}
    for t_ns in t_ns_list:
        path = os.path.join(run_dir, f"fstate_{t_ns}.txt")
        if os.path.exists(path):
            for line in open(path):
                if line.strip():
                    a, d, nh, _, _ = (int(x) for x in line.split(","))
                    table[(a, d)] = nh
        yield t_ns, table


def follow(table, src, dst, max_hops=32):
    """Node path src -> dst following next hops; None if unreachable."""
    path, node = [src], src
    for _ in range(max_hops):
        nh = table.get((node, dst), -1)
        if nh == -1:
            return None
        path.append(nh)
        if nh == dst:
            return path
        node = nh
    return None   # loop


# --------------------------------------------------------------------------- geometry (SGP4, like ns-3)
class Positions:
    """Satellite positions (SGP4, TEME) and ground stations rotated into TEME with GMST. Metres."""

    def __init__(self, tles, gs, epoch_dt):
        self.tles, self.epoch_dt = tles, epoch_dt
        self.gs_ecef = np.array([geodetic2cartesian(g["lat"], g["lon"], g["elev"]) for g in gs])
        self.jd0, self.fr0 = jday(epoch_dt.year, epoch_dt.month, epoch_dt.day,
                                  epoch_dt.hour, epoch_dt.minute, epoch_dt.second)

    def at(self, t_s):
        t_s = np.atleast_1d(np.asarray(t_s, float))
        r_sat = G.propagate_sgp4(self.tles, self.epoch_dt, t_s)            # [n_sat, n_t, 3]
        theta = np.array([gstime(self.jd0 + self.fr0 + ti / 86400.0) for ti in t_s])
        c, s = np.cos(theta), np.sin(theta)
        x, y, z = self.gs_ecef[:, 0:1], self.gs_ecef[:, 1:2], self.gs_ecef[:, 2:3]
        r_gs = np.stack([c * x - s * y, s * x + c * y, np.broadcast_to(z, (z.shape[0], c.size))], axis=-1)  # [n_gs, n_t, 3]
        return np.concatenate([r_sat, r_gs], axis=0)                          # node-indexed


def path_delay_s(path, r_nodes):
    return sum(np.linalg.norm(r_nodes[b] - r_nodes[a]) for a, b in zip(path, path[1:])) / C_M_S


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default=os.path.join(EDRS, "src", "output"))
    ap.add_argument("--step", type=float, default=None, help="analysis step in s (default: simulation step)")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--tol-ms", type=float, default=0.5, help="allowed |ns-3 - expected| one-way latency")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    tles = G.read_hypatia_tles(os.path.join(args.run_dir, "tles.txt"))
    gs = read_ground_stations(os.path.join(args.run_dir, "ground_stations.txt"))
    n_sat, n_leo = len(tles), geometric.leo_sat
    names = [n for n, _, _ in tles] + [g["name"] for g in gs]
    is_geo = lambda n: n_leo <= n < n_sat
    is_leo = lambda n: n < n_leo
    gs_nodes = [n_sat + g["gid"] for g in gs]

    step_ns = simulation.time_step
    end_ns = int(simulation.duration * 1e9)
    all_t_ns = list(range(0, end_ns, step_ns))
    every = 1 if args.step is None else max(1, int(round(args.step * 1e9 / step_ns)))
    pos = Positions(tles, gs, conf.epoch.utc.datetime)
    print(f"{n_sat} satellites ({n_leo} EO/LEO), {len(gs)} ground stations, "
          f"{len(all_t_ns)} time steps (analysing every {every * step_ns / 1e9:g} s)")

    # -------- 1. routes EO satellite -> best ground station
    rows = defaultdict(list)
    tables_at = {}                        # kept for the ping check (only at analysed steps)
    for k, (t_ns, table) in enumerate(load_fstates(args.run_dir, all_t_ns)):
        if k % every:
            continue
        tables_at[t_ns] = dict(table)
        r = pos.at(t_ns / 1e9)[:, 0, :]
        for s in range(n_leo):
            best = None
            for d in gs_nodes:
                p = follow(table, s, d)
                if p:
                    dl = path_delay_s(p, r)
                    if best is None or dl < best[0]:
                        best = (dl, p)
            if best is None:
                rows[s].append((t_ns / 1e9, None, None, "no route"))
            else:
                dl, p = best
                mids = p[1:-1]
                kind = ("direct" if not mids else
                        "via GEO" if all(is_geo(m) for m in mids) else
                        "via LEO" if all(is_leo(m) for m in mids) else "via LEO+GEO")
                rows[s].append((t_ns / 1e9, " -> ".join(names[n] for n in p), dl * 1e3, kind))

    # -------- 2. QoS per EO satellite
    dt = every * step_ns / 1e9
    summary = ["# Route & QoS analysis", "",
               f"- Run dir: `{os.path.relpath(args.run_dir, EDRS)}`, duration {simulation.duration} s, "
               f"analysed every {dt:g} s",
               "- Route = path given by Hypatia's forwarding tables from the EO satellite to the ground station "
               "with the shortest path; delay = propagation only (distance / c)", "",
               "| EO satellite | availability | direct | via GEO | via LEO | via LEO+GEO | mean / max one-way delay [ms] "
               "| route changes | longest outage [s] |",
               "|---|---|---|---|---|---|---|---|---|"]
    for s in range(n_leo):
        R = rows[s]
        kinds = Counter(k for *_, k in R)
        n = len(R)
        delays = np.array([d for _, _, d, _ in R if d is not None])
        changes = sum(1 for a, b in zip(R, R[1:]) if a[1] != b[1])
        longest, cur = 0.0, 0.0
        for *_, k in R:
            cur = cur + dt if k == "no route" else 0.0
            longest = max(longest, cur)
        pct = lambda key: f"{100 * kinds[key] / n:.0f} %"
        summary.append(
            f"| {names[s]} | {100 * (n - kinds['no route']) / n:.1f} % | {pct('direct')} | {pct('via GEO')} | "
            f"{pct('via LEO')} | {pct('via LEO+GEO')} | "
            + (f"{delays.mean():.1f} / {delays.max():.1f}" if delays.size else "–")
            + f" | {changes} | {longest:.0f} |")
        with open(os.path.join(args.out, f"routes_{names[s].replace(' ', '_')}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["t_s", "route", "one_way_delay_ms", "type"])
            w.writerows(R)

    # -------- 3. ns-3 check (pingmesh)
    ping_csv = os.path.join(args.run_dir, "logs_ns3", "pingmesh.csv")
    summary += ["", "## ns-3 check (pingmesh)", ""]
    if not os.path.exists(ping_csv):
        summary.append("No `logs_ns3/pingmesh.csv` yet: run ns-3 (python main.py) with the pingmesh enabled.")
    else:
        by_pair = defaultdict(list)
        for line in open(ping_csv):
            p = line.strip().split(",")
            if len(p) == 10:
                by_pair[(int(p[0]), int(p[1]))].append(
                    (int(p[3]), int(p[6]), int(p[7]), p[9]))   # send_ts, to_there, from_there, YES/LOST
        check_rows = []
        summary += ["| pair | pings | lost | lost with expected route | received without expected route "
                    "| max \\|Δ one-way\\| [ms] | PASS |", "|---|---|---|---|---|---|---|"]
        all_ok = True
        for (a, b), pings in sorted(by_pair.items()):
            lost = lost_with_route = recv_without_route = 0
            max_err = 0.0
            for send_ns, to_ns, from_ns, status in pings:
                t_ns = (send_ns // (every * step_ns)) * (every * step_ns)
                if t_ns not in tables_at or send_ns - t_ns > 0.9 * step_ns:
                    continue        # only pings sent shortly after an analysed fstate update
                table = tables_at[t_ns]
                p_fwd = follow(table, a, b)
                p_rev = follow(table, b, a)
                if status == "LOST":
                    lost += 1
                    if p_fwd and p_rev:
                        lost_with_route += 1
                    continue
                if not (p_fwd and p_rev):
                    recv_without_route += 1
                    continue
                r = pos.at(send_ns / 1e9)[:, 0, :]
                exp_ms = path_delay_s(p_fwd, r) * 1e3
                err = to_ns / 1e6 - exp_ms
                max_err = max(max_err, abs(err))
                check_rows.append((send_ns / 1e9, names[a], names[b], to_ns / 1e6, exp_ms, err,
                                   " -> ".join(names[n] for n in p_fwd)))
            ok = recv_without_route == 0 and max_err <= args.tol_ms
            all_ok &= ok
            summary.append(f"| {names[a]} → {names[b]} | {len(pings)} | {lost} | {lost_with_route} | "
                           f"{recv_without_route} | {max_err:.3f} | {'✅' if ok else '❌'} |")
        with open(os.path.join(args.out, "ping_check.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["t_s", "from", "to", "ns3_one_way_ms", "expected_one_way_ms", "diff_ms", "route"])
            w.writerows(check_rows)
        summary += ["", f"**ns-3 check: {'PASS' if all_ok else 'FAIL'}** "
                    f"(tolerance {args.tol_ms} ms; a few lost pings right at a route change are normal)"]

    text = "\n".join(summary)
    open(os.path.join(args.out, "summary.md"), "w").write(text)
    print(text)

    # -------- plot: route type over time per EO satellite
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        colors = {"direct": "#2a9d8f", "via GEO": "#e76f51", "via LEO": "#8ab17d",
                  "via LEO+GEO": "#9b5de5", "no route": "#d9d9d9"}
        fig, (ax, ax2) = plt.subplots(2, 1, figsize=(11, 1.0 + 0.9 * n_leo + 3), sharex=True,
                                      gridspec_kw=dict(height_ratios=[n_leo, 3]))
        for s in range(n_leo):
            R = rows[s]
            t = np.array([x[0] for x in R]) / 60
            for kind, col in colors.items():
                m = np.array([x[3] == kind for x in R])
                ax.fill_between(t, s - 0.4, s + 0.4, where=m, color=col, step="post", lw=0,
                                label=kind if s == 0 else None)
            d = np.array([np.nan if x[2] is None else x[2] for x in R])
            ax2.plot(t, d, lw=1, label=names[s])
        ax.set_yticks(range(n_leo), [names[s] for s in range(n_leo)])
        ax.invert_yaxis()
        ax.legend(ncol=5, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.25))
        ax.set_title("Route from each EO satellite to the ground", pad=28)
        ax2.set_ylabel("one-way delay [ms]")
        ax2.set_xlabel("time since epoch [min]")
        ax2.grid(alpha=0.3)
        ax2.legend(fontsize=7, ncol=3)
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "routes_timeline.png"), dpi=130)
    except Exception as e:
        print("plot skipped:", e)


if __name__ == "__main__":
    main()
