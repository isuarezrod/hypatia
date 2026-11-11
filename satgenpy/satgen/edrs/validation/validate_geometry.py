"""
Step 2 - Validate geometry.

Compares, pair by pair, the contact windows and ranges produced by the Hypatia
pipeline (ephem propagation + config.radius() + src/visibility.py +
distance_tools) against:

  [ref]      independent propagation: python-sgp4 + validation/geometry.py
  [gmat]     optional: GMAT ReportFile (numerical propagation, see make_gmat_script.py)
  [baseline] optional: external contact windows (e.g. ISRO baseline from the paper)

Usage (from satgenpy/satgen/edrs):
    python -m validation.validate_geometry --duration 86400 --step 10
    python -m validation.validate_geometry --duration 86400 --step 10 \
        --gmat validation/gmat/report_edrs.txt --baseline validation/isro_baseline.csv

Outputs (in --out, default validation/results):
    windows_<source>.csv, window_comparison.csv, ranges.csv, report.md, geometry_validation.png
"""
import argparse
import csv
import os
import sys
from datetime import timedelta

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EDRS = os.path.dirname(HERE)
SATGENPY = os.path.abspath(os.path.join(EDRS, "..", ".."))
for p in (SATGENPY, EDRS):
    if p not in sys.path:
        sys.path.insert(0, p)

from astropy import units as u                                         # noqa: E402

from validation import geometry as G                                   # noqa: E402
from satgen.edrs.config import conf, geometric, visibility, radius     # noqa: E402
from satgen.edrs.src.visibility import contact_leo_geo, contact_leo_leo  # noqa: E402
from satgen.distance_tools.distance_tools import distance_m_between_satellites  # noqa: E402
from satgen.tles.read_tles import read_tles                            # noqa: E402
from src.simple_isls import generate_leo_geo_isls, generate_leo_leo_isls  # noqa: E402


# --------------------------------------------------------------------------- helpers
def pair_kind(a, b, n_leo):
    return "LEO-LEO" if (a < n_leo and b < n_leo) else "LEO-GEO"


def ref_cfg():
    """Reference criteria = same numbers as config.py, converted to SI."""
    return G.VisCfg(
        re_m=geometric.Re * 1e3,
        h_graze_m=geometric.ATM * 1e3,
        theta_leo_max_deg=visibility.theta_leo_max,
        theta_geo_max_deg=visibility.theta_geo_max,
        max_range_m=visibility.max_range * 1e3,
    )


def ref_visibility(kind, ra, rb, cfg):
    if kind == "LEO-GEO":
        return G.visible_leo_geo(ra, rb, cfg)
    return G.visible_leo_leo(ra, rb, cfg)


def read_gmat_report(path, n_sat):
    """
    GMAT ReportFile written by make_gmat_script.py:
    header line, then rows: ElapsedSecs X0 Y0 Z0 X1 Y1 Z1 ...  (km, EarthMJ2000Eq)
    returns t [n_t], r [n_sat, n_t, 3] in metres
    """
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue                       # header
    arr = np.array(rows)
    t = arr[:, 0]
    r = arr[:, 1:1 + 3 * n_sat].reshape(len(t), n_sat, 3).transpose(1, 0, 2) * 1e3
    return t, r


def read_baseline(path):
    """CSV with columns: pair (e.g. '0-2'), start_s, end_s  [optional: range_min_km, range_max_km]"""
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            a, b = (int(x) for x in row["pair"].split("-"))
            out.setdefault((a, b), []).append((float(row["start_s"]), float(row["end_s"])))
    return out


def summarize(rows, tol_edge_s):
    m = [r for r in rows if r["status"] == "matched"]
    missed = sum(r["status"] == "MISSED" for r in rows)
    spurious = sum(r["status"] == "SPURIOUS" for r in rows)
    if m:
        ds = np.array([r["d_start"] for r in m])
        de = np.array([r["d_end"] for r in m])
        worst = float(np.nanmax(np.abs(np.r_[ds, de])))
    else:
        worst = np.nan
    ok = missed == 0 and spurious == 0 and (not m or worst <= tol_edge_s)
    return dict(n_ref=len([r for r in rows if r["ref_idx"] >= 0]), n_matched=len(m),
                missed=missed, spurious=spurious, worst_edge_s=worst, ok=ok)


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tle", default=os.path.join(EDRS, "src", "tles", "tle.txt"))
    ap.add_argument("--duration", type=float, default=86400.0, help="s (use >= 1 day to see several passes)")
    ap.add_argument("--step", type=float, default=10.0, help="s (Hypatia sampling step for this check)")
    ap.add_argument("--tol-range-km", type=float, default=1.0, help="Hypatia vs SGP4 range tolerance")
    ap.add_argument("--tol-edge-s", type=float, default=None, help="default: = step")
    ap.add_argument("--tol-range-km-gmat", type=float, default=10.0)
    ap.add_argument("--tol-edge-s-gmat", type=float, default=30.0)
    ap.add_argument("--gmat", default=None, help="GMAT ReportFile (see make_gmat_script.py)")
    ap.add_argument("--baseline", default=None, help="external windows CSV (e.g. ISRO baseline)")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    args = ap.parse_args()
    tol_edge = args.tol_edge_s if args.tol_edge_s is not None else args.step
    os.makedirs(args.out, exist_ok=True)

    # ---- topology
    tles = G.read_hypatia_tles(args.tle)
    n_sat = len(tles)
    n_leo = geometric.leo_sat
    leo_ids, geo_ids = list(range(n_leo)), list(range(n_leo, n_sat))
    pairs = generate_leo_geo_isls(leo_ids, geo_ids) + generate_leo_leo_isls(leo_ids)
    cfg = ref_cfg()
    t = np.arange(0.0, args.duration + 1e-9, args.step)
    epoch = conf.epoch
    epoch_dt = epoch.utc.datetime
    print(f"{n_sat} satellites ({n_leo} LEO), {len(pairs)} pairs, {len(t)} samples, epoch {epoch.utc.iso} UTC")

    # ---- [ref] SGP4
    r_ref = G.propagate_sgp4(tles, epoch_dt, t)

    def ref_scalar(kind, a, b):
        def f(ts):
            r = G.propagate_sgp4([tles[a], tles[b]], epoch_dt, [ts])
            return bool(ref_visibility(kind, r[0, 0], r[1, 0], cfg)[0])
        return f

    # ---- [hypatia] ephem + radius() + visibility.py + distance_tools
    sats = read_tles(args.tle)["satellites"]
    eps = str(epoch)
    vis_hyp = {p: np.zeros(len(t), bool) for p in pairs}
    rng_hyp = {p: np.zeros(len(t)) for p in pairs}
    print("Running Hypatia geometry path (ephem)...")
    for k, ts in enumerate(t):
        ds = str(epoch + ts * u.s)
        r = [radius(s, eps, ds) for s in sats]
        for (a, b) in pairs:
            kind = pair_kind(a, b, n_leo)
            vis_hyp[(a, b)][k] = contact_leo_geo(r[a], r[b]) if kind == "LEO-GEO" else contact_leo_leo(r[a], r[b])
            rng_hyp[(a, b)][k] = distance_m_between_satellites(sats[a], sats[b], eps, ds)

    # ---- [gmat]
    gmat = None
    if args.gmat:
        tg, rg = read_gmat_report(args.gmat, n_sat)
        gmat = (tg, rg)

    baseline = read_baseline(args.baseline) if args.baseline else None

    # ---- compare
    win_rows, cmp_rows, summary = [], [], []
    rng_csv = open(os.path.join(args.out, "ranges.csv"), "w", newline="")
    rw = csv.writer(rng_csv)
    rw.writerow(["t_s", "pair", "range_ref_km", "range_hypatia_km", "vis_ref", "vis_hypatia"])

    plot_data = {}
    for (a, b) in pairs:
        kind = pair_kind(a, b, n_leo)
        name = f"{a}-{b}"
        m_ref, diag = ref_visibility(kind, r_ref[a], r_ref[b], cfg)
        rng_ref = diag["range_m"]
        w_ref = G.refine_edges(G.windows_from_mask(t, m_ref), t, ref_scalar(kind, a, b))
        w_hyp = G.windows_from_mask(t, vis_hyp[(a, b)])
        sources = {"ref": w_ref, "hypatia": w_hyp}

        for k in range(len(t)):
            rw.writerow([t[k], name, rng_ref[k] / 1e3, rng_hyp[(a, b)][k] / 1e3, int(m_ref[k]), int(vis_hyp[(a, b)][k])])

        # range error over the whole interval (independent of visibility)
        drng = np.abs(rng_hyp[(a, b)] - rng_ref) / 1e3
        rows = G.match_windows(w_ref, w_hyp)
        s = summarize(rows, tol_edge)
        s.update(pair=name, kind=kind, test="hypatia vs sgp4",
                 max_range_err_km=float(drng.max()), range_ok=bool(drng.max() <= args.tol_range_km))
        summary.append(s)
        for r_ in rows:
            cmp_rows.append(dict(pair=name, test="hypatia", **r_))

        if gmat is not None:
            tg, rg = gmat
            m_g, dg = ref_visibility(kind, rg[a], rg[b], cfg)
            w_g = G.windows_from_mask(tg, m_g)
            sources["gmat"] = w_g
            r_ref_g = G.propagate_sgp4(tles, epoch_dt, tg)
            drg = np.abs(dg["range_m"] - G._norm(r_ref_g[b] - r_ref_g[a])) / 1e3
            rows = G.match_windows(w_ref, w_g)
            s = summarize(rows, args.tol_edge_s_gmat)
            s.update(pair=name, kind=kind, test="gmat vs sgp4",
                     max_range_err_km=float(drg.max()), range_ok=bool(drg.max() <= args.tol_range_km_gmat))
            summary.append(s)
            for r_ in rows:
                cmp_rows.append(dict(pair=name, test="gmat", **r_))

        if baseline is not None and (a, b) in baseline:
            w_b = baseline[(a, b)]
            sources["baseline"] = w_b
            rows = G.match_windows(w_ref, w_b)
            s = summarize(rows, args.tol_edge_s_gmat)
            s.update(pair=name, kind=kind, test="baseline vs sgp4", max_range_err_km=np.nan, range_ok=True)
            summary.append(s)
            for r_ in rows:
                cmp_rows.append(dict(pair=name, test="baseline", **r_))

        for src, ws in sources.items():
            for (s0, s1) in ws:
                i0, i1 = np.searchsorted(t, s0), np.searchsorted(t, s1, side="right")
                rr = rng_ref[i0:i1] / 1e3 if i1 > i0 else np.array([np.nan])
                win_rows.append(dict(source=src, pair=name, kind=kind, start_s=s0, end_s=s1,
                                     duration_s=s1 - s0, range_min_km=np.nanmin(rr), range_max_km=np.nanmax(rr),
                                     start_utc=(epoch_dt + timedelta(seconds=float(s0))).isoformat(),
                                     end_utc=(epoch_dt + timedelta(seconds=float(s1))).isoformat()))
        plot_data[name] = (rng_ref, rng_hyp[(a, b)], m_ref, vis_hyp[(a, b)], diag)
    rng_csv.close()

    # ---- write CSVs
    def dump(fname, rows):
        if not rows:
            return
        with open(os.path.join(args.out, fname), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    dump("windows.csv", win_rows)
    dump("window_comparison.csv", cmp_rows)

    # ---- report
    lines = ["# Geometry validation report", "",
             f"- TLE file: `{os.path.relpath(args.tle, EDRS)}`",
             f"- Epoch: {epoch.utc.iso} UTC, duration {args.duration:.0f} s, step {args.step:g} s",
             f"- Criteria (from config.py, SI): h_graze = {cfg.h_graze_m/1e3:.0f} km, "
             f"theta_LEO(zenith) < {cfg.theta_leo_max_deg} deg, theta_GEO(nadir) < {cfg.theta_geo_max_deg} deg, "
             f"max range = {cfg.max_range_m/1e3:.0f} km",
             f"- Tolerances: range {args.tol_range_km} km / edge {tol_edge} s (Hypatia); "
             f"range {args.tol_range_km_gmat} km / edge {args.tol_edge_s_gmat} s (GMAT, baseline)", "",
             "| pair | kind | test | ref windows | matched | missed | spurious | worst edge Δ [s] | max range Δ [km] | PASS |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    all_ok = True
    for s in summary:
        ok = s["ok"] and s["range_ok"]
        all_ok &= ok
        lines.append(f"| {s['pair']} | {s['kind']} | {s['test']} | {s['n_ref']} | {s['n_matched']} | {s['missed']} | "
                     f"{s['spurious']} | {s['worst_edge_s']:.2f} | {s['max_range_err_km']:.3f} | {'✅' if ok else '❌'} |")
    lines += ["", f"**Overall: {'PASS' if all_ok else 'FAIL'}**", ""]
    lines += ["## Reference windows (SGP4, edges refined to 0.01 s)", "",
              "| pair | start UTC | end UTC | duration [min] | range [km] |", "|---|---|---|---|---|"]
    for w in win_rows:
        if w["source"] == "ref":
            lines.append(f"| {w['pair']} | {w['start_utc'][:19]} | {w['end_utc'][:19]} | {w['duration_s']/60:.1f} | "
                         f"{w['range_min_km']:.0f} – {w['range_max_km']:.0f} |")
    report = "\n".join(lines)
    open(os.path.join(args.out, "report.md"), "w").write(report)
    print(report)

    # ---- plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = len(plot_data)
        fig, axs = plt.subplots(n, 1, figsize=(11, 2.3 * n), sharex=True, squeeze=False)
        th = t / 3600
        for ax, (name, (rr, rh, mr, mh, dg)) in zip(axs[:, 0], plot_data.items()):
            ax.plot(th, rr / 1e3, color="#2a6f97", lw=1.2, label="range SGP4")
            ax.plot(th, rh / 1e3, color="#e07a1f", lw=0.8, ls="--", label="range Hypatia (ephem)")
            ymax = np.nanmax(rr) / 1e3
            ax.fill_between(th, 0, ymax, where=mr, color="#2a6f97", alpha=0.15, step="mid", label="visible (SGP4 ref)")
            ax.fill_between(th, 0, ymax * 0.5, where=mh, color="#e07a1f", alpha=0.35, step="mid", label="visible (Hypatia)")
            ax.set_ylabel(f"{name}\nrange [km]")
            ax.grid(alpha=0.3)
        axs[0, 0].legend(loc="upper right", fontsize=8, ncol=4)
        axs[-1, 0].set_xlabel("time since epoch [h]")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "geometry_validation.png"), dpi=130)
    except Exception as e:  # plotting is optional
        print("plot skipped:", e)

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
