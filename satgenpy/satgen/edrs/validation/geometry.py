"""
Independent reference geometry for validating visibility.py.

Deliberately does NOT import anything from satgen / ephem:
  - propagation : python-sgp4 (Vallado's SGP4/SDP4), positions in TEME [m]
  - geometry    : own implementation, everything in SI (metres, radians)

Distances and all visibility tests are invariant under a common rotation of
both position vectors, so TEME, ECEF or GCRS give the same answer as long as
both satellites are expressed in the same frame at the same instant.
"""
from dataclasses import dataclass
import numpy as np
from sgp4.api import Satrec, SatrecArray, jday

RE_M = 6378137.0  # WGS-84 equatorial radius [m]


# --------------------------------------------------------------------------- TLE / propagation
def read_hypatia_tles(path):
    """Hypatia tle.txt: header line '<n_orbits> <n_sats_per_orbit>', then name/L1/L2 triplets."""
    lines = [l.rstrip("\n") for l in open(path) if l.strip()]
    body = lines[1:]
    out = []
    for i in range(0, len(body), 3):
        out.append((body[i].strip(), body[i + 1], body[i + 2]))
    return out


def propagate_sgp4(tles, epoch_utc, t_s):
    """
    tles      : list of (name, l1, l2)
    epoch_utc : datetime (UTC) of t = 0
    t_s       : array of seconds since epoch
    returns   : r [n_sat, n_t, 3] in metres (TEME)
    """
    sats = SatrecArray([Satrec.twoline2rv(l1, l2) for _, l1, l2 in tles])
    jd0, fr0 = jday(epoch_utc.year, epoch_utc.month, epoch_utc.day,
                    epoch_utc.hour, epoch_utc.minute,
                    epoch_utc.second + epoch_utc.microsecond * 1e-6)
    t_s = np.atleast_1d(np.asarray(t_s, dtype=float))
    jd = np.full_like(t_s, jd0)
    fr = fr0 + t_s / 86400.0
    err, r, _ = sats.sgp4(jd, fr)
    if np.any(err):
        bad = np.argwhere(err)
        raise RuntimeError(f"SGP4 error codes at (sat, t) = {bad[:5].tolist()} ...")
    return r * 1e3


# --------------------------------------------------------------------------- geometry
@dataclass
class VisCfg:
    re_m: float = RE_M
    h_graze_m: float = 100e3        # grazing altitude for Earth/atmosphere blockage
    theta_leo_max_deg: float = 90.0  # LEO terminal half-cone around local ZENITH
    theta_geo_max_deg: float = 90.0  # GEO terminal half-cone around local NADIR
    max_range_m: float = np.inf


def _norm(v):
    return np.linalg.norm(v, axis=-1)


def _angle(a, b):
    c = np.sum(a * b, axis=-1) / (_norm(a) * _norm(b))
    return np.arccos(np.clip(c, -1.0, 1.0))


def los_min_distance(r1, r2):
    """Minimum distance from Earth's centre to the SEGMENT r1-r2 (vectorised, metres)."""
    d = r2 - r1                                     # r1 -> r2
    t = -np.sum(r1 * d, axis=-1) / np.sum(d * d, axis=-1)
    t = np.clip(t, 0.0, 1.0)[..., None]
    return _norm(r1 + t * d)


def visible_leo_geo(r_leo, r_geo, cfg: VisCfg):
    """Returns (mask, diagnostics dict). Arrays [..., 3] in metres."""
    d_lg = r_geo - r_leo                            # LEO -> GEO
    rng = _norm(d_lg)
    hmin = los_min_distance(r_leo, r_geo) - cfg.re_m
    th_leo = np.degrees(_angle(d_lg, r_leo))        # off-zenith at LEO
    th_geo = np.degrees(_angle(-d_lg, -r_geo))      # off-nadir at GEO
    mask = ((hmin > cfg.h_graze_m)
            & (th_leo < cfg.theta_leo_max_deg)
            & (th_geo < cfg.theta_geo_max_deg)
            & (rng <= cfg.max_range_m))
    return mask, dict(range_m=rng, hmin_m=hmin, th_leo_deg=th_leo, th_geo_deg=th_geo)


def visible_leo_leo(r_a, r_b, cfg: VisCfg):
    rng = _norm(r_b - r_a)
    hmin = los_min_distance(r_a, r_b) - cfg.re_m
    mask = (hmin > cfg.h_graze_m) & (rng <= cfg.max_range_m)
    return mask, dict(range_m=rng, hmin_m=hmin)


# --------------------------------------------------------------------------- windows
def windows_from_mask(t_s, mask):
    """Contact windows [(start, end)] from a sampled boolean series (edges at sample times)."""
    m = np.asarray(mask, dtype=bool)
    if m.size == 0:
        return []
    dm = np.diff(m.astype(int))
    starts = list(t_s[1:][dm == 1])
    ends = list(t_s[:-1][dm == -1])
    if m[0]:
        starts.insert(0, t_s[0])
    if m[-1]:
        ends.append(t_s[-1])
    return list(zip(starts, ends))


def refine_edges(windows, t_s, f_vis, tol_s=0.01):
    """
    Bisection refinement of window edges with a scalar visibility function f_vis(t)->bool.
    Edges clipped at the simulation boundaries are left untouched.
    """
    dt = t_s[1] - t_s[0] if len(t_s) > 1 else 0.0
    out = []
    for a, b in windows:
        if a > t_s[0]:
            lo, hi = a - dt, a              # invisible at lo, visible at hi
            while hi - lo > tol_s:
                mid = 0.5 * (lo + hi)
                lo, hi = (lo, mid) if f_vis(mid) else (mid, hi)
            a = hi
        if b < t_s[-1]:
            lo, hi = b, b + dt              # visible at lo, invisible at hi
            while hi - lo > tol_s:
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if f_vis(mid) else (lo, mid)
            b = lo
        out.append((a, b))
    return out


def match_windows(ref, test):
    """
    Pair each reference window with the overlapping test window (largest overlap).
    Returns list of dicts with start/end deltas (test - ref); unmatched windows flagged.
    """
    rows, used = [], set()
    for i, (a, b) in enumerate(ref):
        best, best_ov = None, 0.0
        for j, (c, d) in enumerate(test):
            ov = min(b, d) - max(a, c)
            if ov > best_ov:
                best, best_ov = j, ov
        if best is None:
            rows.append(dict(ref_idx=i, ref_start=a, ref_end=b, test_start=np.nan,
                             test_end=np.nan, d_start=np.nan, d_end=np.nan, status="MISSED"))
        else:
            used.add(best)
            c, d = test[best]
            rows.append(dict(ref_idx=i, ref_start=a, ref_end=b, test_start=c, test_end=d,
                             d_start=c - a, d_end=d - b, status="matched"))
    for j, (c, d) in enumerate(test):
        if j not in used:
            rows.append(dict(ref_idx=-1, ref_start=np.nan, ref_end=np.nan, test_start=c,
                             test_end=d, d_start=np.nan, d_end=np.nan, status="SPURIOUS"))
    return rows
