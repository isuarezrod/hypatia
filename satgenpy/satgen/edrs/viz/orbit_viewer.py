"""
3D orbit viewer: Earth + satellites propagated from tle.txt, in a single HTML file.

Usage (WSL / host terminal, from satgenpy/satgen/edrs):
    python -m viz.orbit_viewer                         # 24 h, 30 s step
    python -m viz.orbit_viewer --duration 21600 --step 10
    explorer.exe viz/output/orbit_viewer.html          # open it in your Windows browser

What is drawn
  - Earth (textured, rotating with GMST) in an inertial frame (TEME ~ ECI)
  - every satellite of tle.txt, with its orbit trail
  - LEO-GEO / LEO-LEO lines when the link is geometrically feasible
    (same criteria as src/visibility.py, computed with validation/geometry.py)

Hook for the ns-3 links (later)
  The HTML reads one JSON object, DATA, with:
    DATA.links = [{"a": 0, "b": 2, "kind": "LEO-GEO", "source": "geometry",
                   "windows": [[t_start_s, t_end_s], ...]}, ...]
  To show the ISLs that ns-3 actually used, produce the same structure from the ns-3
  output (e.g. isl_utilization.csv / fstate_*.txt) with source="ns3" and pass it with
  --links-json my_ns3_links.json. Links with source="ns3" are drawn in a different colour.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EDRS = os.path.dirname(HERE)
SATGENPY = os.path.abspath(os.path.join(EDRS, "..", ".."))
for p in (SATGENPY, EDRS):
    if p not in sys.path:
        sys.path.insert(0, p)

from astropy.utils import iers                                          # noqa: E402
iers.conf.auto_download = False

from validation import geometry as G                                    # noqa: E402
from satgen.edrs.config import conf, geometric, visibility              # noqa: E402

EARTH_ROT_RATE = 7.2921150e-5  # rad/s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tle", default=os.path.join(EDRS, "src", "tles", "tle.txt"))
    ap.add_argument("--duration", type=float, default=86400.0)
    ap.add_argument("--step", type=float, default=30.0)
    ap.add_argument("--links-json", default=None, help="extra links (e.g. from ns-3), same format as DATA.links")
    ap.add_argument("--no-geometry-links", action="store_true")
    ap.add_argument("--out", default=os.path.join(HERE, "output", "orbit_viewer.html"))
    args = ap.parse_args()

    tles = G.read_hypatia_tles(args.tle)
    n_leo = geometric.leo_sat
    epoch = conf.epoch.utc
    t = np.arange(0.0, args.duration + 1e-9, args.step)
    r = G.propagate_sgp4(tles, epoch.datetime, t) / 1e3          # km, TEME
    gmst0 = float(epoch.sidereal_time("mean", "greenwich").rad)

    links = []
    if not args.no_geometry_links:
        cfg = G.VisCfg(re_m=geometric.Re * 1e3, h_graze_m=geometric.ATM * 1e3,
                       theta_leo_max_deg=visibility.theta_leo_max,
                       theta_geo_max_deg=visibility.theta_geo_max,
                       max_range_m=visibility.max_range * 1e3)
        for a in range(len(tles)):
            for b in range(a + 1, len(tles)):
                if a >= n_leo:                     # GEO-GEO not modelled
                    continue
                kind = "LEO-LEO" if b < n_leo else "LEO-GEO"
                ra, rb = r[a] * 1e3, r[b] * 1e3
                m = (G.visible_leo_leo(ra, rb, cfg) if kind == "LEO-LEO" else G.visible_leo_geo(ra, rb, cfg))[0]
                links.append(dict(a=a, b=b, kind=kind, source="geometry",
                                  windows=[[float(s), float(e)] for s, e in G.windows_from_mask(t, m)]))
    if args.links_json:
        links += json.load(open(args.links_json))

    data = dict(
        epoch_utc=epoch.isot,
        t=t.tolist(),
        gmst0=gmst0,
        omega_e=EARTH_ROT_RATE,
        re_km=geometric.Re,
        sats=[dict(name=name, kind="LEO" if i < n_leo else "GEO",
                   pos=np.round(r[i], 2).tolist()) for i, (name, _, _) in enumerate(tles)],
        links=links,
    )
    html = TEMPLATE.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    open(args.out, "w").write(html)
    print("Viewer written to", args.out)
    print("Open it in a browser (from WSL: explorer.exe " + os.path.relpath(args.out) + ")")


TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>EDRS / Sentinel orbit viewer</title>
<style>
  html, body { margin: 0; height: 100%; background: #05070d; color: #dde3ee;
               font: 13px/1.4 system-ui, -apple-system, Segoe UI, sans-serif; overflow: hidden; }
  #panel { position: absolute; top: 12px; left: 12px; background: rgba(12,16,26,.85);
           border: 1px solid #263046; border-radius: 8px; padding: 12px 14px; width: 300px; }
  #panel h1 { font-size: 14px; margin: 0 0 8px; }
  #time { font-variant-numeric: tabular-nums; font-size: 15px; margin-bottom: 8px; }
  .row { display: flex; align-items: center; gap: 8px; margin: 6px 0; }
  .row label { flex: 1; }
  input[type=range] { width: 100%; }
  button { background: #1f2a40; color: #dde3ee; border: 1px solid #34425f; border-radius: 5px;
           padding: 4px 10px; cursor: pointer; }
  #legend { margin-top: 8px; border-top: 1px solid #263046; padding-top: 8px; }
  .sw { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 6px; vertical-align: -1px; }
  #links { margin-top: 6px; max-height: 120px; overflow: auto; font-variant-numeric: tabular-nums; }
  .lbl { position: absolute; color: #fff; font-size: 11px; pointer-events: none;
         text-shadow: 0 0 3px #000, 0 0 3px #000; white-space: nowrap; }
  #hint { position: absolute; bottom: 10px; left: 12px; color: #7d8aa5; font-size: 11px; }
</style>
<script type="importmap">
{ "imports": {
  "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
  "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/" } }
</script>
</head>
<body>
<div id="panel">
  <h1>Orbit viewer</h1>
  <div id="time">–</div>
  <div class="row"><button id="play">Pause</button>
    <label>speed <span id="spdv"></span>×</label></div>
  <div class="row"><input id="spd" type="range" min="0" max="4" step="0.01" value="2.48"></div>
  <div class="row"><input id="scrub" type="range" min="0" max="1" step="0.0001" value="0"></div>
  <div class="row"><label><input id="trails" type="checkbox" checked> orbit trails</label>
                   <label><input id="showlinks" type="checkbox" checked> links</label></div>
  <div class="row"><label><input id="labels" type="checkbox" checked> labels</label>
                   <label><input id="ecef" type="checkbox"> Earth-fixed view</label></div>
  <div id="legend"></div>
  <div id="links"></div>
</div>
<div id="hint">drag = rotate · scroll = zoom · right-drag = pan</div>
<script type="module">
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

const DATA = /*__DATA__*/null;
const S = 1 / 1000;                       // scene units: 1 = 1000 km
const T = DATA.t, N = T.length, T_END = T[N - 1];

// ECI/TEME (x, y, z) -> three.js (x, z, -y)  [three.js is y-up]
const toScene = (p, out) => out.set(p[0] * S, p[2] * S, -p[1] * S);

// ---------------------------------------------------------------- scene
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(innerWidth, innerHeight);
document.body.appendChild(renderer.domElement);
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(45, innerWidth / innerHeight, 0.1, 1000);
camera.position.set(45, 25, 65);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const sun = new THREE.DirectionalLight(0xffffff, 1.4); sun.position.set(100, 20, 60); scene.add(sun);

// stars
{ const g = new THREE.BufferGeometry(), n = 2500, a = new Float32Array(3 * n);
  for (let i = 0; i < n; i++) { const v = new THREE.Vector3().randomDirection().multiplyScalar(400 + 200 * Math.random()); a.set([v.x, v.y, v.z], 3 * i); }
  g.setAttribute("position", new THREE.BufferAttribute(a, 3));
  scene.add(new THREE.Points(g, new THREE.PointsMaterial({ color: 0x8892a8, size: 0.6 }))); }

// Earth (texture lon 0 at +x, lon +90E at -z, matching toScene of ECEF)
const earth = new THREE.Group(); scene.add(earth);
const RE = DATA.re_km * S;
const earthMat = new THREE.MeshPhongMaterial({ color: 0x2b5d8a, shininess: 8 });
earth.add(new THREE.Mesh(new THREE.SphereGeometry(RE, 96, 64), earthMat));
new THREE.TextureLoader().load(
  "https://cdn.jsdelivr.net/npm/three-globe@2.31.0/example/img/earth-blue-marble.jpg",
  tex => { tex.colorSpace = THREE.SRGBColorSpace; earthMat.map = tex; earthMat.color.set(0xffffff); earthMat.needsUpdate = true; });
// equator + prime meridian helpers
{ const eq = new THREE.Mesh(new THREE.TorusGeometry(RE * 1.001, 0.012, 6, 128), new THREE.MeshBasicMaterial({ color: 0xffd166 }));
  eq.rotation.x = Math.PI / 2; earth.add(eq);
  const axis = new THREE.Mesh(new THREE.CylinderGeometry(0.03, 0.03, RE * 2.6, 8), new THREE.MeshBasicMaterial({ color: 0x9aa7c2 }));
  earth.add(axis); }

// ---------------------------------------------------------------- satellites
const COLORS = { LEO: [0x4cc9f0, 0x80ffdb, 0xa0c4ff, 0x72efdd], GEO: [0xf72585, 0xff9e00, 0xffd60a, 0xe76f51] };
const cnt = { LEO: 0, GEO: 0 };
const sats = DATA.sats.map(s => {
  const color = COLORS[s.kind][cnt[s.kind]++ % 4];
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(s.kind === "GEO" ? 0.55 : 0.35, 16, 12), new THREE.MeshBasicMaterial({ color }));
  scene.add(mesh);
  // full inertial trail
  const pts = s.pos.map(p => toScene(p, new THREE.Vector3()));
  const trail = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts),
                               new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.35 }));
  scene.add(trail);
  const lbl = document.createElement("div"); lbl.className = "lbl"; lbl.textContent = s.name; document.body.appendChild(lbl);
  return { ...s, color, mesh, trail, lbl, pts };
});
document.getElementById("legend").innerHTML = sats.map(s =>
  `<div><span class="sw" style="background:#${s.color.toString(16).padStart(6, "0")}"></span>${s.name} (${s.kind})</div>`).join("");

// ---------------------------------------------------------------- links
const LINK_COL = { geometry: 0x7CFC00, ns3: 0xff4d6d };
const links = DATA.links.map(l => {
  const g = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3()]);
  const line = new THREE.Line(g, new THREE.LineBasicMaterial({ color: LINK_COL[l.source] ?? 0xffffff, transparent: true, opacity: 0.9 }));
  line.visible = false; scene.add(line);
  return { ...l, line };
});
const active = (l, t) => l.windows.some(([s, e]) => t >= s && t <= e);

// ---------------------------------------------------------------- time
let t = 0, playing = true, last = performance.now();
const spd = document.getElementById("spd"), spdv = document.getElementById("spdv"), scrub = document.getElementById("scrub");
const speed = () => Math.round(Math.pow(10, +spd.value));
spd.oninput = () => spdv.textContent = speed(); spd.oninput();
document.getElementById("play").onclick = e => { playing = !playing; e.target.textContent = playing ? "Pause" : "Play"; };
scrub.oninput = () => { t = +scrub.value * T_END; };
const epoch = new Date(DATA.epoch_utc + "Z");
const tmp = new THREE.Vector3(), tmp2 = new THREE.Vector3();

function posAt(s, t) {                     // linear interpolation between samples
  const f = Math.min(Math.max((t - T[0]) / (T[1] - T[0]), 0), N - 1);
  const i = Math.min(Math.floor(f), N - 2), w = f - i;
  return tmp.copy(s.pts[i]).lerp(s.pts[i + 1], w);
}

function frame(now) {
  const dt = (now - last) / 1000; last = now;
  if (playing) { t += dt * speed(); if (t > T_END) t = 0; scrub.value = t / T_END; }
  const ecef = document.getElementById("ecef").checked;
  const theta = DATA.gmst0 + DATA.omega_e * t;          // Earth rotation angle
  // inertial view: Earth rotates; Earth-fixed view: rotate the whole sky the other way
  earth.rotation.y = ecef ? 0 : theta;
  const frameRot = ecef ? -theta : 0;
  const showTrails = document.getElementById("trails").checked, showLbl = document.getElementById("labels").checked;
  const q = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 1, 0), frameRot);
  for (const s of sats) {
    s.mesh.position.copy(posAt(s, t)).applyQuaternion(q);
    s.trail.visible = showTrails && !ecef; s.trail.quaternion.copy(q);
    const v = tmp2.copy(s.mesh.position).project(camera);
    s.lbl.style.display = showLbl && v.z < 1 ? "block" : "none";
    s.lbl.style.left = ((v.x + 1) / 2 * innerWidth + 8) + "px";
    s.lbl.style.top = ((1 - v.y) / 2 * innerHeight - 6) + "px";
  }
  const showLinks = document.getElementById("showlinks").checked;
  const on = [];
  for (const l of links) {
    const vis = showLinks && active(l, t);
    l.line.visible = vis;
    if (vis) {
      const p = l.line.geometry.attributes.position;
      const A = sats[l.a].mesh.position, B = sats[l.b].mesh.position;
      p.setXYZ(0, A.x, A.y, A.z); p.setXYZ(1, B.x, B.y, B.z); p.needsUpdate = true;
      const d = A.distanceTo(B) / S;
      on.push(`${sats[l.a].name} ↔ ${sats[l.b].name}: ${d.toFixed(0)} km <span style="color:#7d8aa5">(${l.source})</span>`);
    }
  }
  document.getElementById("links").innerHTML = on.length ? "<b>Active links</b><br>" + on.join("<br>") : "<span style='color:#7d8aa5'>no active links</span>";
  const d = new Date(epoch.getTime() + t * 1000);
  document.getElementById("time").textContent = d.toISOString().replace("T", " ").slice(0, 19) + " UTC   (t = " + t.toFixed(0) + " s)";
  controls.update(); renderer.render(scene, camera);
  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);
addEventListener("resize", () => { camera.aspect = innerWidth / innerHeight; camera.updateProjectionMatrix(); renderer.setSize(innerWidth, innerHeight); });
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
