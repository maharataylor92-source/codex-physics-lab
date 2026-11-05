# MERCURY'S PERIHELION (Beyond the Index #5; General
# Relativity, SIM). Newton's gravity gives a closed ellipse.
# Einstein's adds one term to the orbit equation,
#   d2u/dth2 + u = GM/h2 + (3GM/c2) u2 ,   u = 1/r,
# and the ellipse stops closing: each perihelion lands a little
# past the last. For Mercury the accumulated drift is the
# famous 43 arcseconds per century that Newtonian perturbation
# theory could not explain. This integrates the equation
# directly (RK4 in theta) over 200 orbits and measures the
# drift; a control run without the GR term must NOT drift.
# Gates:
#   G1 GR run: total precession within 5% of 42.98 "/century
#   G2 control (Newton): |drift| < 2% of the GR value

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GM = 1.32712440018e20          # m^3/s^2, Sun
C = 299792458.0
A_SEMI = 5.7909e10             # Mercury semi-major axis, m
ECC = 0.20563
P_ORB_YR = 0.240846            # Mercury year
ORBITS = 200
TARGET = 42.98                 # arcsec per century


def run(gr=True, dth=0.00025):
    p = A_SEMI * (1 - ECC ** 2)          # semi-latus rectum
    h2 = GM * p
    u = GM / h2 * (1 + ECC)              # perihelion: u max
    w = 0.0                              # du/dth at perihelion
    coef = 3 * GM / C ** 2 if gr else 0.0

    def acc(uu):
        return GM / h2 + coef * uu * uu - uu

    peri_theta = [0.0]
    th = 0.0
    w_prev = w
    u_hist = []
    total = int(ORBITS * 2 * np.pi / dth)
    for _ in range(total):
        # RK4 on (u, w)
        def f(uu, ww):
            return ww, acc(uu)
        k1u, k1w = f(u, w)
        k2u, k2w = f(u + dth / 2 * k1u, w + dth / 2 * k1w)
        k3u, k3w = f(u + dth / 2 * k2u, w + dth / 2 * k2w)
        k4u, k4w = f(u + dth * k3u, w + dth * k3w)
        u += dth / 6 * (k1u + 2 * k2u + 2 * k3u + k4u)
        w += dth / 6 * (k1w + 2 * k2w + 2 * k3w + k4w)
        th += dth
        u_hist.append(u)
        if w_prev > 0 and w <= 0 and th > 1.0:
            # sub-step perihelion: parabolic peak of u through
            # the last three samples (the first version took the
            # raw grid point; with the coarse dth = 0.002 grid
            # RK4 also accumulated a phase drift -- the Newton
            # control drifted +5.6 "/century, which no closed
            # ellipse may do. dth is now 8x finer (phase error
            # down ~4000x) and peaks are interpolated.)
            y0, y1, y2 = u_hist[-3], u_hist[-2], u_hist[-1]
            denom = (y0 - 2 * y1 + y2)
            off = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
            peri_theta.append(th - dth + off * dth)
        w_prev = w
    peri_theta = np.array(peri_theta)
    k = np.arange(len(peri_theta))
    drift_per_orbit = float(np.polyfit(k, peri_theta - 2 * np.pi * k,
                                       1)[0])
    per_century = drift_per_orbit * (100 / P_ORB_YR) \
        * 206264.806  # rad -> arcsec
    return per_century, len(peri_theta)


gr_val, n1 = run(gr=True)
nw_val, n2 = run(gr=False)
print(f"GR run:      {gr_val:.2f} arcsec/century "
      f"({n1} perihelia; target {TARGET})")
print(f"Newton run:  {nw_val:+.3f} arcsec/century (must be ~0)")
g1 = abs(gr_val - TARGET) / TARGET < 0.05
g2 = abs(nw_val) < 0.02 * gr_val
print(f"GATES: G1 {g1} | G2 {g2}")

fig, ax = plt.subplots(figsize=(7.4, 4.0))
orbits_c = 100 / P_ORB_YR
ax.bar(["Newton\n(control)", "Newton + GR term"],
       [nw_val, gr_val])
ax.axhline(TARGET, ls="--", color="red",
           label="observed anomaly 42.98\"/century")
ax.set_ylabel("perihelion advance (arcsec / century)")
ax.set_title("Mercury: the 43 arcseconds Newton could not explain")
ax.legend()
fig.tight_layout()
fig.savefig("mercury_perihelion.png", dpi=150)

out = {"gr_arcsec_per_century": gr_val,
       "newton_arcsec_per_century": nw_val, "target": TARGET,
       "gates": {"G1_gr_matches": bool(g1), "G2_control": bool(g2)}}
with open("mercury_perihelion_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values())
      else "SOME GATE FAILED")
print("Saved mercury_perihelion_results.json + PNG")
