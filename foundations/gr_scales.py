# GR SCALE NUMBERS (Codex index #1; General Relativity, NUM).
# G = (8 pi G / c^4) T made into arithmetic: Schwarzschild radii
# for real masses, and gravitational time dilation where it is
# measured every day -- GPS.
# Gates (known values, stated as such):
#   G1 Sun r_s = 2.953 km and Earth r_s = 8.87 mm within 1%
#   G2 GPS gravitational term: +45.7 us/day within 3% (the famous
#      figure; with the SR term -7.2 the net is +38.5)
#   G3 neutron-star surface clock rate sqrt(1 - r_s/R) ~ 0.77

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = 6.67430e-11
C = 299792458.0
MSUN = 1.98847e30
MEARTH = 5.9722e24


def rs(m):
    return 2 * G * m / C ** 2


OBJECTS = [
    ("1 kg mass", 1.0),
    ("human (70 kg)", 70.0),
    ("Earth", MEARTH),
    ("Sun", MSUN),
    ("Sagittarius A* (4.3e6 Msun)", 4.3e6 * MSUN),
    ("M87* (6.5e9 Msun)", 6.5e9 * MSUN),
]
rows = []
for name, m in OBJECTS:
    r = rs(m)
    rows.append({"object": name, "mass_kg": m, "rs_m": float(r)})
    print(f"{name:<28}: r_s = {r:.4e} m")

g1 = abs(rs(MSUN) - 2953) / 2953 < 0.01 and abs(rs(MEARTH) - 8.87e-3) / 8.87e-3 < 0.01

# time dilation: clock rate at radius r around Earth
def rate(r, m=MEARTH):
    return float(np.sqrt(1 - rs(m) / r))


R_E = 6.371e6
R_GPS = 26.560e6
# the higher clock runs FASTER: rate increases with r
gain_us_day = (rate(R_GPS) / rate(R_E) - 1) * 86400e6
print(f"\nEarth surface clock rate factor: {rate(R_E):.12f}")
print(f"GPS orbit clock rate factor:     {rate(R_GPS):.12f}")
print(f"GPS gravitational gain: +{gain_us_day:.1f} us/day "
      f"(canonical +45.7; SR term -7.2 -> net +38.5)")
g2 = abs(gain_us_day - 45.7) / 45.7 < 0.03

ns_rate = float(np.sqrt(1 - rs(1.4 * MSUN) / 1.0e4))
print(f"Neutron star (1.4 Msun, R=10 km) surface clock rate: "
      f"{ns_rate:.3f} -- a clock there ticks {100 * (1 - ns_rate):.0f}% slow")
g3 = abs(ns_rate - 0.766) < 0.01
iss_rate = rate(6.371e6 + 420e3)
print(f"ISS clock vs ground: {(iss_rate / rate(R_E) - 1) * 86400e6:+.2f} "
      f"us/day gravitational (before its SR loss)")
print(f"GATES: G1 {g1} G2 {g2} G3 {g3}")

fig, ax = plt.subplots(figsize=(8.2, 3.6))
ax.barh([r["object"] for r in rows],
        [np.log10(r["rs_m"]) for r in rows])
ax.set_xlabel("log10 Schwarzschild radius (m)")
ax.set_title("Schwarzschild radii: from 1 kg to M87*")
fig.tight_layout()
fig.savefig("gr_scales.png", dpi=150)

out = {"rows": rows, "gps_gain_us_day": float(gain_us_day),
       "ns_surface_rate": ns_rate,
       "gates": {"G1_known_radii": bool(g1), "G2_gps": bool(g2),
                 "G3_neutron_star": bool(g3)}}
with open("gr_scales_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values()) else "SOME GATE FAILED")
print("Saved gr_scales_results.json + PNG")
