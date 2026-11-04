# Orbit simulator -- Codex Vol I, Celestial Mechanics: v = sqrt(GM/r).
#
# Newtonian two-body integration (RK4). Dimensionless units GM = 1
# unless stated. Experiments:
#   1. Circular orbit at v_circ = sqrt(GM/r): radius holds, period
#      matches Kepler T = 2*pi*sqrt(r^3/GM).
#   2. Radius sweep: v_circ and T at r = 1, 2, 4, 8 -- Kepler scaling
#      T ~ r^(3/2) checked against integration.
#   3. Elliptical orbit at 0.8 * v_circ: measured period vs Kepler
#      using the semi-major axis from vis-viva / energy.
#   4. Escape: v_esc = sqrt(2) * v_circ; a 1.05 * v_esc launch has
#      positive total energy and never returns.
#   5. Real scale: Earth around the Sun; ISS around the Earth.
#
# Gates: circular radius drift < 1e-3 over 5 periods; energy drift
# < 1e-8 relative; Kepler ratio T(8)/T(1) = 8^(3/2) within 1%;
# ellipse period within 2% of Kepler-from-a.

import json

import numpy as np

GM = 1.0


def accel(r):
    return -GM * r / np.linalg.norm(r) ** 3


def integrate(r0, v0, t_final, dt):
    n = int(t_final / dt)
    rs, vs = np.empty((n + 1, 2)), np.empty((n + 1, 2))
    r, v = np.array(r0, float), np.array(v0, float)
    rs[0], vs[0] = r, v
    for i in range(n):
        def f(r_, v_):
            return v_, accel(r_)
        k1r, k1v = f(r, v)
        k2r, k2v = f(r + k1r * dt / 2, v + k1v * dt / 2)
        k3r, k3v = f(r + k2r * dt / 2, v + k2v * dt / 2)
        k4r, k4v = f(r + k3r * dt, v + k3v * dt)
        r = r + dt / 6 * (k1r + 2 * k2r + 2 * k3r + k4r)
        v = v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        rs[i + 1], vs[i + 1] = r, v
    return rs, vs


def energy(r, v):
    return 0.5 * np.dot(v, v) - GM / np.linalg.norm(r)


def period_kepler(a):
    return 2 * np.pi * np.sqrt(a**3 / GM)


print("=== Gates ===")
# Circular orbit, r = 1
r0, v0 = [1.0, 0.0], [0.0, 1.0]
T1 = period_kepler(1.0)
rs, vs = integrate(r0, v0, 5 * T1, T1 / 4000)
rad = np.linalg.norm(rs, axis=1)
drift = float(np.max(np.abs(rad - 1.0)))
E = np.array([energy(r, v) for r, v in zip(rs, vs)])
edrift = float(np.max(np.abs(E - E[0])) / abs(E[0]))
g1 = drift < 1e-3 and edrift < 1e-8
print(f"Gate 1 (circular r=1, 5 periods): radius drift {drift:.2e}, "
      f"energy drift {edrift:.2e} -> {'PASS' if g1 else 'FAIL'}")

# Radius sweep with measured periods (time one full revolution by
# angle accumulation)
def measured_period(r):
    vc = np.sqrt(GM / r)
    T = period_kepler(r)
    dt = T / 6000
    rs_, _ = integrate([r, 0.0], [0.0, vc], T * 1.05, dt)
    ang = np.unwrap(np.arctan2(rs_[:, 1], rs_[:, 0]))
    idx = np.argmax(ang >= 2 * np.pi)
    return float(idx * dt), vc, T


sweep = {}
for r in (1.0, 2.0, 4.0, 8.0):
    Tm, vc, Tk = measured_period(r)
    sweep[str(r)] = {"v_circ": vc, "T_kepler": Tk, "T_measured": Tm}
    print(f"  r={r}: v_circ = {vc:.4f}, T_kepler = {Tk:.3f}, "
          f"T_measured = {Tm:.3f}")
ratio = sweep["8.0"]["T_measured"] / sweep["1.0"]["T_measured"]
g2 = abs(ratio / 8**1.5 - 1) < 0.01
print(f"Gate 2 (Kepler T(8)/T(1) = {ratio:.2f} vs 8^1.5 = {8**1.5:.2f}) "
      f"-> {'PASS' if g2 else 'FAIL'}")

# Elliptical orbit: v = 0.8 v_circ at r = 1 (this point is apoapsis)
ve = 0.8
Ee = energy([1.0, 0.0], [0.0, ve])
a_ell = -GM / (2 * Ee)
Tk_ell = period_kepler(a_ell)
rs_e, _ = integrate([1.0, 0.0], [0.0, ve], Tk_ell * 1.05, Tk_ell / 8000)
ang_e = np.unwrap(np.arctan2(rs_e[:, 1], rs_e[:, 0]))
idx_e = int(np.argmax(ang_e >= 2 * np.pi))
Tm_ell = idx_e * (Tk_ell / 8000)
rad_e = np.linalg.norm(rs_e, axis=1)
g3 = abs(Tm_ell / Tk_ell - 1) < 0.02
print(f"Gate 3 (ellipse 0.8 v_circ): a = {a_ell:.4f}, r in "
      f"[{rad_e.min():.3f}, {rad_e.max():.3f}], T_kepler = {Tk_ell:.3f}, "
      f"T_measured = {Tm_ell:.3f} -> {'PASS' if g3 else 'FAIL'}")

# Escape
vesc = np.sqrt(2 * GM / 1.0)
Eesc = energy([1.0, 0.0], [0.0, 1.05 * vesc])
rs_x, _ = integrate([1.0, 0.0], [0.0, 1.05 * vesc], 30.0, 0.002)
rad_x = np.linalg.norm(rs_x, axis=1)
g4 = Eesc > 0 and rad_x[-1] > 10 and rad_x[-1] > rad_x[len(rad_x) // 2]
print(f"Gate 4 (1.05 v_esc): E = {Eesc:+.4f} > 0, r_final = "
      f"{rad_x[-1]:.1f} and still growing -> {'PASS' if g4 else 'FAIL'}")

# Real scale
GM_sun, AU = 1.32712440018e20, 1.495978707e11
v_earth = np.sqrt(GM_sun / AU)
T_earth = 2 * np.pi * np.sqrt(AU**3 / GM_sun)
GM_earth, r_iss = 3.986004418e14, 6.771e6
v_iss = np.sqrt(GM_earth / r_iss)
T_iss = 2 * np.pi * np.sqrt(r_iss**3 / GM_earth)
print("\n=== Real scale ===")
print(f"  Earth around Sun: v = {v_earth/1000:.2f} km/s, "
      f"T = {T_earth/86400:.2f} days")
print(f"  ISS around Earth: v = {v_iss/1000:.2f} km/s, "
      f"T = {T_iss/60:.1f} min")

out = {"model": "Newtonian two-body, RK4, GM=1 dimensionless + SI cases",
       "gates": {"circular": bool(g1), "kepler_scaling": bool(g2),
                 "ellipse": bool(g3), "escape": bool(g4)},
       "sweep": sweep,
       "ellipse": {"a": float(a_ell), "r_min": float(rad_e.min()),
                   "r_max": float(rad_e.max()), "T_kepler": float(Tk_ell),
                   "T_measured": float(Tm_ell)},
       "real_scale": {"earth_sun": {"v_km_s": v_earth / 1000,
                                    "T_days": T_earth / 86400},
                      "iss": {"v_km_s": v_iss / 1000, "T_min": T_iss / 60}}}
with open("orbit_results.json", "w") as f:
    json.dump(out, f, indent=2)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
axes[0].plot(rs[::20, 0], rs[::20, 1], label="circular r=1")
axes[0].plot(rs_e[::20, 0], rs_e[::20, 1], label="ellipse 0.8 v_circ")
axes[0].plot([0], [0], "yo", ms=12)
axes[0].set_aspect("equal"); axes[0].legend(); axes[0].grid(alpha=0.3)
axes[0].set_title("Orbits (GM = 1)")
rr = np.array([1.0, 2.0, 4.0, 8.0])
axes[1].loglog(rr, [sweep[str(r)]["v_circ"] for r in rr], "o-",
               label="v_circ ~ r^-1/2")
axes[1].loglog(rr, [sweep[str(r)]["T_measured"] for r in rr], "s-",
               label="T ~ r^3/2 (measured)")
axes[1].legend(); axes[1].grid(alpha=0.3, which="both")
axes[1].set_xlabel("orbital radius r"); axes[1].set_title("Kepler scaling")
fig.tight_layout()
fig.savefig("orbit_simulator.png", dpi=150)
print("\nSaved orbit_results.json and orbit_simulator.png")
