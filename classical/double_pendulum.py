# DOUBLE PENDULUM CHAOS (Beyond the Index #2; Classical
# Mechanics / Chaos, SIM). Two identical double pendulums, one
# released a billionth of a radian from the other. In a chaotic
# regime their separation grows exponentially, |d(t)| ~ e^(Lt) --
# L is the largest Lyapunov exponent, and measuring it turns
# "chaotic" from an adjective into a number. A small-oscillation
# pair is run as a control: there the motion is quasi-periodic
# and L must come out ~ 0, proving the integrator is not
# manufacturing the chaos it reports.
# Gates:
#   G1 chaotic pair fully decorrelates (separation reaches O(1))
#   G2 Lyapunov estimate positive and order-unity: [0.3, 4.0] /s
#   G3 control pair: |L| < 0.05 /s

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G_, L1, L2, M1, M2 = 9.81, 1.0, 1.0, 1.0, 1.0


def deriv(s):
    t1, w1, t2, w2 = s
    d = t1 - t2
    den = 2 * M1 + M2 - M2 * np.cos(2 * d)
    a1 = (-G_ * (2 * M1 + M2) * np.sin(t1) - M2 * G_ * np.sin(t1 - 2 * t2)
          - 2 * np.sin(d) * M2 * (w2 ** 2 * L2 + w1 ** 2 * L1 * np.cos(d))) \
        / (L1 * den)
    a2 = (2 * np.sin(d) * (w1 ** 2 * L1 * (M1 + M2)
          + G_ * (M1 + M2) * np.cos(t1)
          + w2 ** 2 * L2 * M2 * np.cos(d))) / (L2 * den)
    return np.array([w1, a1, w2, a2])


def rk4(s, dt):
    k1 = deriv(s)
    k2 = deriv(s + dt / 2 * k1)
    k3 = deriv(s + dt / 2 * k2)
    k4 = deriv(s + dt * k3)
    return s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def lyapunov(s0, seconds=60.0, dt=0.005, renorm=1.0):
    s = np.array(s0, dtype=float)
    d0 = 1e-9
    pert = np.array([d0, 0, 0, 0])
    sp = s + pert
    logs, steps_per = [], int(renorm / dt)
    sep_final = 0.0
    for n in range(int(seconds / renorm)):
        for _ in range(steps_per):
            s = rk4(s, dt)
            sp = rk4(sp, dt)
        dvec = sp - s
        dist = np.linalg.norm(dvec)
        logs.append(np.log(dist / d0))
        sep_final = dist
        sp = s + dvec / dist * d0
    lam = float(np.sum(logs) / seconds)
    return lam, float(sep_final), np.cumsum(logs)


lam_c, sep_c, cum_c = lyapunov([2.2, 0.0, 2.8, 0.0])
lam_k, sep_k, _ = lyapunov([0.10, 0.0, 0.15, 0.0])
print(f"chaotic pair:  L = {lam_c:.3f} /s, raw separation after "
      f"60 s of renormalized tracking: renorm sum "
      f"{cum_c[-1]:.1f} (e^{cum_c[-1]:.0f} total growth)")
print(f"control pair:  L = {lam_k:+.4f} /s (must be ~0)")
# decorrelation check without renormalization, short window
s = np.array([2.2, 0.0, 2.8, 0.0])
sp = s + np.array([1e-9, 0, 0, 0])
tdec = None
for n in range(int(30 / 0.005)):
    s = rk4(s, 0.005)
    sp = rk4(sp, 0.005)
    if np.linalg.norm(sp - s) > 0.5:
        tdec = (n + 1) * 0.005
        break
print(f"chaotic pair reaches O(1) separation after {tdec} s "
      f"(predictability horizon ~ ln(1e9)/L = "
      f"{np.log(1e9) / lam_c:.1f} s)")
g1 = tdec is not None
g2 = 0.3 <= lam_c <= 4.0
g3 = abs(lam_k) < 0.05
print(f"GATES: G1 {g1} | G2 {g2} | G3 {g3}")

fig, ax = plt.subplots(figsize=(7.6, 4.0))
ax.plot(np.arange(1, 61) * 1.0, cum_c / np.log(10))
ax.set_xlabel("time (s)")
ax.set_ylabel("separation growth (orders of magnitude)")
ax.set_title(f"Double pendulum: a 1e-9 nudge grows e-fold by e-fold"
             f" (L = {lam_c:.2f} /s)")
fig.tight_layout()
fig.savefig("double_pendulum.png", dpi=150)

out = {"lyapunov_chaotic_per_s": lam_c,
       "lyapunov_control_per_s": lam_k,
       "decorrelation_time_s": tdec,
       "predictability_horizon_s": float(np.log(1e9) / lam_c),
       "gates": {"G1_decorrelates": bool(g1), "G2_lambda": bool(g2),
                 "G3_control": bool(g3)}}
with open("double_pendulum_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values())
      else "SOME GATE FAILED")
print("Saved double_pendulum_results.json + PNG")
