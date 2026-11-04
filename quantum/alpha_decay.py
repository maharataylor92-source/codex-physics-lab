# TUNNELING DEEPENED, part 2 -- ALPHA DECAY: the Gamow theory
# against REAL measured half-lives (extends tunneling.py, XIX).
#
# The same tunneling formula, pointed at nature. An alpha particle
# rattles inside the nucleus ~10^21 times a second; each hit on the
# Coulomb barrier transmits with WKB probability T = exp(-2*gamma),
#
#   gamma = (sqrt(2 mu)/hbar) * integral_R^b sqrt(V(r) - E) dr,
#   V(r) = Z1 Z2 e^2 / (4 pi eps0 r),   b = Z1 Z2 e^2/(4 pi eps0 E)
#
# with the closed form integral_R^b sqrt(b/r - 1) dr =
# b * [ arccos(sqrt(R/b)) - sqrt((R/b)(1 - R/b)) ].
# Half-life t1/2 = ln2 / (f * T), assault frequency f = v/(2R),
# v = sqrt(2E/mu), mu = reduced mass, E = Q (total release).
# NO fitted parameters: R = r0 (A_d^(1/3) + 4^(1/3)), r0 = 1.2 fm.
#
# Data: even-even emitters, Q and half-lives from standard nuclear
# data tables (textbook/NNDC values; spot-verified 2026-11-04).
#
# The famous fact to reproduce: a factor ~2 in Q spans ~24 orders
# of magnitude in half-life (Geiger-Nuttall law). Gates:
#   G1 predicted log10(t1/2) correlates with measured at r >= 0.98
#   G2 predicted span within 6 decades of the measured 24.2
#   G3 every prediction within 4 decades of measurement (crude
#      theory tolerance -- residuals are reported, not hidden)

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# constants (SI)
hbar = 1.054571817e-34
e2 = 1.4399764e-9 * 1.602176634e-19   # e^2/(4 pi eps0) in J*m (=1.44 MeV fm)
u = 1.66053906660e-27                  # kg
MeV = 1.602176634e-13                  # J
yr = 365.25 * 86400.0

# name, A_parent, Z_daughter, Q (MeV), measured t1/2 (s)
DATA = [
    ("232Th", 232, 88, 4.083, 1.40e10 * yr),
    ("238U", 238, 90, 4.270, 4.468e9 * yr),
    ("230Th", 230, 88, 4.770, 7.54e4 * yr),
    ("226Ra", 226, 86, 4.871, 1600.0 * yr),
    ("210Po", 210, 82, 5.407, 138.376 * 86400.0),
    ("238Pu", 238, 92, 5.593, 87.7 * yr),
    ("218Po", 218, 82, 6.115, 3.10 * 60.0),
    ("220Rn", 220, 84, 6.405, 55.6),
    ("214Po", 214, 82, 7.833, 1.64e-4),
    ("212Po", 212, 82, 8.954, 2.99e-7),
]

r0 = 1.2e-15
rows = []
for name, A, Zd, Q, t_meas in DATA:
    Ad = A - 4
    R = r0 * (Ad ** (1 / 3) + 4 ** (1 / 3))
    E = Q * MeV
    mu = (4.0 * Ad / (4.0 + Ad)) * u
    b = 2 * Zd * e2 / E                       # Z_alpha = 2
    x = R / b
    integ = b * (np.arccos(np.sqrt(x)) - np.sqrt(x * (1 - x)))
    gamma = np.sqrt(2 * mu) / hbar * np.sqrt(E) * integ
    T = np.exp(-2 * gamma)
    v = np.sqrt(2 * E / mu)
    f = v / (2 * R)
    t_pred = np.log(2) / (f * T)
    rows.append({"isotope": name, "Q_MeV": Q,
                 "gamma": float(gamma), "T_per_hit": float(T),
                 "assault_freq_hz": float(f),
                 "t_pred_s": float(t_pred), "t_meas_s": float(t_meas),
                 "log10_pred": float(np.log10(t_pred)),
                 "log10_meas": float(np.log10(t_meas)),
                 "resid_decades": float(np.log10(t_pred) - np.log10(t_meas))})
    r = rows[-1]
    print(f"{name:>6}: Q={Q:.3f}  gamma={r['gamma']:.1f}  "
          f"log10 t_pred={r['log10_pred']:7.2f}  "
          f"log10 t_meas={r['log10_meas']:7.2f}  "
          f"resid={r['resid_decades']:+.2f} decades")

lp = np.array([r["log10_pred"] for r in rows])
lm = np.array([r["log10_meas"] for r in rows])
corr = float(np.corrcoef(lp, lm)[0, 1])
span_p, span_m = float(lp.max() - lp.min()), float(lm.max() - lm.min())
rms = float(np.sqrt(np.mean((lp - lm) ** 2)))
# Geiger-Nuttall fit of the MEASURED data: log10 t = a + c / sqrt(Q)
Qs = np.array([r["Q_MeV"] for r in rows])
c_m, a_m = np.polyfit(1 / np.sqrt(Qs), lm, 1)
c_p, a_p = np.polyfit(1 / np.sqrt(Qs), lp, 1)
print(f"\ncorrelation r = {corr:.4f} | span pred {span_p:.1f} vs meas "
      f"{span_m:.1f} decades | RMS residual {rms:.2f} decades")
print(f"Geiger-Nuttall slope: predicted {c_p:.0f} vs measured {c_m:.0f} "
      f"(per sqrt(MeV))")
gates = {"G1_correlation>=0.98": corr >= 0.98,
         "G2_span_within_6_decades": abs(span_p - span_m) <= 6,
         "G3_all_within_4_decades": bool(np.all(np.abs(lp - lm) <= 4))}
print(gates)

fig, ax = plt.subplots(figsize=(7.2, 4.0))
ax.plot(1 / np.sqrt(Qs), lm, "o", label="measured")
ax.plot(1 / np.sqrt(Qs), lp, "s", label="Gamow prediction (no fit)")
for r, xx in zip(rows, 1 / np.sqrt(Qs)):
    ax.annotate(r["isotope"], (xx, r["log10_meas"]), fontsize=7,
                xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("1/sqrt(Q)  (MeV^-1/2)")
ax.set_ylabel("log10 half-life (s)")
ax.set_title("Geiger-Nuttall: Gamow tunneling vs measured alpha half-lives")
ax.legend()
fig.tight_layout()
fig.savefig("alpha_decay.png", dpi=150)

out = {"model": "Gamow WKB, r0=1.2 fm, no fitted parameters",
       "rows": rows, "correlation": corr,
       "span_pred_decades": span_p, "span_meas_decades": span_m,
       "rms_residual_decades": rms,
       "GN_slope_pred": float(c_p), "GN_slope_meas": float(c_m),
       "gates": {k: bool(v) for k, v in gates.items()}}
with open("alpha_decay_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(gates.values()) else "SOME GATE FAILED")
print("Saved alpha_decay_results.json + PNG")
