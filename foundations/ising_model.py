# ISING PHASE TRANSITION (Beyond the Index #1; Statistical
# Physics, SIM -- Metropolis Monte Carlo). A lattice of spins
# that each want to align with their neighbours. Heat fights the
# alignment; below a critical temperature order wins SPONTAN-
# EOUSLY -- no field is applied, the magnetization chooses a
# direction by itself. Onsager's 1944 exact solution fixes the
# 2D square-lattice critical temperature at Tc = 2/ln(1+sqrt 2)
# = 2.26919 (J = kB = 1). The simulation knows none of that.
# Gates:
#   G1 susceptibility peak lands within 3% of 2.26919 (for
#      L = 48 the peak is expected a touch ABOVE Tc from the
#      finite-size shift ~ L^-1 -- the gate has room for it)
#   G2 ordered phase: <|m|> > 0.85 at T = 1.6
#   G3 disordered phase: <|m|> < 0.15 at T = 3.2

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20261105)
L = 48
TC_EXACT = 2 / np.log(1 + np.sqrt(2))


def sweep(S, T, color_mask, nbr_sum_pre):
    # checkerboard Metropolis: all spins of one colour at once
    for color in (0, 1):
        nbr = (np.roll(S, 1, 0) + np.roll(S, -1, 0)
               + np.roll(S, 1, 1) + np.roll(S, -1, 1))
        dE = 2 * S * nbr
        accept = (dE <= 0) | (rng.random(S.shape)
                              < np.exp(-dE / T))
        flip = accept & (color_mask == color)
        S = np.where(flip, -S, S)
    return S


ii, jj = np.indices((L, L))
COLOR = (ii + jj) % 2

Ts = np.concatenate([np.arange(1.6, 2.0, 0.1),
                     np.arange(2.0, 2.14, 0.04),
                     np.arange(2.15, 2.42, 0.02),
                     np.arange(2.44, 2.62, 0.04),
                     np.arange(2.7, 3.25, 0.1)])
rows = []
S = rng.choice([-1, 1], size=(L, L))
for T in Ts:
    # critical slowing down: near Tc the lattice decorrelates
    # slowly, so the first version's 400+600 sweeps put the chi
    # peak 4.1% high and G1 failed. Near Tc, run 6x longer.
    eq, meas = (3000, 3000) if 2.1 <= T <= 2.45 else (400, 600)
    for _ in range(eq):
        S = sweep(S, T, COLOR, None)
    ms = []
    for _ in range(meas):
        S = sweep(S, T, COLOR, None)
        ms.append(abs(S.mean()))
    ms = np.array(ms)
    chi = L * L * (np.mean(ms ** 2) - np.mean(ms) ** 2) / T
    rows.append({"T": float(T), "m": float(ms.mean()),
                 "chi": float(chi)})
    print(f"T = {T:.2f}: <|m|> = {ms.mean():.3f}  chi = {chi:.1f}")

ch = np.array([r["chi"] for r in rows])
Tt = np.array([r["T"] for r in rows])
k = int(np.argmax(ch))
# parabolic refinement of the peak
if 0 < k < len(ch) - 1:
    x = Tt[k - 1:k + 2]
    yv = ch[k - 1:k + 2]
    a, b, _ = np.polyfit(x, yv, 2)
    tc_est = float(-b / (2 * a))
else:
    tc_est = float(Tt[k])
print(f"\nsusceptibility peak at T = {tc_est:.4f} "
      f"(Onsager exact {TC_EXACT:.5f}; error "
      f"{100 * abs(tc_est - TC_EXACT) / TC_EXACT:.2f}%)")
m_lo = rows[0]["m"]
m_hi = rows[-1]["m"]
g1 = abs(tc_est - TC_EXACT) / TC_EXACT < 0.03
g2 = m_lo > 0.85
g3 = m_hi < 0.15
print(f"GATES: G1 {g1} | G2 ordered {g2} ({m_lo:.3f}) | "
      f"G3 disordered {g3} ({m_hi:.3f})")

fig, ax = plt.subplots(figsize=(7.6, 4.2))
ax.plot(Tt, [r["m"] for r in rows], "o-", label="<|m|>")
ax.axvline(TC_EXACT, ls="--", color="gray",
           label=f"Onsager Tc = {TC_EXACT:.3f}")
ax.axvline(tc_est, ls=":", color="red",
           label=f"simulated peak {tc_est:.3f}")
ax.set_xlabel("temperature (J/kB)")
ax.set_ylabel("magnetization per spin")
ax.set_title("Ising model: spontaneous order below Tc")
ax.legend()
fig.tight_layout()
fig.savefig("ising_model.png", dpi=150)

out = {"Tc_exact": float(TC_EXACT), "Tc_estimate": tc_est,
       "m_at_1.6": m_lo, "m_at_3.2": m_hi, "rows": rows,
       "gates": {"G1_Tc": bool(g1), "G2_ordered": bool(g2),
                 "G3_disordered": bool(g3)}}
with open("ising_model_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values())
      else "SOME GATE FAILED")
print("Saved ising_model_results.json + PNG")
