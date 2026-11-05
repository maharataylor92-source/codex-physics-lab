# PERCOLATION THRESHOLD (Beyond the Index #7; Statistical
# Physics, SIM). Fill each cell of a grid independently with
# probability p. At low p the filled cells form small islands;
# at high p they merge into a continent. Somewhere between,
# a single cluster first SPANS the grid from top to bottom --
# a geometric phase transition. For site percolation on the
# square lattice the threshold is p_c = 0.592746... (known
# numerically to high precision; no closed form exists).
# Estimated here as the p where the spanning probability
# crosses 1/2, averaged over random realizations (clusters by
# connected-component labelling, 4-connectivity).
# Gates:
#   G1 p_c estimate within 0.01 of 0.5927
#   G2 P_span(0.50) < 0.05 and G3 P_span(0.68) > 0.95

import json
import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20261105)
L = 128
TRIALS = 200
PC_KNOWN = 0.5927


def spans(p):
    grid = rng.random((L, L)) < p
    lab, _ = ndimage.label(grid)
    top = set(np.unique(lab[0])) - {0}
    bot = set(np.unique(lab[-1])) - {0}
    return len(top & bot) > 0


ps = np.arange(0.50, 0.685, 0.01)
P = []
for p in ps:
    frac = np.mean([spans(p) for _ in range(TRIALS)])
    P.append(float(frac))
    print(f"p = {p:.2f}: P_span = {frac:.3f}")
P = np.array(P)
k = int(np.argmax(P >= 0.5))
pc_est = float(np.interp(0.5, P[k - 1:k + 1], ps[k - 1:k + 1]))
print(f"\nestimated p_c = {pc_est:.4f} (known {PC_KNOWN}; "
      f"error {abs(pc_est - PC_KNOWN):.4f})")
g1 = abs(pc_est - PC_KNOWN) < 0.01
g2 = P[0] < 0.05
g3 = P[-1] > 0.95
print(f"GATES: G1 {g1} | G2 {g2} | G3 {g3}")

fig, ax = plt.subplots(figsize=(7.4, 4.2))
ax.plot(ps, P, "o-")
ax.axvline(PC_KNOWN, ls="--", color="red",
           label=f"known p_c = {PC_KNOWN}")
ax.axhline(0.5, ls=":", color="gray")
ax.set_xlabel("occupation probability p")
ax.set_ylabel("spanning probability")
ax.set_title(f"Percolation: the spanning transition "
             f"(estimated p_c = {pc_est:.4f})")
ax.legend()
fig.tight_layout()
fig.savefig("percolation.png", dpi=150)

out = {"pc_estimate": pc_est, "pc_known": PC_KNOWN,
       "P_at_0.50": float(P[0]), "P_at_0.68": float(P[-1]),
       "ps": ps.tolist(), "P": P.tolist(),
       "gates": {"G1_pc": bool(g1), "G2_below": bool(g2),
                 "G3_above": bool(g3)}}
with open("percolation_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values())
      else "SOME GATE FAILED")
print("Saved percolation_results.json + PNG")
