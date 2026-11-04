# AMDAHL'S LAW (Codex index #17; Computing, NUM).
# If a fraction p of a workload parallelizes and the rest is
# serial, s processors give speedup S = 1 / ((1 - p) + p/s),
# capped at 1/(1 - p) no matter how many processors arrive.
# The serial fraction is the wall -- the computing cousin of the
# lab's other ceilings (T2 for circuits, resonance width for
# tunneling): every system in this repo has one.
# Gates (identities, stated as such):
#   G1 S(1) = 1 exactly for every p
#   G2 S(inf) = 1/(1 - p) (evaluated at s = 10^8, relative tol)
#   G3 the textbook anchor: p = 0.95 caps at 20x; at s = 20 the
#      speedup is 10.26x (half the ceiling, double the machine
#      count away from 19x -- diminishing returns made numeric)

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PS = [0.5, 0.75, 0.9, 0.95, 0.99]
SS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 1024]


def S(p, s):
    return 1.0 / ((1 - p) + p / s)


table = {str(p): {str(s): float(S(p, s)) for s in SS} for p in PS}
for p in PS:
    row = " ".join(f"{S(p, s):7.2f}" for s in SS)
    print(f"p={p:<5} ceiling {1 / (1 - p):6.1f}x | S(s): {row}")

g1 = all(abs(S(p, 1) - 1.0) < 1e-12 for p in PS)
g2 = all(abs(S(p, 10 ** 8) - 1 / (1 - p)) / (1 / (1 - p)) < 1e-4
         for p in PS)  # relative: approach is slow, err ~ ceiling^2 p/s
g3 = abs(S(0.95, 20) - 10.2564) < 0.001
print(f"G1 S(1)=1: {g1} | G2 ceiling identity: {g2} | "
      f"G3 p=.95,s=20 -> {S(0.95, 20):.4f}: {g3}")

# the wall in people terms: to reach 90% of the ceiling for
# p = 0.95 you need s ~ 171 processors; for p = 0.99, s ~ 891.
def s_for_fraction(p, frac):
    target = frac / (1 - p)
    return p / (1 / target - (1 - p))


walls = {str(p): float(s_for_fraction(p, 0.9)) for p in PS}
for p in PS:
    print(f"p={p}: 90% of ceiling needs s = {s_for_fraction(p, 0.9):.0f} "
          f"processors")

fig, ax = plt.subplots(figsize=(7.4, 4.0))
ss = np.logspace(0, 3.5, 100)
for p in PS:
    ax.plot(ss, [S(p, s) for s in ss], label=f"p = {p}")
    ax.axhline(1 / (1 - p), color="gray", ls=":", lw=0.7)
ax.set_xscale("log")
ax.set_xlabel("processors s")
ax.set_ylabel("speedup S")
ax.set_title("Amdahl's law: the serial fraction is the wall")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("amdahl.png", dpi=150)

out = {"table": table, "ceilings": {str(p): 1 / (1 - p) for p in PS},
       "s_for_90pct_of_ceiling": walls,
       "gates": {"G1_S1": bool(g1), "G2_ceiling": bool(g2),
                 "G3_anchor": bool(g3)}}
with open("amdahl_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values()) else "SOME GATE FAILED")
print("Saved amdahl_results.json + PNG")
