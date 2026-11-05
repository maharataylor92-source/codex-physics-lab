# FEIGENBAUM'S CONSTANT (Beyond the Index #6; Chaos Theory,
# NUM). The logistic map x -> r x (1 - x) is one line of code,
# but as r grows its long-term behaviour doubles its period
# again and again -- 1, 2, 4, 8, 16... -- and the doublings
# crowd together geometrically. The ratio of successive gaps
# tends to delta = 4.669201609..., the SAME constant for every
# unimodal map (and measured in real fluids and circuits).
#
# Method note (two earlier attempts kept honestly): locating
# the bifurcation points r_n by period classification failed
# twice -- near a doubling, the new cycle attracts so weakly
# that finite transients misclassify the period at the 1e-4
# level, and the tiny late gaps amplify that into bad ratios
# (runs gave 4.599-4.691 for the third ratio). The robust
# route is the SUPERSTABLE points s_n, where the critical
# point x = 1/2 itself has period 2^n: s_n solves
# f^(2^n)(1/2) = 1/2, a finite polynomial equation in r with
# NO transients and NO classification -- bisection is exact.
# The s_n gaps shrink with the same universal ratio delta.
# Gates:
#   G1 s_0 = 2 and s_1 = 1 + sqrt(5) = 3.236068 (exact values)
#   G2 last ratio within 1% of 4.669201609
#   G3 net convergence: |last ratio - delta| < |first - delta|

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DELTA = 4.669201609


def F(r, n):
    """f^(2^n)(1/2) - 1/2, exact finite iteration."""
    x = 0.5
    for _ in range(2 ** n):
        x = r * x * (1 - x)
    return x - 0.5


def superstable(n, lo, hi):
    flo = F(lo, n)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (F(mid, n) > 0) == (flo > 0):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# brackets need care: f^(2^n)(1/2) = 1/2 is ALSO satisfied at
# every lower superstable point s_k (k < n), since a period
# dividing 2^n qualifies. The first version's fixed brackets
# made n=4 re-find s_3 and n=5 re-find s_4. Sequential brackets
# (s_{n-1}, r_inf) exclude all lower roots.
s = [superstable(0, 1.5, 2.5), superstable(1, 3.0, 3.4)]
for n in range(2, 6):
    lo = s[-1] + 1e-9
    assert F(lo, n) * F(3.5699, n) < 0, f"no bracket for n={n}"
    s.append(superstable(n, lo, 3.5699))
print("superstable points:", [f"{v:.6f}" for v in s])
deltas = [(s[i] - s[i - 1]) / (s[i + 1] - s[i])
          for i in range(1, 5)]
print("ratios:", [f"{d:.4f}" for d in deltas],
      f"-> Feigenbaum {DELTA}")
g1 = abs(s[0] - 2.0) < 1e-9 and abs(s[1] - (1 + np.sqrt(5))) < 1e-4
g2 = abs(deltas[-1] - DELTA) / DELTA < 0.01
g3 = abs(deltas[-1] - DELTA) < abs(deltas[0] - DELTA)
print(f"GATES: G1 exact anchors {g1} | G2 delta {g2} | "
      f"G3 net convergence {g3}")

fig, ax = plt.subplots(figsize=(7.6, 4.2))
rr = np.linspace(2.9, 3.57, 900)
for r in rr:
    x = 0.5
    for _ in range(600):
        x = r * x * (1 - x)
    xs = []
    for _ in range(24):
        x = r * x * (1 - x)
        xs.append(x)
    ax.plot([r] * len(xs), xs, ",", color="black", alpha=0.5)
ax.set_xlabel("r")
ax.set_ylabel("long-term x")
ax.set_title("Logistic map bifurcation cascade -- the gaps "
             "shrink by Feigenbaum's ratio")
fig.tight_layout()
fig.savefig("feigenbaum.png", dpi=150)

out = {"superstable_points": s, "ratios": deltas,
       "delta_exact": DELTA,
       "gates": {"G1_anchors": bool(g1), "G2_delta": bool(g2),
                 "G3_convergence": bool(g3)}}
with open("feigenbaum_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values())
      else "SOME GATE FAILED")
print("Saved feigenbaum_results.json + PNG")
