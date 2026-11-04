# Radical-pair field sweep -- extension of radical_pair.py (Vol II §XVI).
#
# The original model fixed Earth's field at 50 uT and swept the angle.
# This sweeps BOTH: field strength across Earth's real range
# (25-65 uT, per Vol II §XVI) x field angle, using the original model's
# own functions (source loaded up to its verification gates, so the
# physics is byte-identical: omega scales linearly with B).
#
# Sanity check built in: at B = 50 uT the sweep must reproduce the
# original run (Phi_S(0 deg) ~ 0.5075, Phi_S(90 deg) ~ 0.6719).

import json

import numpy as np

src = open("radical_pair.py").read()
head = src.split("# ---------- verification gates ----------")[0]
ns = {}
exec(head, ns)
phi_s = ns["phi_s"]

OMEGA_50 = 2 * np.pi * 1.4  # rad/us at 50 uT (model's OMEGA)


def phi_at(B_uT, theta_deg):
    return phi_s(theta_deg, omega=OMEGA_50 * B_uT / 50.0)


print("=== Sanity: reproduce the original 50 uT run ===")
p0, p90 = phi_at(50, 0), phi_at(50, 90)
ok = abs(p0 - 0.5075) < 0.01 and abs(p90 - 0.6719) < 0.01
print(f"  Phi_S(0) = {p0:.4f} (orig 0.5075), Phi_S(90) = {p90:.4f} "
      f"(orig 0.6719) -> {'PASS' if ok else 'FAIL'}")

Bs = [25, 35, 45, 50, 55, 65]
thetas = list(range(0, 181, 15))
print("\n=== Phi_S(B, theta) ===")
table = {}
for B in Bs:
    row = [phi_at(B, t) for t in thetas]
    table[str(B)] = {str(t): round(v, 4) for t, v in zip(thetas, row)}
    arr = np.array(row)
    aniso = float(arr.max() - arr.min())
    print(f"  B={B:>2} uT: Phi(0)={row[0]:.4f}  Phi(90)={row[6]:.4f}  "
          f"range [{arr.min():.4f}, {arr.max():.4f}]  "
          f"anisotropy {aniso:.4f} ({100 * aniso / arr.mean():.1f}% of mean)")

results = {
    "model": "radical_pair.py toy model (2 electrons + 1 nucleus), "
             "omega scaled linearly with B; k_rec = 1/us",
    "sanity_50uT": {"phi_0deg": p0, "phi_90deg": p90, "pass": bool(ok)},
    "phi_s_table": table,
    "note": "Toy model: one nucleus, Haberkorn recombination with "
            "k_S = k_T. Yields are model values, not biological data.",
}
with open("radical_pair_fieldsweep_results.json", "w") as f:
    json.dump(results, f, indent=2)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5))
for B in Bs:
    ys = [table[str(B)][str(t)] for t in thetas]
    ax.plot(thetas, ys, "o-", label=f"B = {B} µT")
ax.set_xlabel("field angle θ (degrees)")
ax.set_ylabel("singlet yield Φ_S")
ax.set_title("Radical-pair singlet yield vs angle, across Earth's field range")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("radical_pair_fieldsweep.png", dpi=150)
print("\nSaved radical_pair_fieldsweep_results.json and "
      "radical_pair_fieldsweep.png")
