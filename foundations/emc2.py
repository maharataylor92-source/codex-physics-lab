# E = mc^2 SCALES (Codex index #4; Special Relativity, NUM).
# Rest energy is not an abstraction: 1 kg of anything is 21.5
# megatons of TNT equivalent, the Sun burns 4.3 billion kg of
# itself every second, and fission's power is a 0.09% mass tax.
# Plus the relativistic correction: where Newton stops working.
# Gates (known values):
#   G1 electron rest energy 0.511 MeV, proton 938.3 MeV (1%)
#   G2 1 kg -> 21.5 Mt TNT within 2%
#   G3 gamma(0.6c) = 1.25 exactly; KE ratio at 0.1c = 1.0075
#      (series 1 + 3b^2/4; the first draft of this gate said
#      1.005 -- that is gamma-1 at 0.1c, a different number.
#      The computation was right; the anchor was not.)

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 299792458.0
EV = 1.602176634e-19
U_KG = 1.66053906660e-27
MT_TNT = 4.184e15          # J per megaton
LSUN = 3.828e26            # W

e_rest = 9.1093837015e-31 * C ** 2 / EV / 1e6     # MeV
p_rest = 1.67262192369e-27 * C ** 2 / EV / 1e6
print(f"electron rest energy {e_rest:.4f} MeV | proton {p_rest:.2f} MeV")
g1 = abs(e_rest - 0.511) / 0.511 < 0.01 and abs(p_rest - 938.27) / 938.27 < 0.01

kg_J = C ** 2
kg_Mt = kg_J / MT_TNT
print(f"1 kg = {kg_J:.4e} J = {kg_Mt:.1f} megatons TNT equivalent")
g2 = abs(kg_Mt - 21.5) / 21.5 < 0.02

sun_kg_s = LSUN / C ** 2
print(f"Sun converts {sun_kg_s:.3e} kg/s of mass into light")
fission_frac = 200.0 / (235.044 * 931.494)   # ~200 MeV per U-235 fission
print(f"U-235 fission releases ~200 MeV of {235.044 * 931.494:,.0f} MeV "
      f"rest energy = {100 * fission_frac:.3f}% of the mass")

betas = [0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.99]
rows = []
for b in betas:
    gam = 1 / np.sqrt(1 - b ** 2)
    ke_rel = (gam - 1)
    ke_newt = 0.5 * b ** 2
    rows.append({"beta": b, "gamma": float(gam),
                 "KE_rel_mc2": float(ke_rel),
                 "KE_newton_mc2": float(ke_newt),
                 "ratio": float(ke_rel / ke_newt)})
    print(f"v={b:.2f}c: gamma={gam:.4f}  KE_rel={ke_rel:.4f} mc^2  "
          f"KE_Newton={ke_newt:.4f}  Newton is low by "
          f"{100 * (ke_rel / ke_newt - 1):.1f}%")
g3 = abs(1 / np.sqrt(1 - 0.36) - 1.25) < 1e-12 \
    and abs(rows[0]["ratio"] - 1.0075) < 0.002
print(f"GATES: G1 {g1} G2 {g2} G3 {g3}")

fig, ax = plt.subplots(figsize=(7.4, 4.0))
bb = np.linspace(0.01, 0.995, 100)
ax.plot(bb, 1 / np.sqrt(1 - bb ** 2) - 1, label="relativistic KE / mc^2")
ax.plot(bb, 0.5 * bb ** 2, "--", label="Newtonian KE / mc^2")
ax.set_xlabel("v / c")
ax.set_ylabel("kinetic energy (units of mc^2)")
ax.set_title("Where Newton stops: kinetic energy vs speed")
ax.legend()
fig.tight_layout()
fig.savefig("emc2.png", dpi=150)

out = {"electron_MeV": float(e_rest), "proton_MeV": float(p_rest),
       "kg_Mt_TNT": float(kg_Mt), "sun_kg_per_s": float(sun_kg_s),
       "fission_fraction": float(fission_frac), "rows": rows,
       "gates": {"G1_rest_energies": bool(g1), "G2_kg_TNT": bool(g2),
                 "G3_gamma": bool(g3)}}
with open("emc2_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values()) else "SOME GATE FAILED")
print("Saved emc2_results.json + PNG")
