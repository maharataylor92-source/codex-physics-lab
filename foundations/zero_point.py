# ZERO-POINT ENERGY (Codex index #21; Quantum, NUM).
# E_n = (n + 1/2) hbar omega. The ground state of an oscillator
# is NOT zero: E_0 = (1/2) hbar omega, and it is not removable
# (it is the uncertainty floor of index #20, wearing a number).
# Scale table across real oscillators; the comparison energy is
# k_B T at 300 K (0.02585 eV): wherever E_0 >> k_B T, the
# oscillator lives in its ground state at room temperature and
# its zero-point motion is a room-temperature fact.
# Gates (arithmetic/known-value checks, stated as such):
#   G1 H2 molecule: E_0 within 5% of the known 0.27 eV
#   G2 level spacing E_1 - E_0 = hbar omega exactly in the table
#   G3 every molecular row's level SPACING hbar*omega = 2 E_0
#      exceeds 5 k_B T (first excited state thermally empty at
#      300 K -- the correct freeze criterion is the spacing, not
#      E_0 itself; O2's E_0 alone is only 3.8 kT)

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HBAR = 1.054571817e-34
EV = 1.602176634e-19
KB_T_300 = 1.380649e-23 * 300 / EV      # 0.025852 eV

# name, vibrational wavenumber (cm^-1) for molecules, or angular
# frequency omega (rad/s) directly; source: standard spectroscopy
ROWS = []
SPECS = [
    ("H2 molecule", 4401.2),
    ("N2 molecule", 2358.6),
    ("O2 molecule", 1580.2),
    ("CO molecule", 2169.8),
    ("H2O bend", 1595.0),
    ("C-H stretch", 2960.0),
]
for name, wn in SPECS:
    omega = 2 * np.pi * 2.99792458e10 * wn      # rad/s
    E0 = 0.5 * HBAR * omega / EV
    ROWS.append({"system": name, "omega_rad_s": float(omega),
                 "E0_eV": float(E0), "E1_eV": float(3 * E0),
                 "spacing_eV": float(2 * E0),
                 "E0_over_kT300": float(E0 / KB_T_300)})
    print(f"{name:<14}: omega={omega:.3e} rad/s  E0={E0:.4f} eV  "
          f"(= {E0 / KB_T_300:.1f} kT at 300 K)")

# macroscopic contrast: a 1 Hz pendulum-scale oscillator and a
# quartz watch crystal (32.768 kHz)
for name, f_hz in [("pendulum 1 Hz", 1.0),
                   ("quartz crystal 32.768 kHz", 32768.0),
                   ("quartz oscillator 10 MHz", 1.0e7)]:
    omega = 2 * np.pi * f_hz
    E0 = 0.5 * HBAR * omega / EV
    ROWS.append({"system": name, "omega_rad_s": float(omega),
                 "E0_eV": float(E0), "E1_eV": float(3 * E0),
                 "spacing_eV": float(2 * E0),
                 "E0_over_kT300": float(E0 / KB_T_300)})
    print(f"{name:<26}: E0={E0:.3e} eV (= {E0 / KB_T_300:.2e} kT)")

h2 = ROWS[0]
g1 = abs(h2["E0_eV"] - 0.27) / 0.27 < 0.05
g2 = all(abs(r["E1_eV"] - r["E0_eV"] - r["spacing_eV"]) < 1e-12
         for r in ROWS)
g3 = all(r["spacing_eV"] / KB_T_300 > 5 for r in ROWS[:6])
print(f"G1 H2 E0 ~ 0.27 eV: {g1} | G2 spacing identity: {g2} | "
      f"G3 molecular spacings > 5 kT at 300 K: {g3}")

fig, ax = plt.subplots(figsize=(8.4, 3.8))
names = [r["system"] for r in ROWS]
vals = [r["E0_eV"] for r in ROWS]
ax.barh(names, np.log10(vals))
ax.axvline(np.log10(KB_T_300), color="r", ls="--",
           label="kT at 300 K")
ax.set_xlabel("log10 E_0 (eV)")
ax.set_title("Zero-point energy across oscillators")
ax.legend()
fig.tight_layout()
fig.savefig("zero_point.png", dpi=150)

out = {"kT_300K_eV": KB_T_300, "rows": ROWS,
       "gates": {"G1_H2_known_value": bool(g1),
                 "G2_spacing_identity": bool(g2),
                 "G3_molecules_frozen": bool(g3)}}
with open("zero_point_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values()) else "SOME GATE FAILED")
print("Saved zero_point_results.json")
