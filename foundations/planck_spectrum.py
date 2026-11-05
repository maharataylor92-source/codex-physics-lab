# PLANCK'S BLACKBODY SPECTRUM (Beyond the Index #8; Thermal
# Radiation, NUM). One formula -- Planck's law -- must pay out
# two older laws at once:
#   * Stefan-Boltzmann: total radiated power = sigma T^4,
#     sigma = 5.670374419e-8 W m^-2 K^-4
#   * Wien displacement: lambda_peak * T = 2.897771955e-3 m K
# Both are EXTRACTED from the numerically integrated spectrum
# here, not assumed: sigma falls out of the area under the
# curve, Wien's constant out of its peak position.
# Gates:
#   G1 recovered sigma within 1% at every temperature
#   G2 recovered Wien constant within 0.5% at every T
#   G3 power ratio M(2T)/M(T) = 16 within 1% (the T^4 law
#      itself, checked as a scaling, not a constant)

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = 6.62607015e-34
C = 299792458.0
KB = 1.380649e-23
SIGMA_EXACT = 5.670374419e-8
WIEN_EXACT = 2.897771955e-3

lam = np.logspace(-8, -2, 400000)          # 10 nm .. 1 cm


def spectrum(T):
    x = H * C / (lam * KB * T)
    B = (2 * H * C ** 2 / lam ** 5) / np.expm1(np.minimum(x, 700))
    return B                                  # W sr^-1 m^-3


rows = []
for T in (300.0, 1000.0, 3000.0, 5778.0):
    B = spectrum(T)
    M = np.pi * np.trapezoid(B, lam)          # W/m^2 (Lambert)
    sigma_rec = M / T ** 4
    lam_peak = float(lam[int(np.argmax(B))])
    wien_rec = lam_peak * T
    rows.append({"T": T, "sigma_recovered": float(sigma_rec),
                 "sigma_err_frac": float(abs(sigma_rec - SIGMA_EXACT)
                                         / SIGMA_EXACT),
                 "wien_recovered": float(wien_rec),
                 "wien_err_frac": float(abs(wien_rec - WIEN_EXACT)
                                        / WIEN_EXACT),
                 "M": float(M)})
    print(f"T = {T:6.0f} K: sigma = {sigma_rec:.6e} "
          f"(err {rows[-1]['sigma_err_frac'] * 100:.3f}%), "
          f"lambda_peak*T = {wien_rec:.6e} "
          f"(err {rows[-1]['wien_err_frac'] * 100:.3f}%)")

# M(3000)/M(1000) should be 3^4 = 81:
ratio = rows[2]['M'] / rows[1]['M']
print(f"M(3000)/M(1000) = {ratio:.2f} (T^4 law predicts 81)")
g1 = all(r["sigma_err_frac"] < 0.01 for r in rows)
g2 = all(r["wien_err_frac"] < 0.005 for r in rows)
g3 = abs(ratio - 81) / 81 < 0.01
print(f"GATES: G1 sigma {g1} | G2 Wien {g2} | G3 T^4 {g3}")

fig, ax = plt.subplots(figsize=(7.6, 4.2))
for T in (300.0, 1000.0, 3000.0, 5778.0):
    B = spectrum(T)
    ax.plot(lam * 1e6, np.pi * B, label=f"{T:.0f} K")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("wavelength (um)")
ax.set_ylabel("exitance per unit wavelength (W m^-2 m^-1)")
ax.set_title("Planck spectra: one curve, two laws inside it")
ax.legend()
fig.tight_layout()
fig.savefig("planck_spectrum.png", dpi=150)

out = {"rows": rows, "sigma_exact": SIGMA_EXACT,
       "wien_exact": WIEN_EXACT, "T4_ratio": float(ratio),
       "gates": {"G1_sigma": bool(g1), "G2_wien": bool(g2),
                 "G3_T4": bool(g3)}}
with open("planck_spectrum_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values())
      else "SOME GATE FAILED")
print("Saved planck_spectrum_results.json + PNG")
