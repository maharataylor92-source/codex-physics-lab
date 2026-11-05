# NAVIER-STOKES, SIMPLIFIED (Codex index #14; Fluid Dynamics,
# SIM via lattice-Boltzmann D2Q9). Flow past a cylinder in a
# channel. Below Re ~ 47 the wake is steady and symmetric; above
# it, the flow itself starts alternating -- the von Karman vortex
# street. The solver knows nothing about vortices; if they
# appear, Navier-Stokes put them there.
# Cases: Re = 20 (expect steady) and Re = 80 (expect shedding,
# Strouhal St = f D / U in the known cylinder range ~0.15-0.21).
# Gates:
#   G1 mass conserved: mean density drift < 1% in both cases
#   G2 Re=20 wake probe is steady (no spectral peak; tiny RMS)
#   G3 Re=80 probe has a clear shedding peak with St in
#      [0.12, 0.30]

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NX, NY = 220, 72
CX, CY, RAD = 55, 36, 12
D_CYL = 2 * RAD
U0 = 0.05
CXQ = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
CYQ = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
W = np.array([4 / 9] + [1 / 9] * 4 + [1 / 36] * 4)
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])
yy, xx = np.mgrid[0:NY, 0:NX]
SOLID = (xx - CX) ** 2 + (yy - CY) ** 2 <= RAD ** 2
SOLID[0, :] = SOLID[-1, :] = True          # channel walls


def equilibrium(rho, ux, uy):
    feq = np.empty((9, NY, NX))
    usq = ux ** 2 + uy ** 2
    for i in range(9):
        cu = CXQ[i] * ux + CYQ[i] * uy
        feq[i] = rho * W[i] * (1 + 3 * cu + 4.5 * cu ** 2 - 1.5 * usq)
    return feq


def run_case(Re, steps, sample_from):
    """Full LBM with proper bounce-back using a pre-stream copy."""
    nu = U0 * D_CYL / Re
    tau = 3 * nu + 0.5
    rho = np.ones((NY, NX))
    ux = np.full((NY, NX), U0)
    uy = np.zeros((NY, NX))
    ux[SOLID] = 0
    f = equilibrium(rho, ux, uy)
    probe = []
    for n in range(steps):
        rho = f.sum(axis=0)
        ux = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho
        uy = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / rho
        ux[SOLID] = 0
        uy[SOLID] = 0
        feq = equilibrium(rho, ux, uy)
        f = f - (f - feq) / tau
        fpre = f.copy()
        for i in range(9):
            f[i] = np.roll(fpre[i], (CYQ[i], CXQ[i]), axis=(0, 1))
        # bounce-back: populations that streamed INTO solid take
        # the opposite direction's pre-stream value at the solid
        for i in range(9):
            f[i][SOLID] = fpre[OPP[i]][SOLID]
        # inlet (x=0): impose equilibrium with U0, rho from local
        # (column helper: equilibrium() is shaped for full grids;
        # the first launch crashed here on a (NY,) vs (NY,NX)
        # broadcast -- fixed by building the column directly)
        rho0c = f[:, :, 0].sum(axis=0)
        for i in range(9):
            cu = CXQ[i] * U0
            f[i, :, 0] = rho0c * W[i] * (1 + 3 * cu + 4.5 * cu ** 2
                                         - 1.5 * U0 ** 2)
        # outlet: zero-gradient
        f[:, :, -1] = f[:, :, -2]
        if n >= sample_from:
            probe.append(float(uy[CY + 10, CX + 70]))
    vort = ux  # placeholder, vorticity computed in main if needed
    return np.array(probe), float(rho.mean()), tau, (ux, uy)


results = {}
series = {}
for Re, steps, sfrom in [(20, 9000, 6000), (80, 22000, 12000)]:
    probe, rhomean, tau, _ = run_case(Re, steps, sfrom)
    series[Re] = probe
    drift = abs(rhomean - 1.0)
    if len(probe) > 100:
        spec = np.abs(np.fft.rfft(probe - probe.mean()))
        freqs = np.fft.rfftfreq(len(probe), d=1.0)
        k = int(np.argmax(spec[1:]) + 1)
        f_peak = float(freqs[k])
        St = f_peak * D_CYL / U0
        peak_prom = float(spec[k] / (np.median(spec[1:]) + 1e-12))
    else:
        f_peak, St, peak_prom = 0.0, 0.0, 0.0
    results[Re] = {"tau": tau, "rho_drift": float(drift),
                   "probe_rms": float(np.std(probe)),
                   "f_peak_per_step": f_peak, "strouhal": float(St),
                   "peak_prominence": peak_prom}
    print(f"Re={Re}: tau={tau:.3f} rho drift {drift:.4f} probe RMS "
          f"{np.std(probe):.5f} f_peak {f_peak:.5f}/step St {St:.3f} "
          f"prominence {peak_prom:.1f}")

g1 = all(results[r]["rho_drift"] < 0.01 for r in results)
g2 = results[20]["probe_rms"] < 0.004
g3 = 0.12 <= results[80]["strouhal"] <= 0.30 \
    and results[80]["peak_prominence"] > 10
print(f"GATES: G1 mass {g1} | G2 Re20 steady {g2} | G3 Re80 sheds {g3}")

fig, axes = plt.subplots(2, 1, figsize=(9, 4.6), sharex=False)
for ax, Re in zip(axes, (20, 80)):
    ax.plot(series[Re])
    ax.set_ylabel(f"Re={Re}: wake uy")
    ax.set_xlabel("time step (sampling window)")
axes[0].set_title("Wake probe: steady at Re=20, shedding at Re=80 "
                  f"(St = {results[80]['strouhal']:.3f})")
fig.tight_layout()
fig.savefig("navier_stokes.png", dpi=150)

out = {"cases": {str(k): v for k, v in results.items()},
       "gates": {"G1_mass": bool(g1), "G2_steady_low_Re": bool(g2),
                 "G3_shedding_St": bool(g3)}}
with open("navier_stokes_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values()) else "SOME GATE FAILED")
print("Saved navier_stokes_results.json + PNG")
