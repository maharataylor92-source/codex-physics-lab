# Heat flow -- Codex Vol I, Thermodynamics: Fourier's law q = -k * grad(T)
# (heat equation dT/dt = alpha * d2T/dx2, alpha = k/(rho*c)) with the
# thermal-energy law DQ = m*c*DT as the energy-balance gate.
#
# Setup: 1D rod, left end held at 100 C, right end held at 0 C, interior
# starts at 0 C. Explicit FTCS integration (stability alpha*dt/dx^2 <= 0.5).
#
# Gates:
#   1. Steady state is the exact linear profile; flux q is uniform.
#   2. Transient center temperature matches the analytic Fourier-series
#      solution (50 terms) within 1.5 C at sampled times.
#   3. Energy balance: heat stored in the rod (integral of rho*c*T dx)
#      equals time-integrated net boundary flux (Fourier at both ends)
#      within 2% -- this is DQ = mcDT checked against q = -k grad T.
#   4. Materials order by diffusivity: copper equilibrates fastest,
#      then aluminium, then iron (alpha values from standard tables).

import json

import numpy as np

L = 0.1                      # rod length, m
N = 101
x = np.linspace(0, L, N)
dx = x[1] - x[0]
T_HOT, T_COLD = 100.0, 0.0

# material: (k W/mK, rho kg/m3, c J/kgK)
MATERIALS = {
    "copper": (401.0, 8960.0, 385.0),
    "aluminium": (237.0, 2700.0, 900.0),
    "iron": (80.0, 7870.0, 450.0),
}


def analytic_center(t, alpha, terms=50):
    """Center temperature from the Fourier-series solution."""
    # T(x,t) = T_ss(x) + series; T_ss linear 100 -> 0
    s = 0.0
    for n in range(1, terms + 1):
        # initial deviation from steady state: -T_ss(x)
        bn = (2 / L) * np.trapezoid(
            -T_HOT * (1 - x / L) * np.sin(n * np.pi * x / L), x)
        s += (bn * np.sin(n * np.pi * 0.5)
              * np.exp(-alpha * (n * np.pi / L) ** 2 * t))
    return T_HOT * 0.5 + s


def run_rod(alpha, t_final, safety=0.45):
    dt = safety * dx**2 / alpha
    steps = int(t_final / dt)
    T = np.zeros(N)
    T[0], T[-1] = T_HOT, T_COLD
    r = alpha * dt / dx**2
    snaps, center = {}, []
    snap_times = [t_final * f for f in (0.05, 0.2, 0.5, 1.0)]
    stored_flux = 0.0
    for i in range(steps):
        Tn = T.copy()
        Tn[1:-1] = T[1:-1] + r * (T[2:] - 2 * T[1:-1] + T[:-2])
        T = Tn
        # net influx per unit area: q_in(left) - q_out(right), q=-k dT/dx
        # in alpha form with rho*c = k/alpha: flux/area = -(k) dT/dx
        q_left = -(1.0) * (T[1] - T[0]) / dx      # units of k * C/m
        q_right = -(1.0) * (T[-1] - T[-2]) / dx
        stored_flux += (q_left - q_right) * dt    # k * C * s / m
        t_now = (i + 1) * dt
        center.append((t_now, T[N // 2]))
        for st in snap_times:
            if abs(t_now - st) < dt / 2:
                snaps[round(st, 3)] = T.copy()
    return T, snaps, center, stored_flux, dt


print("=== Gates (copper rod, L = 0.1 m) ===")
k_cu, rho_cu, c_cu = MATERIALS["copper"]
alpha_cu = k_cu / (rho_cu * c_cu)
tau = L**2 / alpha_cu
T_fin, snaps, center, flux_int, dt = run_rod(alpha_cu, 3 * tau)
ss_exact = T_HOT * (1 - x / L)
err_ss = float(np.max(np.abs(T_fin - ss_exact)))
q_prof = -k_cu * np.gradient(T_fin, dx)
q_spread = float((q_prof[5:-5].max() - q_prof[5:-5].min())
                 / abs(q_prof[5:-5].mean()))
g1 = err_ss < 0.5 and q_spread < 0.02
print(f"Gate 1 (steady state linear, uniform flux): max profile err "
      f"{err_ss:.3f} C, flux spread {100 * q_spread:.2f}% -> "
      f"{'PASS' if g1 else 'FAIL'}")

# transient check at t = 0.2 tau
t_chk = 0.2 * tau
_, _, center2, _, _ = run_rod(alpha_cu, t_chk)
num_c = center2[-1][1]
ana_c = analytic_center(t_chk, alpha_cu)
g2 = abs(num_c - ana_c) < 1.5
print(f"Gate 2 (center T at 0.2 tau): numeric {num_c:.2f} C vs series "
      f"{ana_c:.2f} C -> {'PASS' if g2 else 'FAIL'}")

# energy balance at t = tau (mid-transient)
_, _, _, flux_tau, _ = run_rod(alpha_cu, tau)
T_tau, _, _, _, _ = run_rod(alpha_cu, tau)
stored = rho_cu * c_cu * float(np.trapezoid(T_tau, x))   # J/m2
flux_J = k_cu * flux_tau                                  # J/m2
bal = abs(stored - flux_J) / stored
g3 = bal < 0.02
print(f"Gate 3 (DQ = mcDT vs integrated Fourier flux): stored "
      f"{stored:.4g} J/m2, flux {flux_J:.4g} J/m2, diff "
      f"{100 * bal:.2f}% -> {'PASS' if g3 else 'FAIL'}")

print("\n=== Materials: time for the rod center to reach 45 C ===")
# (45 C, not 50: the center tends to exactly 50 C only as t -> infinity,
# since 50 is the steady-state midpoint temperature.)
mat_out = {}
curves = {}
for name, (k, rho, c) in MATERIALS.items():
    alpha = k / (rho * c)
    t45 = None
    _, _, cen, _, _ = run_rod(alpha, 4 * L**2 / alpha)
    curves[name] = cen
    for tt, Tc in cen:
        if Tc >= 45.0:
            t45 = tt
            break
    mat_out[name] = {"k": k, "alpha_m2_s": alpha, "t_center_45C_s": t45}
    print(f"  {name:<10} alpha = {alpha:.3e} m2/s, "
          f"center hits 45 C at t = {t45:.1f} s")
g4 = (mat_out["copper"]["t_center_45C_s"]
      < mat_out["aluminium"]["t_center_45C_s"]
      < mat_out["iron"]["t_center_45C_s"])
print(f"Gate 4 (ordering copper < aluminium < iron) -> "
      f"{'PASS' if g4 else 'FAIL'}")

out = {"model": "1D heat equation (Fourier), FTCS, rod 0.1 m, ends "
                "100 C / 0 C",
       "gates": {"steady_state": bool(g1), "series_match": bool(g2),
                 "energy_balance": bool(g3), "material_order": bool(g4)},
       "materials": mat_out,
       "steady_flux_W_m2": float(abs(q_prof.mean()))}
with open("heat_flow_results.json", "w") as f:
    json.dump(out, f, indent=2)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for st, prof in snaps.items():
    axes[0].plot(x * 100, prof, label=f"t = {st:.0f} s")
axes[0].plot(x * 100, ss_exact, "k:", label="steady state")
axes[0].set_xlabel("position along rod (cm)")
axes[0].set_ylabel("temperature (°C)")
axes[0].set_title("Copper rod: temperature profiles")
axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)
for name, cen in curves.items():
    tt = np.array([p[0] for p in cen]); Tc = np.array([p[1] for p in cen])
    axes[1].plot(tt, Tc, label=name)
axes[1].axhline(45, color="gray", ls=":")
axes[1].set_xlabel("time (s)"); axes[1].set_ylabel("center temperature (°C)")
axes[1].set_title("Center warming by material")
axes[1].legend(); axes[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig("heat_flow.png", dpi=150)
print("\nSaved heat_flow_results.json and heat_flow.png")
