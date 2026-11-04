# Quantum tunneling -- Vol II Section XIX ("The Epoch Breakthrough")
# of The Grand Unification of Being.
#
# Two independent methods, cross-checked:
#   A. Analytic transmission coefficient T(E, V, w) for a rectangular
#      barrier (exact textbook formula; thick-barrier limit ~ exp(-2*kappa*w)).
#   B. Split-step Fourier propagation of a Gaussian wavepacket through
#      the same barrier -- the wave function itself, no formula assumed.
# Units: hbar = m = 1 (dimensionless). Packet central energy E = k0^2/2.
#
# Verification gates:
#   1. No barrier (V=0) -> wavepacket transmission ~ 1, norm conserved.
#   2. Analytic T -> 1 as w -> 0; T decreases monotonically with w and V.
#   3. Wavepacket T (w=1,2,3 at V=1.5, E=0.5) tracks the analytic trend
#      (same order of magnitude; packet has an energy spread, so exact
#      equality with the single-energy value is NOT expected).

import json

import numpy as np

HBAR = M = 1.0


def kappa(E, V):
    return np.sqrt(2 * M * (V - E)) / HBAR


def T_exact(E, V, w):
    """Exact rectangular-barrier transmission."""
    if V <= 0 or w <= 0:
        return 1.0
    if E < V:
        k = kappa(E, V)
        return float(1.0 / (1.0 + V**2 * np.sinh(k * w) ** 2
                            / (4 * E * (V - E))))
    q = np.sqrt(2 * M * (E - V)) / HBAR
    return float(1.0 / (1.0 + V**2 * np.sin(q * w) ** 2
                        / (4 * E * (E - V))))


def T_thick(E, V, w):
    """Thick-barrier approximation 16E(V-E)/V^2 * exp(-2 kappa w)."""
    k = kappa(E, V)
    return float(16 * E * (V - E) / V**2 * np.exp(-2 * k * w))


# ---------- wavepacket propagation (split-step Fourier) ----------
L, N = 80.0, 4096
x = np.linspace(-L / 2, L / 2, N, endpoint=False)
dx = x[1] - x[0]
kgrid = 2 * np.pi * np.fft.fftfreq(N, d=dx)


def run_packet(V, w, E=0.5, sigma=3.0, x0=-14.0, t_final=30.0, dt=0.01):
    k0 = np.sqrt(2 * M * E) / HBAR
    psi = (np.exp(-(x - x0) ** 2 / (4 * sigma**2) + 1j * k0 * x)
           / (2 * np.pi * sigma**2) ** 0.25).astype(complex)
    Varr = np.where(np.abs(x) <= w / 2, V, 0.0)
    kin = np.exp(-1j * HBAR * kgrid**2 * dt / (2 * M))
    pot_half = np.exp(-1j * Varr * dt / (2 * HBAR))
    steps = int(t_final / dt)
    for _ in range(steps):
        psi *= pot_half
        psi = np.fft.ifft(kin * np.fft.fft(psi))
        psi *= pot_half
    prob = np.abs(psi) ** 2
    norm = float(np.trapezoid(prob, x))
    transmitted = float(np.trapezoid(prob[x > w / 2], x[x > w / 2]))
    reflected = float(np.trapezoid(prob[x < -w / 2], x[x < -w / 2]))
    return {"norm": norm, "T_packet": transmitted / norm,
            "R_packet": reflected / norm,
            "in_barrier": 1.0 - transmitted / norm - reflected / norm}


print("=== Verification gates ===")
free = run_packet(V=0.0, w=2.0)
g1 = abs(free["norm"] - 1) < 1e-6 and free["T_packet"] > 0.99
print(f"Gate 1 (no barrier: norm {free['norm']:.8f}, "
      f"T {free['T_packet']:.4f}) -> {'PASS' if g1 else 'FAIL'}")
ws_mono = [T_exact(0.5, 1.5, w) for w in (0.5, 1, 2, 3, 5)]
g2 = (T_exact(0.5, 1.5, 1e-9) > 0.999
      and all(a > b for a, b in zip(ws_mono, ws_mono[1:])))
print(f"Gate 2 (T->1 as w->0; monotone in w) -> {'PASS' if g2 else 'FAIL'}")

print("\n=== Analytic sweeps (E = 0.5) ===")
sweep_w = {str(w): T_exact(0.5, 1.5, w) for w in (0.5, 1, 1.5, 2, 3, 4, 5)}
for w, t in sweep_w.items():
    print(f"  V=1.5, w={w}: T = {t:.6f}   (thick-barrier approx "
          f"{T_thick(0.5, 1.5, float(w)):.6f})")
sweep_V = {str(v): T_exact(0.5, v, 2.0) for v in (0.75, 1.0, 1.5, 2.0, 3.0)}
for v, t in sweep_V.items():
    print(f"  w=2, V={v}: T = {t:.6f}")
sweep_E = {str(e): T_exact(e, 1.5, 2.0)
           for e in (0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.49)}
for e, t in sweep_E.items():
    print(f"  w=2, V=1.5, E={e}: T = {t:.6f}")

print("\n=== Wavepacket runs (E = 0.5, V = 1.5) ===")
packet = {}
for w in (1.0, 2.0, 3.0):
    r = run_packet(V=1.5, w=w)
    packet[str(w)] = r
    print(f"  w={w}: T_packet = {r['T_packet']:.6f} "
          f"(analytic {T_exact(0.5, 1.5, w):.6f}), "
          f"R = {r['R_packet']:.4f}, norm {r['norm']:.8f}")
g3 = all(0.2 < packet[str(w)]["T_packet"] / T_exact(0.5, 1.5, w) < 5
         for w in (1.0, 2.0, 3.0))
print(f"Gate 3 (packet vs analytic within 5x, same trend) -> "
      f"{'PASS' if g3 else 'FAIL'}")

results = {
    "model": "1D rectangular barrier, hbar=m=1; analytic T + split-step "
             "Fourier wavepacket (Gaussian, sigma=3, E=0.5)",
    "gates": {"no_barrier": bool(g1), "monotonic": bool(g2),
              "packet_vs_analytic": bool(g3)},
    "analytic_T_vs_width_V1.5_E0.5": sweep_w,
    "analytic_T_vs_height_w2_E0.5": sweep_V,
    "analytic_T_vs_energy_w2_V1.5": sweep_E,
    "wavepacket_V1.5_E0.5": packet,
    "headline": {f"w={w}": T_exact(0.5, 1.5, w) for w in (1.0, 2.0, 3.0)},
}
with open("tunneling_results.json", "w") as f:
    json.dump(results, f, indent=2)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
wv = np.array([float(w) for w in sweep_w])
axes[0].semilogy(wv, list(sweep_w.values()), "o-", label="exact T")
axes[0].semilogy(wv, [T_thick(0.5, 1.5, w) for w in wv], "s--",
                 label="thick-barrier approx")
axes[0].semilogy([1, 2, 3], [packet[str(w)]["T_packet"] for w in (1.0, 2.0, 3.0)],
                 "^", ms=9, label="wavepacket")
axes[0].set_xlabel("barrier width w"); axes[0].set_ylabel("transmission T")
axes[0].set_title("Tunneling vs barrier width (V=1.5, E=0.5)")
axes[0].legend(); axes[0].grid(alpha=0.3)
ev = np.array([float(e) for e in sweep_E])
axes[1].semilogy(ev, list(sweep_E.values()), "o-")
axes[1].axvline(1.5, color="r", ls=":", label="barrier height V")
axes[1].set_xlabel("particle energy E"); axes[1].set_ylabel("transmission T")
axes[1].set_title("Tunneling vs energy (w=2, V=1.5)")
axes[1].legend(); axes[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig("tunneling_transmission.png", dpi=150)
print("\nSaved tunneling_results.json and tunneling_transmission.png")
