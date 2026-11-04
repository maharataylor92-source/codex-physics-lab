# UNCERTAINTY WAVEPACKET (Codex index #20; Quantum Foundations,
# NUM/SIM). A free Gaussian wavepacket spreads; Heisenberg says
# dx*dp >= hbar/2 at every instant. With hbar = m = 1 and initial
# width s0: dx(t) = s0 sqrt(1 + (t/(2 s0^2))^2), dp = 1/(2 s0)
# constant, so the product starts AT the floor (a Gaussian is a
# minimum-uncertainty state) and rises. Split-step FFT evolution
# on a grid, moments measured from |psi|^2 and |phi(p)|^2.
# Gates:
#   G1 measured dx(t) tracks the analytic curve within 2%
#   G2 product never dips below hbar/2 (0.5) by more than 0.5%
#   G3 dp constant within 1% (free particle: momentum frozen)

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N, L = 4096, 240.0
x = np.linspace(-L / 2, L / 2, N, endpoint=False)
dx = x[1] - x[0]
p = 2 * np.pi * np.fft.fftfreq(N, d=dx)
s0 = 3.0
psi = (2 * np.pi * s0 ** 2) ** (-0.25) * np.exp(-x ** 2 / (4 * s0 ** 2))
psi = psi.astype(complex)
dt = 0.05
T = 60.0
steps = int(T / dt)
kin = np.exp(-0.5j * p ** 2 * dt)


def moments(psi):
    prob = np.abs(psi) ** 2 * dx
    mx = np.sum(x * prob)
    dx_ = np.sqrt(np.sum((x - mx) ** 2 * prob))
    phi = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(psi))) * dx \
        / np.sqrt(2 * np.pi)
    ps = np.fft.fftshift(p)
    pp = np.abs(phi) ** 2 * (ps[1] - ps[0])
    pp /= pp.sum()
    mp = np.sum(ps * pp)
    dp_ = np.sqrt(np.sum((ps - mp) ** 2 * pp))
    return dx_, dp_


rows = []
check_every = int(5.0 / dt)
for n in range(steps + 1):
    t = n * dt
    if n % check_every == 0:
        dxm, dpm = moments(psi)
        dxa = s0 * np.sqrt(1 + (t / (2 * s0 ** 2)) ** 2)
        rows.append({"t": float(t), "dx_meas": float(dxm),
                     "dx_analytic": float(dxa), "dp_meas": float(dpm),
                     "product": float(dxm * dpm)})
        print(f"t={t:5.1f}: dx {dxm:.3f} (analytic {dxa:.3f})  "
              f"dp {dpm:.4f}  product {dxm * dpm:.4f}")
    if n < steps:
        psi = np.fft.ifft(kin * np.fft.fft(psi))

err = max(abs(r["dx_meas"] - r["dx_analytic"]) / r["dx_analytic"]
          for r in rows)
prod_min = min(r["product"] for r in rows)
dp_drift = (max(r["dp_meas"] for r in rows)
            - min(r["dp_meas"] for r in rows)) / rows[0]["dp_meas"]
gates = {"G1_dx_tracks_analytic_2pct": err < 0.02,
         "G2_product_never_below_floor": prod_min >= 0.5 * 0.995,
         "G3_dp_constant_1pct": dp_drift < 0.01}
print(f"max dx error {100 * err:.2f}% | min product {prod_min:.4f} "
      f"| dp drift {100 * dp_drift:.2f}%")
print(gates)

fig, ax = plt.subplots(figsize=(7.2, 3.8))
ts = [r["t"] for r in rows]
ax.plot(ts, [r["dx_meas"] for r in rows], "o-", label="dx measured")
ax.plot(ts, [r["dx_analytic"] for r in rows], "--", label="dx analytic")
ax2 = ax.twinx()
ax2.plot(ts, [r["product"] for r in rows], "s-", color="tab:red",
         label="dx*dp")
ax2.axhline(0.5, color="tab:red", ls=":", label="hbar/2 floor")
ax.set_xlabel("time (hbar = m = 1)")
ax.set_ylabel("dx")
ax2.set_ylabel("dx * dp")
ax.legend(loc="upper left", fontsize=8)
ax2.legend(loc="lower right", fontsize=8)
ax.set_title("Free Gaussian packet: spreading on the uncertainty floor")
fig.tight_layout()
fig.savefig("uncertainty_wavepacket.png", dpi=150)

out = {"rows": rows, "max_dx_rel_err": float(err),
       "min_product": float(prod_min), "dp_drift": float(dp_drift),
       "gates": {k: bool(v) for k, v in gates.items()}}
with open("uncertainty_wavepacket_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(gates.values()) else "SOME GATE FAILED")
print("Saved uncertainty_wavepacket_results.json + PNG")
