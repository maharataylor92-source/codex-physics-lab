# MEMBRANE EIGENMODES -- CHLADNI PATTERNS (Beyond the Index #3;
# Wave Physics, SIM). A square membrane clamped at its edges can
# only vibrate in special shapes: the eigenmodes sin(m pi x)
# sin(n pi y), whose squared frequencies form the exact ladder
# lambda_mn = pi^2 (m^2 + n^2) (unit square, c = 1). Sand on a
# real plate collects on the nodal lines -- Chladni figures.
# Here the modes are NOT assumed: the discrete Laplacian's
# smallest eigenvalues are computed numerically (sparse solver),
# then one mode is excited in a time-domain wave simulation and
# its oscillation frequency is measured independently.
# Gates:
#   G1 lowest eigenvalues match pi^2(m^2+n^2) within 2.5%
#      (the discrete Laplacian errs low, O(h^2) -- sign expected)
#   G2 the (1,2)/(2,1) degeneracy holds within 0.5%
#   G3 dynamics: (2,1) mode excited in the wave solver rings at
#      sqrt(lambda_21) within 3%

import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 59                       # interior points per side
h = 1 / (N + 1)
# 2D Laplacian (negative definite form), Kronecker sum
T = sp.diags([np.ones(N - 1), -2 * np.ones(N), np.ones(N - 1)],
             [-1, 0, 1]) / h ** 2
Lap = sp.kron(T, sp.eye(N)) + sp.kron(sp.eye(N), T)
vals, vecs = spla.eigsh(-Lap, k=8, which="SM")
vals = np.sort(vals)
expected = sorted([np.pi ** 2 * (m * m + n * n)
                   for m in (1, 2, 3) for n in (1, 2, 3)])[:8]
print("mode  numeric lambda   pi^2(m^2+n^2)   error")
pairs = [(1, 1), (1, 2), (2, 1), (2, 2), (1, 3), (3, 1), (2, 3), (3, 2)]
for i, ((m, n), v, e) in enumerate(zip(pairs, vals, expected)):
    print(f"({m},{n})  {v:10.3f}      {e:10.3f}     "
          f"{100 * (v - e) / e:+.2f}%")
errs = np.abs(vals - expected) / expected
g1 = bool(np.max(errs) < 0.025)
g2 = abs(vals[1] - vals[2]) / vals[1] < 0.005
print(f"GATES so far: G1 ladder {g1} (max err "
      f"{100 * np.max(errs):.2f}%) | G2 degeneracy {g2}")

# ---- dynamics: excite (2,1), measure the ring frequency ----
xs = np.linspace(h, 1 - h, N)
X, Y = np.meshgrid(xs, xs)
u = np.sin(2 * np.pi * X) * np.sin(1 * np.pi * Y)
u_prev = u.copy()
c, dt = 1.0, 0.4 * h
probe = []
steps = 6000
for n in range(steps):
    lap = (np.roll(u, 1, 0) + np.roll(u, -1, 0)
           + np.roll(u, 1, 1) + np.roll(u, -1, 1) - 4 * u) / h ** 2
    lap[[0, -1], :] = 0
    lap[:, [0, -1]] = 0
    u_next = 2 * u - u_prev + (c * dt) ** 2 * lap
    u_next[[0, -1], :] = 0
    u_next[:, [0, -1]] = 0
    u_prev, u = u, u_next
    probe.append(u[N // 2, N // 4])
# (probe post-mortem: the first run probed (y=N//3, x=N//2);
#  x = 0.5 is a NODE of sin(2 pi x), so it heard only numerical
#  leakage into the (1,1) fundamental, 0.725 Hz, and G3 failed
#  at 35%. The antinode of (2,1) is x = 0.25, y = 0.5.)
probe = np.array(probe)
spec = np.abs(np.fft.rfft(probe - probe.mean()))
freqs = np.fft.rfftfreq(len(probe), d=dt)
f_meas = float(freqs[int(np.argmax(spec[1:]) + 1)])
f_pred = float(np.sqrt(expected[2]) / (2 * np.pi))
print(f"(2,1) mode rings at {f_meas:.3f} Hz; sqrt(lambda)/2pi = "
      f"{f_pred:.3f} Hz (err {100 * abs(f_meas - f_pred) / f_pred:.2f}%)")
g3 = abs(f_meas - f_pred) / f_pred < 0.03
print(f"GATES: G1 {g1} | G2 {g2} | G3 dynamics {g3}")

fig, axes = plt.subplots(1, 4, figsize=(11, 3.0))
for ax, (m, n) in zip(axes, [(1, 1), (1, 2), (2, 2), (1, 3)]):
    Z = np.sin(m * np.pi * X) * np.sin(n * np.pi * Y)
    ax.contour(X, Y, Z, levels=[0], colors="k", linewidths=1.2)
    ax.contourf(X, Y, Z, levels=24, cmap="RdBu_r")
    ax.set_title(f"mode ({m},{n})")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle("Membrane eigenmodes: black lines are the nodes "
             "(where Chladni sand collects)")
fig.tight_layout()
fig.savefig("membrane_modes.png", dpi=150)

out = {"eigenvalues": vals.tolist(), "expected": expected,
       "max_err_frac": float(np.max(errs)),
       "ring_freq_measured": f_meas, "ring_freq_predicted": f_pred,
       "gates": {"G1_ladder": g1, "G2_degeneracy": bool(g2),
                 "G3_dynamics": bool(g3)}}
with open("membrane_modes_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values())
      else "SOME GATE FAILED")
print("Saved membrane_modes_results.json + PNG")
