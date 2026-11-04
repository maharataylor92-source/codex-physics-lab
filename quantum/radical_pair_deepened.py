# RADICAL PAIR DEEPENED -- robustness of the section XVI compass.
# Uses the ORIGINAL model's exact conventions (radical_pair.py):
# anisotropic hyperfine tensor (a_perp/2pi = 0.5 MHz, a_par/2pi =
# 1.5 MHz on electron 1), omega = 2pi*1.4 rad/us at 50 uT, singlet
# birth, nuclei unpolarized, k = 1/us, t grid 0-12 us.
#
# (A) MORE NUCLEI: add nuclei 2 and 3, each with the same tensor
#     shape scaled down (0.3/0.9 and 0.2/0.6 MHz). Does the
#     anisotropy survive contact with a busier spin bath?
# (B) UNEQUAL RATES: Haberkorn dynamics with k_S != k_T, solved
#     exactly through the Liouvillian: Phi_S = k_S Tr[P_S (-L^-1
#     rho0)]. The original model set k_S = k_T = 1.
#
# Gates:
#   G1 part A at n=1 reproduces the original: yields 0.5075/0.6719,
#      anisotropy 0.174
#   G2 part B at k_T = k_S reproduces the same numbers
#   G3 (report) anisotropy trend across n and across k_T/k_S

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
I2 = np.eye(2, dtype=complex)
SPIN = (sx, sy, sz)
OMEGA = 2 * np.pi * 1.4
T_GRID = np.linspace(0, 12, 1201)
K = 1.0
# (a_perp, a_par) per nucleus, MHz x 2pi
TENSORS = [(2 * np.pi * 0.5, 2 * np.pi * 1.5),
           (2 * np.pi * 0.3, 2 * np.pi * 0.9),
           (2 * np.pi * 0.2, 2 * np.pi * 0.6)]


def kron_all(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def factor_op(factor, comp, nf):
    ops = [I2] * nf
    ops[factor] = SPIN[comp]
    return kron_all(ops)


def build_H(theta, n_nuc):
    nf = 2 + n_nuc
    n = np.array([np.sin(theta), 0.0, np.cos(theta)])
    H = sum(OMEGA * n[j] * (factor_op(0, j, nf) + factor_op(1, j, nf))
            for j in range(3))
    for k in range(n_nuc):
        a_perp, a_par = TENSORS[k]
        for c in range(3):
            a = a_par if c == 2 else a_perp
            ops = [I2] * nf
            ops[0] = SPIN[c]
            ops[2 + k] = SPIN[c]
            H = H + a * kron_all(ops)
    return H


def singlet_PS(nf):
    s = np.zeros(4, dtype=complex)
    s[1] = 1 / np.sqrt(2)
    s[2] = -1 / np.sqrt(2)
    return kron_all([np.outer(s, s.conj())] + [I2] * (nf - 2))


def yield_equal(theta, n_nuc):
    nf = 2 + n_nuc
    dim = 2 ** nf
    H = build_H(theta, n_nuc)
    PS = singlet_PS(nf)
    evals, V = np.linalg.eigh(H)
    s = np.zeros(4, dtype=complex)
    s[1] = 1 / np.sqrt(2)
    s[2] = -1 / np.sqrt(2)
    rho0 = kron_all([np.outer(s, s.conj())]
                    + [I2 / 2] * n_nuc)          # singlet x mixed nuclei
    M = (V.conj().T @ PS @ V) * (V.conj().T @ rho0 @ V).T
    dE = evals[:, None] - evals[None, :]
    ph = np.exp(-1j * T_GRID[:, None, None] * dE[None])
    ps_t = np.real(np.einsum("ij,tij->t", M, ph))
    w = np.exp(-K * T_GRID)
    return float(np.trapezoid(w * ps_t, T_GRID) * K)


print("=== Part A: anisotropy vs number of nuclei ===")
angles_deg = [0, 15, 30, 45, 60, 75, 90]
partA = {}
for n in (1, 2, 3):
    ys = [yield_equal(np.radians(a), n) for a in angles_deg]
    aniso = max(ys) - min(ys)
    partA[n] = {"yields": [float(y) for y in ys],
                "anisotropy": float(aniso),
                "aniso_fraction": float(aniso / np.mean(ys))}
    print(f"  n={n}: min {min(ys):.4f} max {max(ys):.4f} "
          f"anisotropy {aniso:.4f} ({100 * aniso / np.mean(ys):.1f}%)")
g1 = (abs(partA[1]["yields"][0] - 0.5075) < 0.005
      and abs(partA[1]["yields"][-1] - 0.6719) < 0.005)
print(f"G1 n=1 reproduces original yields/anisotropy: {g1}")

print("=== Part B: unequal rates k_T/k_S (1 nucleus) ===")
nf, dim = 3, 8
H0 = build_H(0.0, 1)
PS = singlet_PS(nf)
PT = np.eye(dim) - PS
s = np.zeros(4, dtype=complex)
s[1] = 1 / np.sqrt(2)
s[2] = -1 / np.sqrt(2)
rho0 = kron_all([np.outer(s, s.conj()), I2 / 2])
Id = np.eye(dim)


def yield_unequal(theta, kS, kT):
    H = build_H(theta, 1)
    L = (-1j * (np.kron(H, Id) - np.kron(Id, H.T))
         - (kS / 2) * (np.kron(PS, Id) + np.kron(Id, PS.T))
         - (kT / 2) * (np.kron(PT, Id) + np.kron(Id, PT.T)))
    rho_int = (-np.linalg.inv(L) @ rho0.reshape(-1)).reshape(dim, dim)
    return float(kS * np.real(np.trace(PS @ rho_int)))


partB = {}
for ratio in [0.1, 0.3, 1.0, 3.0, 10.0]:
    y0 = yield_unequal(0.0, 1.0, ratio)
    y90 = yield_unequal(np.pi / 2, 1.0, ratio)
    partB[ratio] = {"yield_0deg": y0, "yield_90deg": y90,
                    "anisotropy": float(y90 - y0)}
    print(f"  kT/kS={ratio:>5}: Phi_S(0)={y0:.4f} Phi_S(90)={y90:.4f} "
          f"anisotropy {y90 - y0:+.4f}")
eq = partB[1.0]
g2 = abs(eq["yield_0deg"] - 0.5075) < 0.005 and abs(eq["yield_90deg"] - 0.6719) < 0.005
print(f"G2 equal-rate Liouvillian matches original: {g2}")

fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
for n, d in partA.items():
    axes[0].plot(angles_deg, d["yields"], "o-", label=f"{n} nucleus")
axes[0].set_xlabel("field angle (deg)")
axes[0].set_ylabel("singlet yield")
axes[0].set_title("A: anisotropy vs nuclei count")
axes[0].legend()
ratios = sorted(partB, key=float)
axes[1].bar([str(r) for r in ratios],
            [partB[r]["anisotropy"] for r in ratios])
axes[1].set_xlabel("k_T / k_S")
axes[1].set_ylabel("anisotropy (90 - 0 deg)")
axes[1].set_title("B: reaction-rate asymmetry")
fig.tight_layout()
fig.savefig("radical_pair_deepened.png", dpi=150)

out = {"partA_nuclei": {str(k): v for k, v in partA.items()},
       "partB_rates": {str(k): v for k, v in partB.items()},
       "gates": {"G1_n1_reproduces": bool(g1),
                 "G2_liouvillian_matches": bool(g2)}}
with open("radical_pair_deepened_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values()) else "SOME GATE FAILED")
print("Saved radical_pair_deepened_results.json + PNG")
