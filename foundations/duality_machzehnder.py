# WHICH-PATH DUALITY (Codex index #8; Quantum Foundations, QSIM
# on Aer). A Mach-Zehnder interferometer on one qubit:
# H - phase(phi) - H. Unmarked: full fringes in phi.
# Mark the path: entangle the path qubit with a marker qubit
# (CNOT after the first H). Full marking erases the fringes --
# the information exists, whether or not anyone reads it.
# Partial marking (controlled-Ry(theta) instead of CNOT) trades
# visibility against distinguishability along V^2 + D^2 = 1,
# the Englert-Greenberger duality relation.
#   D = (1/2) || rho_marker(path0) - rho_marker(path1) ||_1
#     = sin(theta/2) for this marking; V = |cos(theta/2)|.
# Gates:
#   G1 unmarked fringe visibility >= 0.99 (Aer, exact statevector)
#   G2 fully marked visibility <= 0.01
#   G3 V^2 + D^2 = 1 within 0.01 across the theta sweep

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, Operator


def mz_probs(theta_mark, phis):
    """P(0) vs phase for a marking strength theta_mark (0 = none,
    pi = full CNOT-strength marking via CRY(pi) up to phases)."""
    out = []
    for phi in phis:
        qc = QuantumCircuit(2)
        qc.h(0)
        if theta_mark > 0:
            qc.cry(theta_mark, 0, 1)
        qc.p(phi, 0)
        qc.h(0)
        sv = Statevector.from_instruction(qc)
        probs = sv.probabilities([0])   # trace over marker
        out.append(float(probs[0]))
    return np.array(out)


def marker_D(theta_mark):
    """Trace distance between marker states conditioned on path."""
    # path 0 leaves marker |0>; path 1 rotates it by theta around Y
    rho0 = np.array([[1, 0], [0, 0]], dtype=complex)
    c, s = np.cos(theta_mark / 2), np.sin(theta_mark / 2)
    psi1 = np.array([c, s], dtype=complex)
    rho1 = np.outer(psi1, psi1.conj())
    return float(0.5 * np.sum(np.abs(np.linalg.eigvalsh(rho0 - rho1))))


phis = np.linspace(0, 2 * np.pi, 41)
rows = []
for th in [0.0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2,
           2 * np.pi / 3, np.pi]:
    P = mz_probs(th, phis)
    V = float((P.max() - P.min()) / (P.max() + P.min()))
    D = marker_D(th)
    rows.append({"theta": float(th), "V": V, "D": D,
                 "V2_plus_D2": float(V ** 2 + D ** 2)})
    print(f"theta={np.degrees(th):6.1f} deg: V={V:.4f} D={D:.4f} "
          f"V^2+D^2={V ** 2 + D ** 2:.4f}")

g1 = rows[0]["V"] >= 0.99
g2 = rows[-1]["V"] <= 0.01
g3 = all(abs(r["V2_plus_D2"] - 1) < 0.01 for r in rows)
print(f"G1 unmarked V>=0.99: {g1} | G2 marked V<=0.01: {g2} | "
      f"G3 duality holds: {g3}")

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for th, lab in [(0.0, "unmarked"), (np.pi / 2, "half-marked"),
                (np.pi, "fully marked")]:
    axes[0].plot(np.degrees(phis), mz_probs(th, phis), label=lab)
axes[0].set_xlabel("phase phi (deg)")
axes[0].set_ylabel("P(path qubit = 0)")
axes[0].set_title("Mach-Zehnder fringes vs path marking")
axes[0].legend(fontsize=8)
ths = [np.degrees(r["theta"]) for r in rows]
axes[1].plot(ths, [r["V"] for r in rows], "o-", label="V")
axes[1].plot(ths, [r["D"] for r in rows], "s-", label="D")
axes[1].plot(ths, [r["V2_plus_D2"] for r in rows], "^--",
             label="V^2 + D^2")
axes[1].set_xlabel("marking angle theta (deg)")
axes[1].set_title("Englert-Greenberger: V^2 + D^2 = 1")
axes[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig("duality_machzehnder.png", dpi=150)

out = {"rows": rows,
       "gates": {"G1_unmarked": bool(g1), "G2_marked": bool(g2),
                 "G3_duality": bool(g3)}}
with open("duality_machzehnder_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values()) else "SOME GATE FAILED")
print("Saved duality_machzehnder_results.json + PNG")
