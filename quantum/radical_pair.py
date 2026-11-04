# Radical-pair chemical compass -- toy model for Vol II Section XVI
# ("The Avian Compass Theorem") of The Grand Unification of Being.
#
# Model (see RADICAL_PAIR_SPEC.md): two electron spins born in a singlet
# + one spin-1/2 nucleus hyperfine-coupled to electron 1. Singlet yield
# Phi_S(theta) vs the angle theta between Earth's field and the molecular
# (hyperfine) axis. Exact dynamics by eigendecomposition (8x8 Hilbert
# space); Qiskit circuit evolution is used for verification gate 3.
#
# Units: microseconds and rad/us. B0 = 50 uT -> omega/2pi = 1.4 MHz
# (gamma_e/2pi ~ 28 GHz/T). Hyperfine: a_perp/2pi = 0.5 MHz,
# a_par/2pi = 1.5 MHz. Recombination k = 1 /us (Haberkorn, k_S = k_T).
#
# Verification gates (must pass before any number is quoted):
#   1. No field, no hyperfine -> P_S(t) = 1 at all times.
#   2. Isotropic hyperfine -> Phi_S independent of theta (flat).
#   3. Qiskit PauliEvolutionGate (exact and Trotter) reproduces the
#      numpy P_S at theta = 90 deg, t = 1 us.
#   4. Anisotropy requires BOTH field and anisotropic hyperfine:
#      field-only -> Phi_S = 1 for all theta; no-field -> flat.

import json

import numpy as np
from scipy.linalg import expm  # noqa: F401  (kept for ad-hoc checks)

# ---------- spin algebra (basis |e1> x |e2> x |n>) ----------
sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
I2 = np.eye(2, dtype=complex)


def op(a, b, c):
    return np.kron(np.kron(a, b), c)


S1 = [op(s, I2, I2) for s in (sx, sy, sz)]
S2 = [op(I2, s, I2) for s in (sx, sy, sz)]
IN = [op(I2, I2, s) for s in (sx, sy, sz)]

# singlet projector on the electron pair (x identity on nucleus)
singlet = np.zeros(4, dtype=complex)
singlet[1] = 1 / np.sqrt(2)   # |01>
singlet[2] = -1 / np.sqrt(2)  # |10>
PI_S = np.kron(np.outer(singlet, singlet.conj()), I2)

OMEGA = 2 * np.pi * 1.4    # rad/us, 50 uT
A_PERP = 2 * np.pi * 0.5   # rad/us
A_PAR = 2 * np.pi * 1.5    # rad/us
K_REC = 1.0                # /us


def hamiltonian(theta, omega=OMEGA, a_perp=A_PERP, a_par=A_PAR):
    n = np.array([np.sin(theta), 0.0, np.cos(theta)])
    H = sum(omega * n[j] * (S1[j] + S2[j]) for j in range(3))
    H = H + a_perp * (IN[0] @ S1[0] + IN[1] @ S1[1]) + a_par * (IN[2] @ S1[2])
    return H


def rho0(nucleus_up=True):
    psi = np.zeros(8, dtype=complex)
    # (e1,e2) singlet, nucleus in basis state: index = e1*4 + e2*2 + n
    n = 0 if nucleus_up else 1
    psi[0 * 4 + 1 * 2 + n] = 1 / np.sqrt(2)
    psi[1 * 4 + 0 * 2 + n] = -1 / np.sqrt(2)
    return np.outer(psi, psi.conj())


T_GRID = np.linspace(0, 12, 1201)  # us; e^-12 tail is negligible


def singlet_probabilities(theta, **hkw):
    """P_S(t) on T_GRID, averaged over the two nucleus polarizations."""
    H = hamiltonian(theta, **hkw)
    evals, evecs = np.linalg.eigh(H)
    dE = evals[:, None] - evals[None, :]          # E_i - E_j
    phase = np.exp(-1j * T_GRID[:, None, None] * dE[None])  # (t, i, j)
    out = np.zeros(len(T_GRID))
    for up in (True, False):
        r0 = evecs.conj().T @ rho0(up) @ evecs    # rho(0) in eigenbasis
        rt = (evecs @ (r0[None] * phase) @ evecs.conj().T)  # rho(t), lab frame
        out += np.real(np.einsum("tij,ji->t", rt, PI_S)) / 2
    return out


def yield_from_ps(ps):
    w = np.exp(-K_REC * T_GRID)
    return float(np.trapezoid(w * ps, T_GRID) * K_REC)


def phi_s(theta_deg, **hkw):
    return yield_from_ps(singlet_probabilities(np.radians(theta_deg), **hkw))


# ---------- verification gates ----------
print("=== Verification gates ===")
ps = singlet_probabilities(0.7, omega=0.0, a_perp=0.0, a_par=0.0)
g1 = float(np.max(np.abs(ps - 1.0)))
print(f"Gate 1 (no field, no hyperfine, P_S=1): max dev {g1:.2e} ->",
      "PASS" if g1 < 1e-9 else "FAIL")

iso = [phi_s(t, a_perp=2 * np.pi * 1.0, a_par=2 * np.pi * 1.0)
       for t in (0, 45, 90, 135, 180)]
g2 = max(iso) - min(iso)
print(f"Gate 2 (isotropic hyperfine -> flat): spread {g2:.2e} ->",
      "PASS" if g2 < 1e-9 else "FAIL")

field_only = [phi_s(t, a_perp=0.0, a_par=0.0) for t in (0, 90, 180)]
no_field = [phi_s(t, omega=0.0) for t in (0, 90, 180)]
g4a = max(field_only) - min(field_only)
g4b = max(no_field) - min(no_field)
print(f"Gate 4a (field only -> Phi_S=1 flat): values "
      f"{[round(v, 6) for v in field_only]} spread {g4a:.2e} ->",
      # absolute tolerance is integration-aware: Phi_S is a trapezoid
      # integral (dt=0.01), whose own error is ~1e-5; the value 1.000002
      # is the integrator, not the physics (P_S(t) itself is exact).
      "PASS" if g4a < 1e-9 and abs(field_only[0] - 1) < 1e-4 else "FAIL")
print(f"Gate 4b (no field -> flat): spread {g4b:.2e} ->",
      "PASS" if g4b < 1e-9 else "FAIL")

# Gate 3: Qiskit circuit evolution vs numpy, theta=90, t=1 us.
from qiskit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.synthesis import LieTrotter, SuzukiTrotter

th = np.pi / 2
om, ap, az = OMEGA, A_PERP, A_PAR
# Qiskit qubit order: qubit2=e1, qubit1=e2, qubit0=n (labels read q2 q1 q0).
H_op = SparsePauliOp.from_list([
    ("XII", om * np.sin(th) / 2),   # S1x  (e1 = leftmost char)
    ("IXI", om * np.sin(th) / 2),   # S2x
    ("ZII", om * np.cos(th) / 2),   # S1z
    ("IZI", om * np.cos(th) / 2),   # S2z
    ("XIX", ap / 4),                # I_x S1x  (n = rightmost char)
    ("YIY", ap / 4),                # I_y S1y
    ("ZIZ", az / 4),                # I_z S1z
]).simplify()


def prep_circuit(nucleus_up=True):
    qc = QuantumCircuit(3)
    e1, e2, n = 2, 1, 0
    if not nucleus_up:
        qc.x(n)
    qc.h(e1); qc.cx(e1, e2); qc.x(e2); qc.z(e1)  # singlet on (e1, e2)
    return qc


def ps_from_statevector(sv):
    amps = sv.data.reshape(2, 2, 2)  # [e1][e2][n]
    p = 0.0
    for n in range(2):
        pair = amps[:, :, n].reshape(4)
        p += abs(np.vdot(singlet, pair)) ** 2
    return float(p)


t_ev = 1.0
ps_num = float(singlet_probabilities(th)[np.argmin(np.abs(T_GRID - t_ev))])
vals = {}
for name, gate in [("exact", PauliEvolutionGate(H_op, time=t_ev)),
                   ("trotter", PauliEvolutionGate(
                       H_op, time=t_ev,
                       synthesis=LieTrotter(reps=100))),
                   ("suzuki2", PauliEvolutionGate(
                       H_op, time=t_ev,
                       synthesis=SuzukiTrotter(reps=50)))]:
    ps_acc = 0.0
    for up in (True, False):
        qc = prep_circuit(up)
        qc.append(gate, [0, 1, 2])
        ps_acc += ps_from_statevector(Statevector(qc)) / 2
    vals[name] = ps_acc
g3e = abs(vals["exact"] - ps_num)
g3t = abs(vals["trotter"] - ps_num)
# Informational: even a coarse Trotter setting (10 reps) reproduces the
# exact P_S at this system size -- Qiskit's synthesis evaluates this
# 3-qubit evolution essentially exactly, so Gate 3 validates the circuit
# Hamiltonian and qubit mapping, and says nothing about Trotter error
# at hardware scale (where step count vs depth remains the open cost).
qc10 = prep_circuit(True)
qc10.append(PauliEvolutionGate(H_op, time=t_ev,
                               synthesis=LieTrotter(reps=10)), [0, 1, 2])
qc10b = prep_circuit(False)
qc10b.append(PauliEvolutionGate(H_op, time=t_ev,
                                synthesis=LieTrotter(reps=10)), [0, 1, 2])
ps10 = (ps_from_statevector(Statevector(qc10))
        + ps_from_statevector(Statevector(qc10b))) / 2
print(f"Gate 3 (Qiskit vs numpy P_S at 90 deg, 1 us): numpy {ps_num:.6f}, "
      f"exact {vals['exact']:.6f} (dev {g3e:.2e}), "
      f"trotter100 {vals['trotter']:.6f} (dev {g3t:.2e}), "
      f"suzuki2x50 {vals['suzuki2']:.6f}; "
      f"trotter10 dev {abs(ps10 - ps_num):.2e} (also exact at 3 qubits) ->",
      "PASS" if g3e < 1e-6 and g3t < 5e-3 else "FAIL")

gates_pass = (g1 < 1e-9 and g2 < 1e-9 and g4a < 1e-9 and g4b < 1e-9
              and g3e < 1e-6 and g3t < 5e-3)
print("ALL GATES:", "PASS" if gates_pass else "FAIL -- do not quote results")

# ---------- the sweep ----------
thetas = list(range(0, 181, 15))
yields = [phi_s(t) for t in thetas]
print("\n=== Phi_S(theta) ===")
for t, y in zip(thetas, yields):
    print(f"  theta {t:>3} deg : {y:.4f}")
y = np.array(yields)
aniso = float(y.max() - y.min())
rel = aniso / float(y.mean())
print(f"Anisotropy: max-min = {aniso:.4f} "
      f"({100 * rel:.1f}% of mean yield {y.mean():.4f}); "
      f"max at {thetas[int(y.argmax())]} deg, min at {thetas[int(y.argmin())]} deg")

# sample P_S(t) curves for the record
curves = {str(t): singlet_probabilities(np.radians(t)).tolist()
          for t in (0, 90)}
with open("radical_pair_results.json", "w") as f:
    json.dump({"model": "3-spin toy: 2 electrons (singlet-born) + 1 nucleus",
               "B0_uT": 50, "omega_MHz": 1.4, "a_perp_MHz": 0.5,
               "a_par_MHz": 1.5, "k_rec_per_us": K_REC,
               "gates_pass": bool(gates_pass),
               "thetas_deg": thetas, "phi_s": yields,
               "anisotropy_abs": aniso, "anisotropy_rel": rel,
               "t_grid_us": T_GRID.tolist(), "ps_curves": curves},
              f, indent=2)
print("Saved radical_pair_results.json")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(thetas, yields, "o-", color="#0f4c5c", lw=2)
ax.set_xlabel("Field angle theta (degrees from molecular axis)")
ax.set_ylabel("Singlet yield Phi_S")
ax.set_title("Radical-pair compass: singlet yield vs field angle\n"
             "toy model - B0 50 uT, one nucleus (a_perp 0.5, a_par 1.5 MHz)")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("radical_pair_yield.png", dpi=160)
print("Saved radical_pair_yield.png")
