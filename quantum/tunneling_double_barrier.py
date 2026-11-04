# TUNNELING DEEPENED, part 1 -- the DOUBLE barrier and resonant
# tunneling (extends tunneling.py, Codex section XIX).
#
# Two identical rectangular barriers (height V, width w) separated
# by a well of width L. A particle whose energy matches a quasi-
# bound level of the well transmits (nearly) PERFECTLY, even though
# each barrier alone would reflect almost all of it -- the cir-
# cuit QED / resonant-tunneling-diode effect, in 1D.
#
# Part 1 -- exact T(E) by transfer matrix (same propagator as the
# single barrier, applied segment by segment). Symmetric double
# barriers have T = 1 exactly at resonance.
# Part 2 -- split-step wavepackets tuned ON and OFF resonance.
# Units hbar = m = 1, as in tunneling.py.
#
# Gates:
#   G1 the resonance nearest E ~ 0.8 has exact T >= 0.98
#   G2 far off resonance exact T is tiny (< 0.02)
#   G3 wavepacket T(on resonance) >= 4x T(off resonance)

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V, w, L = 1.5, 1.5, 2.0


def T_exact_double(E):
    """Transfer matrix across barrier-well-barrier.  k in free
    regions, q = sqrt(2(E-V)) (imaginary below V) in barriers,
    k again in the well."""
    k = np.sqrt(2 * E + 0j)
    q = np.sqrt(2 * (E - V) + 0j)

    def interface(k1, k2):
        t = 2 * k1 / (k1 + k2)
        r = (k1 - k2) / (k1 + k2)
        return np.array([[1, r], [r, 1]], dtype=complex) / t

    def prop(kk, d):
        return np.diag([np.exp(1j * kk * d), np.exp(-1j * kk * d)])

    Mtot = interface(q, k)                     # exit face, barrier 2
    Mtot = prop(q, w) @ Mtot
    Mtot = interface(k, q) @ Mtot
    Mtot = prop(k, L) @ Mtot
    Mtot = interface(q, k) @ Mtot
    Mtot = prop(q, w) @ Mtot
    Mtot = interface(k, q) @ Mtot
    t_amp = 1.0 / Mtot[0, 0]
    return float(abs(t_amp) ** 2)


Es = np.linspace(0.05, 1.45, 1401)
Ts = np.array([T_exact_double(e) for e in Es])
# resonances = local maxima above 0.5
res = []
for i in range(1, len(Es) - 1):
    if Ts[i] > 0.5 and Ts[i] >= Ts[i - 1] and Ts[i] >= Ts[i + 1]:
        res.append((float(Es[i]), float(Ts[i])))
E_res, T_res = max(res, key=lambda p: abs(p[0] - 0.8) * -1 if False else p[1])
E_res, T_res = min(res, key=lambda p: abs(p[0] - 0.8))
E_off = E_res + 0.35
T_off = T_exact_double(E_off)
print(f"resonances (E, T): {[(round(a,3), round(b,3)) for a,b in res]}")
print(f"chosen resonance E={E_res:.3f} T={T_res:.4f}; "
      f"off-resonance E={E_off:.3f} T={T_off:.5f}")
g1 = T_res >= 0.98
g2 = T_off < 0.02
print(f"G1 resonance T>=0.98: {g1} | G2 off-res T<0.02: {g2}")


# ---- wavepackets, split-step (machinery from tunneling.py) ----
def packet_T(E0, label):
    N, Lx = 16384, 200.0
    x = np.linspace(-Lx / 2, Lx / 2, N, endpoint=False)
    dx = x[1] - x[0]
    sigma = 8.0
    x0 = -34.0
    k0 = np.sqrt(2 * E0)
    psi = (np.exp(-(x - x0) ** 2 / (4 * sigma ** 2))
           * np.exp(1j * k0 * x)).astype(complex)
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * dx)
    pot = np.zeros(N)
    b1 = (x >= 0) & (x < w)
    well = (x >= w) & (x < w + L)
    b2 = (x >= w + L) & (x < 2 * w + L)
    pot[b1 | b2] = V
    kk = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    dt = 0.02
    steps = int(80.0 / dt)
    half = np.exp(-0.5j * pot * dt)
    kin = np.exp(-0.5j * kk ** 2 * dt)
    for _ in range(steps):
        psi *= half
        psi = np.fft.ifft(kin * np.fft.fft(psi))
        psi *= half
    dens = np.abs(psi) ** 2 * dx
    xedge = 2 * w + L
    T = float(np.sum(dens[x > xedge + 1]))
    R = float(np.sum(dens[x < -1]))
    print(f"  packet {label}: E={E0:.3f}  T={T:.4f} R={R:.4f}")
    return T


Tp_on = packet_T(E_res, "ON resonance")
Tp_off = packet_T(E_off, "OFF resonance")
g3 = Tp_on >= 4 * max(Tp_off, 1e-6)
print(f"G3 packet T_on >= 4x T_off: {g3} ({Tp_on:.4f} vs {Tp_off:.4f})")

fig, ax = plt.subplots(figsize=(7.2, 3.6))
ax.plot(Es, Ts, lw=1.2)
ax.axvline(E_res, color="g", ls="--", lw=0.9)
ax.axvline(E_off, color="r", ls="--", lw=0.9)
ax.set_xlabel("energy E (hbar=m=1)")
ax.set_ylabel("exact transmission T(E)")
ax.set_title(f"Double barrier (V={V}, w={w}, gap={L}): resonant tunneling")
fig.tight_layout()
fig.savefig("tunneling_double_barrier.png", dpi=150)

out = {"barrier": {"V": V, "width": w, "gap": L},
       "resonances": [{"E": a, "T": b} for a, b in res],
       "chosen": {"E_res": E_res, "T_exact_res": T_res,
                  "E_off": E_off, "T_exact_off": T_off,
                  "T_packet_on": Tp_on, "T_packet_off": Tp_off},
       "gates": {"G1_resonance_unity": bool(g1),
                 "G2_off_resonance_tiny": bool(g2),
                 "G3_packet_contrast": bool(g3)}}
with open("tunneling_double_barrier_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("ALL GATES PASS" if all(out["gates"].values()) else "SOME GATE FAILED")
print("Saved tunneling_double_barrier_results.json + PNG")
