# Codex Physics Lab

**Author:** Mahara Taylor (ManyMoons), independent researcher

A working laboratory built from the formula codex of *The Grand
Unification of Being*: every experiment here starts from a formula in
that codex and turns it into a runnable, gate-checked computation.
The lab's scope is deliberately wide — any physics in the codex, by
whatever method honestly fits. It is **not** a quantum-teleportation
project (that program lives, with its IBM hardware data and DOI, in
[quantum-teleportation-benchmark](https://github.com/maharataylor92-source/quantum-teleportation-benchmark)).

## Method labels — read these first

Every result in this repo carries its method. A simulation result is
a simulation result, never a measurement.

- **SIM** — classical numerical simulation (a solver integrates the physics)
- **NUM** — numerical/scale computation from exact formulas
- **DATA** — analysis of real, previously collected data
- **QSIM / QPU** — quantum circuit simulator / quantum hardware

## Lab rules

1. Every result states its method wherever it is quoted.
2. Gates before headlines: conservation laws, analytic limits, and
   known special cases must pass before any number is reported.
   Failed gates and their post-mortems stay in the record.
3. The codex's philosophical anchors (e.g. Wheeler–DeWitt, Ω₀) are
   never presented as experimentally validated. They have no entry
   here, because they have no honest experiment.

## Experiments (v1.0.0, November 4, 2026)

### quantum/

| Experiment | Field | Method | Headline result |
|---|---|---|---|
| `tunneling.py` | Quantum mechanics | SIM | Wavepacket vs a classically forbidden barrier (E=0.5, V=1.5): transmission 19.9% / 1.5% / 0.13% at widths 1 / 2 / 3, tracking the exact analytic formula (19.2% / 1.2% / 0.07%). All gates pass. |
| `radical_pair.py` + `RADICAL_PAIR_SPEC.md` | Quantum biology / spin chemistry | SIM | Radical-pair chemical compass toy model (2 electrons + 1 nucleus, exact dynamics, Qiskit cross-checked): singlet yield Φ_S swings 0.5075 → 0.6719 with field angle at 50 µT; anisotropy 0.174 (≈32% of mean). Toy model — one nucleus, k_S = k_T. |
| `radical_pair_fieldsweep.py` | Quantum biology / spin chemistry | SIM | The same model across Earth's real field range (25–65 µT) × angle: anisotropy peaks at ≈45 µT (0.177, 33% of mean) and fades at both ends (25 µT: 24%). Sanity gate reproduces the 50 µT run exactly. |

### classical/

| Experiment | Field | Method | Headline result |
|---|---|---|---|
| `orbit_simulator.py` | Celestial mechanics | SIM | Newtonian RK4 orbits from v = √(GM/r): circular-orbit radius drift 2×10⁻¹³ over 5 periods; Kepler T ∝ r^3/2 exact; ellipse period matches Kepler; escape at 1.05 v_esc. Real scale: Earth 29.78 km/s, 365.26 d; ISS 7.67 km/s, 92.4 min. |
| `heat_flow.py` | Thermodynamics | SIM | Fourier heat flow in a rod (q = −k∇T): steady state exactly linear; transient center temperature matches the analytic Fourier series to 0.01 °C; ΔQ = mcΔT energy balance vs integrated flux within 1.9%. Center reaches 45 °C: copper 22 s, aluminium 26 s, iron 114 s. |
| `wave_interference.py` | Wave physics / optics | SIM | 2D wave-equation double-source interference: visibility 0.98, node lines at 0.8% of peak (P₁ + P₂ = 0), fringe spacing 0.370 m vs exact geometry 0.349 m (the paraxial λL/d = 0.320 m runs ~9% narrow at sin θ = 0.4, as theory says). The file's header keeps the v1–v2 post-mortem: a transposed-indexing bug and cavity ringing that the gates caught. |

## The index

`CODEX_PHYSICS_LAB_INDEX.md` maps **every** formula in the codex to
an experiment, a method, and a status — including the entries that
live in the teleportation repo (entropy ledger, decoherence fit),
the ones armed for hardware (Ramsey precession, quantum speed
limit), and the ones honestly marked as having no experiment.

## Running

Python 3 with numpy / scipy / matplotlib (`radical_pair.py` also
uses Qiskit for one verification gate). Each script is
self-contained, prints its verification gates, and writes its
`*_results.json` and chart next to itself.

## Licenses

- **Code:** MIT License — see `LICENSE`.
- **Data and documentation:** CC BY 4.0 (Mahara Taylor).

*Independent research. The Grand Unification of Being is the
author's book project; this lab is its experimental companion —
formulas are quoted from its codex, results stand or fall on
their gates.*
