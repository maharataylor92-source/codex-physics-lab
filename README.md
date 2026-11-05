# Codex Physics Lab

**Author:** Mahara Taylor (ManyMoons), independent researcher

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23964708.svg)](https://doi.org/10.5281/zenodo.23964708)

## How to cite

Taylor, M. (2026). *Codex Physics Lab: gate-checked experiments
from the Grand Unification of Being formula codex* (v1.1.0).
Zenodo. https://doi.org/10.5281/zenodo.23964708

The concept DOI https://doi.org/10.5281/zenodo.23963470 always
resolves to the latest version of this record.

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

## Experiments (v1.0.0 + the November 4, 2026 deepening batches)

### quantum/

| Experiment | Field | Method | Headline result |
|---|---|---|---|
| `tunneling.py` | Quantum mechanics | SIM | Wavepacket vs a classically forbidden barrier (E=0.5, V=1.5): transmission 19.9% / 1.5% / 0.13% at widths 1 / 2 / 3, tracking the exact analytic formula (19.2% / 1.2% / 0.07%). All gates pass. |
| `radical_pair.py` + `RADICAL_PAIR_SPEC.md` | Quantum biology / spin chemistry | SIM | Radical-pair chemical compass toy model (2 electrons + 1 nucleus, exact dynamics, Qiskit cross-checked): singlet yield Φ_S swings 0.5075 → 0.6719 with field angle at 50 µT; anisotropy 0.174 (≈32% of mean). Toy model — one nucleus, k_S = k_T. |
| `radical_pair_fieldsweep.py` | Quantum biology / spin chemistry | SIM | The same model across Earth's real field range (25–65 µT) × angle: anisotropy peaks at ≈45 µT (0.177, 33% of mean) and fades at both ends (25 µT: 24%). Sanity gate reproduces the 50 µT run exactly. |
| `tunneling_double_barrier.py` | Quantum mechanics | SIM | Resonant tunneling: exact T = 0.9955 at the E = 0.471 resonance of a symmetric double barrier; wavepacket contrast 8.4× on/off resonance (absolute packet T is small — a narrow resonance accepts only a thin spectral slice, stated in the file). |
| `alpha_decay.py` | Nuclear physics | NUM vs DATA | Gamow's 1928 theory with NO fitted parameters (r₀ = 1.2 fm) against 10 measured half-lives: correlation 0.9987, predicted span 24.8 vs measured 24.2 decades, RMS residual 0.61 decades, Geiger–Nuttall slope 161 vs 157. Tunneling, checked against nature. |
| `radical_pair_deepened.py` | Quantum biology / spin chemistry | SIM | Robustness of the compass: anisotropy washes out with nuclei count (0.174 → 0.105 → 0.055 for n = 1→2→3); reaction-rate asymmetry k_T/k_S = 0.1 erases and slightly inverts it (−0.011), while 3.0 amplifies it (0.255). |

### classical/

| Experiment | Field | Method | Headline result |
|---|---|---|---|
| `orbit_simulator.py` | Celestial mechanics | SIM | Newtonian RK4 orbits from v = √(GM/r): circular-orbit radius drift 2×10⁻¹³ over 5 periods; Kepler T ∝ r^3/2 exact; ellipse period matches Kepler; escape at 1.05 v_esc. Real scale: Earth 29.78 km/s, 365.26 d; ISS 7.67 km/s, 92.4 min. |
| `heat_flow.py` | Thermodynamics | SIM | Fourier heat flow in a rod (q = −k∇T): steady state exactly linear; transient center temperature matches the analytic Fourier series to 0.01 °C; ΔQ = mcΔT energy balance vs integrated flux within 1.9%. Center reaches 45 °C: copper 22 s, aluminium 26 s, iron 114 s. |
| `wave_interference.py` | Wave physics / optics | SIM | 2D wave-equation double-source interference: visibility 0.98, node lines at 0.8% of peak (P₁ + P₂ = 0), fringe spacing 0.370 m vs exact geometry 0.349 m (the paraxial λL/d = 0.320 m runs ~9% narrow at sin θ = 0.4, as theory says). The file's header keeps the v1–v2 post-mortem: a transposed-indexing bug and cavity ringing that the gates caught. |
| `orbit_transfer.py` | Celestial mechanics | SIM | Hohmann transfer r = 1 → 4 flown by the RK4 integrator: arrival radius exact, circularized to 0.021% over a full period, transfer on the analytic schedule. Real scale: Earth → Mars in 259 days, departure Δv 2.95 km/s. |
| `heat_flow_extended.py` | Thermodynamics | SIM | 2D plate (matches the 1D linear profile to 2.7×10⁻⁸ °C) plus a layered insulation stack: R = 1.583 m²K/W, the copper layer drops 0.01 °C (a spreader, not an insulator), and a plain fiberglass wall of equal thickness beats the stack — reported as found. |
| `wave_whichpath.py` | Wave physics / optics | SIM | The which-path analog: attenuate source B by β and the visibility follows V = 2β/(1+β²) within 0.01 at every β; fully blocked, V = 0.045; V² + D² = 1 within 0.011. Two measurement post-mortems (estimator artifacts, edge reflections) are documented in the file. |
| `divergence_b.py` | Electromagnetism | SIM | Maxwell's no-monopole law, checked three ways on a point dipole: net flux through any closed box is zero to machine precision (~10⁻¹⁶ of gross); a traced field line rides the exact curve r = L sin²θ to 0.02%; numerical ∇·B is 0.18% of the field's own scale (median) and shrinks with textbook Δ² convergence — the first global-max gate failed on core-boundary truncation and the fix is documented in the file. |
| `navier_stokes.py` | Fluid dynamics | SIM | Lattice-Boltzmann (D2Q9) channel flow past a cylinder: the wake is dead steady at Re = 20 (probe RMS 7×10⁻⁷) and the von Kármán vortex street emerges on its own at Re = 80 (spectral prominence ~21,000). St = 0.288 on inlet speed — high vs the free-stream ~0.16 because the 33% channel blockage speeds the flow past the cylinder; stated plainly in the file. |
| `double_pendulum.py` | Classical mechanics / chaos | SIM | Chaos measured, not asserted: two double pendulums released 10⁻⁹ rad apart decorrelate in ~16 s; the Benettin Lyapunov exponent is L = 1.51 /s, while a small-oscillation control pair gives L ≈ 0.03 /s — the integrator does not manufacture the chaos it reports. |
| `membrane_modes.py` | Wave physics | SIM | Chladni patterns from first principles: the discrete Laplacian's lowest eigenvalues match the analytic ladder π²(m²+n²) to 0.19%, the (1,2)/(2,1) degeneracy holds, and a (2,1) mode excited in the time-domain solver rings at the predicted frequency. (A probe placed on a node failed the dynamics gate first — post-mortem in the file.) |
| `mercury_perihelion.py` | General relativity | SIM | The 43-arcsecond anomaly, reproduced: integrating the GR orbit equation over 200 Mercury orbits gives 43.05 arcsec/century of perihelion advance vs the observed 42.98, while the Newtonian control run drifts just 0.065. (The first pass read 57.7 with a +5.6 control drift — RK4 phase error on a coarse grid; the fix is documented in the file.) |
| `feigenbaum.py` | Chaos theory | NUM | Feigenbaum's δ from the logistic map's superstable points (located by exact finite iteration — period classification failed twice at the 10⁻⁴ level, post-mortem in the file): ratios 4.709, 4.681, 4.663, 4.668 → δ = 4.669202. |

### foundations/

| Experiment | Field | Method | Headline result |
|---|---|---|---|
| `uncertainty_wavepacket.py` | Quantum foundations | SIM | Free Gaussian packet: Δx·Δp starts at exactly ħ/2 (0.5000) and only rises; Δx(t) tracks the analytic spreading curve to 0.00%; Δp constant. The uncertainty floor, watched live. |
| `duality_machzehnder.py` | Quantum foundations | QSIM | Mach–Zehnder with a marker qubit: visibility falls 1 → 0 as path marking strengthens; V² + D² = 1.0000 at all seven marking angles (Englert–Greenberger), exact statevector on Aer. |
| `zero_point.py` | Quantum mechanics | NUM | Zero-point scale table: H₂ E₀ = 0.2728 eV (known 0.27); every molecular vibration's spacing exceeds 5 kT at 300 K (frozen in the ground state); a 1 Hz pendulum's E₀ is ~10⁻¹⁴ kT by contrast. |
| `amdahl.py` | Computing | NUM | Amdahl's law ceilings verified: S∞ = 1/(1−p); reaching 90% of the ceiling takes s = 171 processors at p = 0.95 and s = 891 at p = 0.99. The serial fraction is the wall. |
| `gr_scales.py` | General relativity | NUM | Schwarzschild radii from 1 kg to M87*, and gravitational time dilation where it is measured daily: the computed GPS gain is +45.7 µs/day, exactly the canonical figure (net +38.5 with the SR term); a neutron-star surface clock ticks at 0.766×. |
| `emc2.py` | Special relativity | NUM | Rest-energy scales: 1 kg = 21.5 Mt TNT equivalent; the Sun converts 4.26×10⁹ kg/s; U-235 fission releases 0.091% of the mass. Relativistic vs Newtonian KE: Newton is low by 0.76% at 0.1c and 1,143% at 0.99c. (A wrong gate anchor was caught and corrected mid-build — kept in the file.) |
| `ising_model.py` | Statistical physics | SIM | Spontaneous magnetization by Metropolis Monte Carlo: the susceptibility peak lands at T = 2.313 vs Onsager's exact Tc = 2.269 (1.9% — the expected finite-size shift for a 48×48 lattice); ordered at T = 1.6 (⟨|m|⟩ = 0.94), disordered at 3.2 (0.05). The first run's peak sat 4.1% high from critical slowing down; the 6× longer near-Tc sampling is documented in the file. |
| `percolation.py` | Statistical physics | SIM | The spanning transition: site percolation on a 128×128 lattice switches from islands to a continent at p_c = 0.5909, vs the known 0.5927. |
| `planck_spectrum.py` | Thermal radiation | NUM | Two laws inside one curve: numerically integrating Planck's spectrum recovers σ to 0.000% and Wien's constant to 0.001% at four temperatures, and M(3000)/M(1000) = 81.00 exactly as T⁴ demands. |

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
