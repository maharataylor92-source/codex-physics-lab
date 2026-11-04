# CODEX PHYSICS LAB — INDEX

*Created 2026-11-04 at Mahara's direction. Scope (his words, 2026-11-04): the
Codex is a physics lab menu — any physics formula in the Codex / Formulas
files may seed an experiment of any type (classical simulation, numerical
experiment, data analysis, quantum hardware where it fits). NOT limited to
quantum teleportation or IBM hardware. Every entry is labeled honestly by
method: a simulation result is a simulation result, never a measurement.*

## Method legend
- **SIM** — classical numerical simulation (solver integrates the physics)
- **NUM** — numerical/scale computation (exact formulas, real constants)
- **DATA** — analysis of already-collected real data
- **QPU** — quantum hardware experiment (IBM, when account is active)
- **QSIM** — quantum circuit simulator (Aer)
- **—** — no honest experiment at our scale; stays a philosophical anchor

## The index

| # | Codex source | Field | Formula | Experiment | Method | Status |
|---|---|---|---|---|---|---|
| 1 | Vol I | General Relativity | G_μν = (8πG/c⁴)T_μν | Schwarzschild radius + gravitational time-dilation scales for real masses | NUM | open (quick win) |
| 2 | Vol I | Electromagnetism | ∇·B = 0 | Dipole field solver: field lines close, numerical ∇·B ≈ 0 on a grid (no monopoles) | SIM | open |
| 3 | Vol I | Statistical Physics | S = k·log W | Entropy ledger of the Bell pair + all saved count distributions | DATA | **DONE 2026-11-04** — S_AB 0.155 bits, I(A:B) 1.845/2, E_F 0.935 bits (teleportation repo, `notebooks/entropy_analysis.py`, v1.3) |
| 4 | Vol I | Special Relativity | E = mc² | Mass–energy scale table + relativistic vs Newtonian kinetic energy vs v/c | NUM | open (quick win) |
| 5 | Vol I | Celestial Mechanics | v = √(GM/r) | Orbit simulator: circular + elliptical orbits, Kepler III check, escape velocity, Earth/Sun + ISS real-scale cases | SIM | **DONE 2026-11-04** — all gates PASS: circular radius drift 2e-13, Kepler ratio exact, ellipse period matches, Earth 29.78 km/s / 365.26 d, ISS 7.67 km/s / 92.4 min (`classical/orbit_simulator.py`) |
| 6 | Vol I | Thermodynamics | q = −k∇T | Heat-flow simulation: rod with held ends, transient vs analytic series, steady-state flux | SIM | **DONE 2026-11-04** — steady state exact-linear, center T matches Fourier series (41.16 °C both), center hits 45 °C: copper 22 s, aluminium 26 s, iron 114 s (`classical/heat_flow.py`) |
| 7 | Vol I | Thermal Energy | ΔQ = mcΔT | Energy-balance gate inside the heat sim: stored heat (mcΔT integral) vs time-integrated boundary flux | SIM | **DONE 2026-11-04** — balance within 1.9% (with #6) |
| 8 | Vol II §XI | Quantum Information | \|Ψ⟩ = α\|0⟩ + β\|1⟩ | The substrate of every quantum run; interference/which-path duality is the open experiment on it | QPU/QSIM | partial (in all quantum work) |
| 9 | Vol II §XII | Quantum Matrix | S_bh = kA/(4ℓ_p²) | Black-hole entropy scale computation only — not a test | NUM | scale-only, do not oversell |
| 10 | Vol II §XVIII | Quantum Limits | ΔE·Δt ≥ ħ/2 | Mandelstam–Tamm quantum speed limit: segmented π rotation vs the bound | QSIM→QPU | **BUILT, Aer PASS** (`qsl_phoenix.py`); hardware armed, blocked on IBM verification |
| 11 | Vol II | Quantum Gravity | Ĥ\|Ψ⟩ = 0 (Wheeler–DeWitt) | None — frontier equation, no bench version | — | philosophy anchor |
| 12 | Vol II §XX | Quantum Teleportation | \|Φ⁺⟩ = (\|00⟩+\|11⟩)/√2 | Six-state teleportation benchmark, 3 IBM devices + Phoenix follow-up arc | QPU | **DONE** — repo v1.2, DOI 10.5281/zenodo.23934913 |
| 13 | Vol II (String Field) + Matrix | Waves | ∂²ψ/∂t² = v² ∂²ψ/∂x² | Double-source interference: fringe pattern, node lines (P₁+P₂=0), fringe spacing vs λL/d | SIM | **DONE 2026-11-04** — all gates PASS after a v1–v2 post-mortem (transposed source indexing + cavity ringing): visibility 0.98, nodes at 0.8% of maxima, spacing 0.370 m vs exact-geometry 0.349 m (paraxial λL/d = 0.320 m runs ~9% narrow at sinθ=0.4, as theory says) (`classical/wave_interference.py`) |
| 14 | Vol III | Fluid Dynamics | Navier–Stokes | Simplified 2D flow past an obstacle (lattice-Boltzmann or projection solver); vortex shedding vs viscosity | SIM | open (heavier build) |
| 15 | Vol III + §XVI | Quantum Biology | τ = m·B·sin(θ) | Radical-pair singlet yield vs field angle and Earth-range field strength (25–65 µT) | SIM | **DONE** — original + 2026-11-04 field sweep: anisotropy peaks ~45 µT (33%) (`quantum/radical_pair.py`, `quantum/radical_pair_fieldsweep.py`) |
| 16 | Vol III (torque's quantum face) | Spin Physics | Larmor precession | Ramsey interferometry on ibm_phoenix: precession frequency + T2* | QSIM→QPU | **BUILT, Aer PASS** (`ramsey_phoenix.py`); hardware armed, blocked on IBM verification |
| 17 | Vol III | Parallel Computing | Amdahl S(p) | Speedup curves vs parallel fraction — exact computation | NUM | open (quick win) |
| 18 | Vol II §XIX | Quantum Mechanics | Tunneling T ≈ e^(−2κw) | Wavepacket vs rectangular barrier: transmission vs width/height/energy, packet vs exact formula | SIM | **DONE 2026-11-04** — w=1/2/3 → 19.9%/1.5%/0.13% (`quantum/tunneling.py`) |
| 19 | Vol II §XXV | Quantum Mechanics | ρ_ij(t) decoherence | Decoherence-rate fit to the Phoenix hold sweep: precession + fringe-collapse envelope | DATA | **DONE 2026-11-04** — rotation required; ≈π in ~17–20 µs; T2* ~10–25 µs (teleportation repo, `notebooks/decoherence_fit.py`, v1.3) |
| 20 | Repository §2 | Quantum Limits | Δx·Δp ≥ ħ/2 | Gaussian wavepacket Fourier tradeoff: narrow in x ⇒ broad in p, product → ħ/2 | SIM | open |
| 21 | Vol II §XVII | Quantum Mechanics | E₀ = ½ħω | Harmonic-oscillator spectrum / zero-point scale computation | NUM/SIM | open |
| 22 | Vol III | System Dynamics | Ω₀ = 1.3689910 | None — no derivation exists; philosophical anchor by standing accuracy position | — | philosophy anchor |

## Rules of the lab
1. Every result states its method (SIM / NUM / DATA / QPU / QSIM) wherever it is quoted.
2. Simulation gates before headlines: conservation laws, analytic limits, and known special cases must pass before any number is reported.
3. Philosophy anchors (#11, #22, and the interpretive layers everywhere) are never presented as experimentally validated.

**Published 2026-11-04:** the six simulation experiments live in the public repo github.com/maharataylor92-source/codex-physics-lab, archived as Zenodo DOI 10.5281/zenodo.23963471 (concept 10.5281/zenodo.23963470). The entropy ledger (#15) and decoherence fit (#16) were published the same day in the teleportation repo v1.3 (DOI 10.5281/zenodo.23963444; concept 10.5281/zenodo.23809803).
