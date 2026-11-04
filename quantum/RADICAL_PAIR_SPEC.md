# §XVI Radical-Pair Compass — Simulation Spec

*Second research line, from The Grand Unification of Being, Vol II §XVI
"The Avian Compass Theorem" (cryptochrome radical pairs). Spec written
2026-11-02; simulation build is the next work block after the Phoenix
post-selection run. Status: SPEC ONLY — not yet built or run.*

## The claim being modeled

In cryptochrome (the bird-retina protein), light creates a pair of
radicals — two molecules each holding one unpaired electron spin, born
together in the **singlet** state. The two spins precess under (a) Earth's
magnetic field (~50 µT) and (b) internal **hyperfine** fields from nearby
nuclei. Singlet and triplet configurations interconvert at a rate that
depends on the field's **direction** relative to the molecule, and only
the singlet channel recombines to the original chemical product. The
singlet yield Φ_S therefore varies with field angle θ — a chemical
compass. This is the leading hypothesis for avian magnetoreception
(Ritz et al., 2000) — a hypothesis, not settled fact (book Flag 7
discipline applies to anything published from this).

## The toy model (v1)

Three spins, the smallest system that shows the effect:

- **e1, e2** — the two electron spins (spin-1/2), initial state singlet
  (|↑↓⟩ − |↓↑⟩)/√2.
- **n1** — one nucleus (spin-1/2 stand-in for a nitrogen), hyperfine-coupled
  to e1 only, initially unpolarized (run |↑⟩ and |↓⟩, average).

Hamiltonian (angular-frequency units, ħ = 1):

    H(θ) = ω · (S1 + S2)  +  a⊥ (I_x S1x + I_y S1y) + a∥ I_z S1z

- Zeeman: field direction n̂(θ) = (sinθ, 0, cosθ) in the molecular frame;
  ω = γ_e B0, with B0 = 50 µT → ω/2π ≈ 1.4 MHz.
- Hyperfine: axial, a⊥/2π ≈ 0.5 MHz, a∥/2π ≈ 1.5 MHz (order-of-magnitude
  nitrogen values; v1 uses round numbers, sensitivity sweep later).

Observable: singlet probability over time, P_S(t) = ⟨S|ρ(t)|S⟩, and the
Haberkorn-style yield with equal recombination rates k_S = k_T = k:

    Φ_S(θ) = k ∫₀^∞ e^(−kt) P_S(t; θ) dt,   k ≈ 1/µs

Expected signature: Φ_S varies with θ (the anisotropy **is** the result);
toy-model anisotropy typically comes out cleaner/larger than in vivo —
report the number the model gives, labeled as the toy model it is.

## Qiskit implementation plan

- **Qubits:** 3. Singlet prep on (e1, e2): H + CNOT + phase (standard
  (|01⟩ − |10⟩)/√2 construction).
- **Evolution:** Trotterized exp(−iHt) (second-order Suzuki), t sampled
  to ~4 µs in scaled steps; cross-check against exact evolution via
  `PauliEvolutionGate` on statevector. Per-step terms: Zeeman rotations
  on e1, e2; hyperfine Pauli gadgets XX/YY/ZZ on (n1, e1).
- **Singlet readout:** apply CNOT(e1→e2) then H(e1); the singlet maps to
  |11⟩ on the electron pair (verified identity), so P_S = fraction of
  shots with electrons = 11, nucleus ignored.
- **Sweep:** θ ∈ {0°, 15°, …, 180°} → Φ_S(θ) table + chart (the chart is
  the publishable figure).
- **Verification gates (all must pass before any result is quoted):**
  1. No field, no hyperfine → P_S(t) = 1 at all t (state never leaves
     the singlet).
  2. Isotropic hyperfine (a⊥ = a∥) → Φ_S independent of θ (flat line).
  3. Trotter vs exact evolution agree within tolerance at the chosen dt.
  4. Anisotropy appears only when BOTH field and anisotropic hyperfine
     are present.
- **Runs where:** statevector/Aer locally or in Colab. Cost $0.

## Hardware feasibility (honest)

Full time evolution is deep: tens of gates per Trotter step × 50–200
steps → hundreds to thousands of gates, past what current devices return
cleanly for a percent-level yield difference. **v1 is simulation-only.**
A hardware version is a later stretch goal (few-step short-time dynamics
or heavy step reduction), decided only after the simulation exists and
its numbers are worth chasing.

## Extensions queued behind v1

- More nuclei (2–3), distinct k_S ≠ k_T kinetics, spin relaxation via
  Aer noise models, hyperfine parameter sensitivity sweep.
- A short write-up + chart tying the result back to §XVI (the chapter
  already exists — this gives it the same treatment §XX got: the
  author's own numbers behind the book's physics).
- Possible r/QuantumComputing post as the second research line, DOI'd
  via Zenodo with the teleportation repo pattern.

## Limits (state these wherever the result is shown)

One nucleus, spin-1/2, no molecular tumbling, no relaxation in v1.
This models the *mechanism* qualitatively; it is not a cryptochrome-
accurate prediction, and it does not prove birds use it — it shows the
physics the hypothesis rests on actually computes, in the author's own
simulation.
