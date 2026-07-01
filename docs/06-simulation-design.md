# Simulation design

How we plan to simulate the measurement. Reflects the approach discussed
2026-06-18: simulate τ decays all along the beam from ν_τ interactions and find the
sensitivity, assuming the detector has sufficient angular resolution.

## Two milestones (channels)
- **Milestone 1 (technique demo): ν_μ → ν_τ.** The dominant tau-appearance channel
  (large, "atmospheric" amplitude). Establishes whether off-axis muons can reveal
  tau appearance at all. This is the current focus.
- **Milestone 2 (the README goal): ν_e → ν_τ + δ_CP.** Subdominant but carries the
  CP phase; the harder follow-on once the technique is shown.

## Observable (the money plot) — 2D
For muons arriving at the detector plane, the joint distribution of:
- **transverse position** x (offset from beam axis), and
- **arrival angle** θ relative to the beam axis.

They are correlated via the production distance d: a muon born distance d upstream,
emitted at angle θ, lands at x ≈ d·θ with arrival angle θ. Plotting both (and their
correlation) is more discriminating than either alone.

## The comparison (hypothesis test)
**With** ν_μ→ν_τ oscillation vs **without** it.
- The detector muon sample has two sources, both produced all along the beam:
  1. **Direct ν_μ CC** → muons ~along axis (only DIS p_T) = bright **on-axis core**.
  2. **τ→μ** (from appearance) → muons with the decay kick = **off-axis halo**.
- Turning ν_μ→ν_τ off: the core gets **brighter** (no ν_μ disappearance) and the
  **halo vanishes**. The signal is thus a core deficit + halo appearance — both
  captured by the 2D (x, θ) distribution.

## The Monte Carlo chain (Milestone 1: ν_μ→ν_τ)
1. **Sample ν_μ energy** from the beam flux (`flux.py`, `numu_flux`), ~collinear with
   the axis (beam angular spread ~ m_μ/E_μ ~ few mrad).
2. **Oscillation** P(ν_μ→ν_τ; E, L) (`osc.py`, initial flavor μ → nuind=1, τ
   component index 2). Also track the direct ν_μ survival for the core.
3. **ν_τ CC interaction** weight by σ(E) (`cross_sections.py`). Produce a τ with
   energy/angle. v0: τ collinear, E_τ ≈ E_ν. v1: add inelasticity y so
   E_τ = (1−y)E_ν and a τ production p_T.
4. **τ → μ ν_μ ν_τ**: apply BR (≈17.4%); sample the leptonic (Michel) muon spectrum
   in the τ rest frame; boost to lab → muon energy and angle relative to the τ/beam.
5. **Geometry + range acceptance**: pick the τ production point along the beam;
   accept the muon only if it can reach the detector (range cut, see below); project
   its direction → transverse position x and arrival angle θ at the detector plane.
6. **Histogram** the 2D (x, θ) distribution; also build the direct-ν_μ-CC core the
   same way (steps 3–5 with a prompt muon instead of a τ).
7. **Sensitivity**: χ²/likelihood across (x, θ) bins, with-osc vs no-ν_μ→ν_τ.

## v1 simplifications (explicit)
- Detector = an **idealized detection plane** (`detector.py`): the sim records truth-level
  θ, r, φ, E_μ for every muon crossing the plane. A downstream `apply_detector` hook adds
  an optional footprint + Gaussian θ/E smearing + threshold (pheno-level, **off by
  default** for now — "record all, cut later"). `sigma_theta` is meant to carry the TPC
  resolution **and** rock MCS in quadrature.
- **Multiple Coulomb scattering** now available (analytical Highland, `mcs=True`, off by
  default) — smears the arrival angle and landing position over the rock path; see the
  progress log. Full Moliere tails / range straggling remain a PROPOSAL-level upgrade.
- **Muon range IS kept** even in v1 — it is what bounds "all along the beam." A muon
  is accepted only if its production point is within range of the detector
  (~19 m at 10 GeV, ~87 m at 50 GeV in rock). Without it the off-axis halo gets
  unphysical contributions from arbitrarily far upstream. A simple range formula
  (CSDA / dE/dx) suffices; no MCS needed yet.

## Key physics caveat — why "enough angular resolution" is an assumption
The signal is the muon's transverse kick from the τ decay. Competing effects:
- **Muon range** limits how far upstream a τ can be and still deliver a muon to the
  detector: ~19 m in rock at 10 GeV, ~87 m at 50 GeV. So "along the beam" ≈ the last
  few tens of meters of rock before the detector (the lever arm).
- **Multiple Coulomb scattering** smears the muon angle in rock.

Back-of-envelope (to be checked properly):
- τ decay opening angle ~ m_τ/E_τ ~ 0.1–0.2 rad.
- MCS over ~20 m rock at ~10 GeV ~ 0.02 rad.
- Decay angle beats MCS by ~×few, and the ratio is ~**energy-independent** (both
  ∝ 1/E). Encouraging, not doomed — but MCS is THE check that validates the
  "enough angular resolution" assumption. First realism item after v1.

Note on lever arm: decay-angle displacement ∝ d (lever arm), MCS displacement ∝
d^{3/2}; MCS *angle* ∝ √d. So shorter production-to-detector distance helps the
angular measurement — the optimum production region is the near-detector rock.

## Reusable vs to-build
- Reusable: `flux.py` (angle-resolved), `osc.py`/`oschelper` (P, matter, δ_CP),
  `cross_sections.py` (σ_τ).
- To build: τ→μ decay kinematics module; geometry/projection; δ_CP sensitivity
  layer. Optional: ν_τ CC differential kinematics (inelasticity); ν_μ CC background
  sample (reuses `numu_flux`); muon transport (range + MCS) for the realism pass.

## Open modeling choices
See [05-open-questions.md](05-open-questions.md) — esp. exact observable definition
(arrival angle vs transverse offset), whether the target is detector-volume-only or
upstream rock + detector, and τ production kinematics fidelity.
