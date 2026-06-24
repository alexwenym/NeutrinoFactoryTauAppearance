# Open questions

Unresolved physics & feasibility questions. Resolve → move the answer into the
relevant doc and note the decision in [04-decisions.md](04-decisions.md).

## Feasibility (make-or-break)
- ❓ Can off-axis muon **angular resolution** statistically separate τ-daughter muons
  from direct ν_μ CC muons? This is the central question the whole study must answer.
- ❓ How much do the two unseen neutrinos in τ→μν_μν_τ smear the decay angle — does
  enough information survive to be useful?

## Simulation geometry — RESOLVED 2026-06-18
- ✅ Observable: **both** transverse position x **and** arrival angle θ at the detector
  plane (the 2D distribution; they're correlated via production distance).
- ✅ Target region: τ's produced **all along the beam** (upstream rock + detector),
  bounded in practice by muon range.
- ✅ Comparison: with ν_μ→ν_τ oscillation vs without (Milestone 1 channel = ν_μ→ν_τ).
- ❓ τ production kinematics fidelity: collinear E_τ≈E_ν (v0) → inelasticity model
  (v1) → full generator (GENIE)? (still open)

## Modeling choices to settle
- ❓ Where do the τ's interact/decay — in the detector volume, or along a region of
  the baseline? Determines the production-point distribution for off-axis mapping.
- ❓ Detector assumptions: size, angular resolution, energy threshold for muons,
  on-/off-axis placement. What benchmark detector do we target?
- ❓ Which beam normalization (`accelerator` / `dump` / `baseline`) and E_μ (25 vs 50
  GeV) and baseline L do we anchor the headline result to?
- ❓ Do we need τ polarization in the decay spectrum, or is the unpolarized Michel
  spectrum sufficient at first pass?

## Validation
- ❓ Does our `oschelper`-based P(ν_e→ν_τ) agree with the legacy nuSQuIDS results
  over the relevant (E, L) range?
