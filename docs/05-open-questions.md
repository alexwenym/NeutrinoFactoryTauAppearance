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
- ✅ τ production kinematics — RESOLVED 2026-06-24: full DIS production from **GENIE
  d²σ/dx dy tables** (`resources/diffxsec`), same source for the ν_μ core and σ. The real
  τ production angle is ~93 mrad (τ-mass Q² floor), not ultra-forward. ❓ remaining:
  QE/RES below the DIS region (unmodeled, low weight at the ~12–24 GeV bulk).

## Modeling choices to settle
- ❓ Where do the τ's interact/decay — in the detector volume, or along a region of
  the baseline? Determines the production-point distribution for off-axis mapping.
- ❓ Detector assumptions: size, angular resolution, energy threshold for muons,
  on-/off-axis placement. What benchmark detector do we target?
- ❓ Which beam normalization (`accelerator` / `dump` / `baseline`) and E_μ (25 vs 50
  GeV) and baseline L do we anchor the headline result to?
- ✅ τ polarization — RESOLVED 2026-06-24: `tau_decay.py` uses alouette with **fixed
  helicity** (−1 τ⁻ / +1 τ⁺, V−A). It shifts muons forward (lab ⟨E_μ⟩/E_τ 0.354→0.405,
  median θ 59→53 mrad) so it is kept on by default. ❓ remaining: **energy-dependent**
  CC polarization (Hagiwara et al.) instead of the fixed-helicity approximation.

## Validation
- ❓ Does our `oschelper`-based P(ν_e→ν_τ) agree with the legacy nuSQuIDS results
  over the relevant (E, L) range?
