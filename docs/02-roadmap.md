# Roadmap — gaps & planned work

Legend: ✅ done · 🚧 in progress · ⬜ not started · 🐛 bug

## Environment
- Cluster `.venv` (**py3.8**) with numpy + matplotlib + **alouette 1.0.1**; `oschelper`
  rebuilt into it with `LDSHARED="g++ -shared" CC=g++ pip install --no-build-isolation
  ./oschelper` (a plain build links without libstdc++ and import fails). Run code via
  `.venv/bin/python`.

## Known bugs (fix before trusting any rate)
- ✅ `cross_sections.py` — missing `import warnings` — FIXED 2026-06-22.
- 🐛 `rates.py:43` — `np.abs(oscillate(...))**2` keeps all 4 flavor amplitudes;
  must select τ component (index 2). (`simulate.py` already does this correctly;
  `rates.py` itself still unfixed but is not on the critical path.)

## Missing physics (the core deliverable)
The flux/osc/σ scaffolding exists, but the **measurement itself is unmodeled**.

### A. τ → μ decay kinematics module ✅ (`tau_decay.py`, alouette/TAUOLA)
1. ✅ τ → μ ν_μ ν_τ branching ratio (~17.4%) — realized by keeping muon-PID decays.
2. ✅ Boosted leptonic spectrum → daughter-μ energy and angle, via **alouette**
   (full kinematics: massive μ, radiative corrections, τ-spin↔μ correlation).
3. ✅ τ polarization = fixed helicity (−1 τ⁻ / +1 τ⁺, V−A). ⬜ energy-dependent CC
   polarization (Hagiwara et al.) is a later refinement.
4. ✅ τ production-point distribution along the beam (uniform, range-bounded) →
   off-axis displacement (`simulate.py`). v1 still assumes τ collinear (E_τ=E_ν).

### B. Detector geometry for off-axis muons ✅ v1 (`simulate.py`)
- ✅ (position r, angle θ) at the detector plane, range-bounded acceptance.
- ⬜ Detector model (size, threshold, resolution); 2D binning for the fit.

### C. Cross section integration ✅
- ✅ ν_τ CC σ now folded into `simulate.py` production weight.

### D. Background model ✅ v1 (`simulate_numu_cc_muons`, `compare.py`)
- ✅ Direct ν_μ CC muon (r, θ) distribution from DIS production; numeric core-vs-halo
  separability (S/B ratio, KS, Asimov q0). Both channels use full DIS production.
- ✅ Driven by **GENIE d²σ/dx dy tables** (`resources/diffxsec`, isoscalar CC).
- ⬜ Absolute normalization (n_target, area, exposure) for an absolute significance.
- ⬜ ν_e→ν_τ (δ_CP) and ν_e→ν_μ cross-flavor channels (tables available, deferred).
- ⬜ QE/RES (sub-DIS, low weight) — tables are DIS-only.

## Proposed build order (to discuss)
1. Fix the two bugs (quick, unblocks trustworthy rates).
2. Wire σ into `rates.py` → first real ν_τ CC event-rate estimate vs δ_CP.
3. Build the τ→μ decay kinematics module (A) — standalone + tested against known
   Michel-spectrum limits.
4. Add off-axis geometry / angular binning (B).
5. Add the ν_μ CC background model (D) and do the signal-vs-background study.

> Order is a proposal, not a decision — see [04-decisions.md](04-decisions.md) once we
> agree. Open feasibility questions live in [05-open-questions.md](05-open-questions.md).
