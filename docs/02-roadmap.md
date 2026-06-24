# Roadmap — gaps & planned work

Legend: ✅ done · 🚧 in progress · ⬜ not started · 🐛 bug

## Environment
- Project `.venv` (py3.9) with numpy + matplotlib; `oschelper` rebuilt into it
  (`pip install --no-build-isolation ./oschelper`). Run code via `.venv/bin/python`.
  The committed `oschelper/build/*.so` was a dead py3.12 build.

## Known bugs (fix before trusting any rate)
- ✅ `cross_sections.py` — missing `import warnings` — FIXED 2026-06-22.
- 🐛 `rates.py:43` — `np.abs(oscillate(...))**2` keeps all 4 flavor amplitudes;
  must select τ component (index 2). (`simulate.py` already does this correctly;
  `rates.py` itself still unfixed but is not on the critical path.)

## Missing physics (the core deliverable)
The flux/osc/σ scaffolding exists, but the **measurement itself is unmodeled**.

### A. τ → μ decay kinematics module ✅ (`tau_decay.py`)
1. ✅ τ → μ ν_μ ν_τ branching ratio (~17.4%).
2. ✅ Boosted leptonic (Michel) decay spectrum → daughter-μ energy and angle.
3. ✅ τ production-point distribution along the beam (uniform, range-bounded) →
   off-axis displacement (`simulate.py`). v1: unpolarized, massless-μ, τ collinear.

### B. Detector geometry for off-axis muons ✅ v1 (`simulate.py`)
- ✅ (position r, angle θ) at the detector plane, range-bounded acceptance.
- ⬜ Detector model (size, threshold, resolution); 2D binning for the fit.

### C. Cross section integration ✅
- ✅ ν_τ CC σ now folded into `simulate.py` production weight.

### D. Background model ⬜
- ⬜ Direct ν_μ CC muon angular distribution, to test whether off-axis μ angular
  resolution can separate signal from background (the make-or-break question).

## Proposed build order (to discuss)
1. Fix the two bugs (quick, unblocks trustworthy rates).
2. Wire σ into `rates.py` → first real ν_τ CC event-rate estimate vs δ_CP.
3. Build the τ→μ decay kinematics module (A) — standalone + tested against known
   Michel-spectrum limits.
4. Add off-axis geometry / angular binning (B).
5. Add the ν_μ CC background model (D) and do the signal-vs-background study.

> Order is a proposal, not a decision — see [04-decisions.md](04-decisions.md) once we
> agree. Open feasibility questions live in [05-open-questions.md](05-open-questions.md).
