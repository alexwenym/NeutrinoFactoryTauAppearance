# Decisions log (ADR-lite)

Record physics/design decisions and their rationale as they're made. Each entry:
context → decision → rationale → consequences. Newest first.

## Template
```
## YYYY-MM-DD — <short title>
**Context:** what prompted the decision.
**Decision:** what we chose.
**Rationale:** why.
**Consequences:** what it implies / what it rules out.
```

## 2026-06-24 — DIS production from d²σ/dx dy splines; both channels upgraded
**Context:** the background core was unbuilt and the τ signal was collinear (E_τ=E_ν, no
production angle). The make-or-break observable is the muon angle, whose background width
comes entirely from DIS p_T — so it must be modeled from real DIS kinematics.
**Decision:** new `dis_kinematics.py` samples (x, y) from **double-differential
d²σ/dx dy** (user-supplied splines; analytic stand-in until delivered) and computes the
outgoing lepton energy + angle by exact massive-lepton kinematics; σ(E) by integrating
the same differential. **Both** the ν_μ core and the τ signal use full DIS production
(τ: E_τ=(1−y)E_ν + production angle composed with the decay kick). Normalization is
**consistent-relative** (core:halo ratio fixed by ∫flux·P·σ·BR·acc; absolute K deferred).
**Rationale:** the angular distribution is the deliverable, so model it from DIS, not a
parametrized kick; symmetric DIS treatment makes signal vs background apples-to-apples;
relative normalization answers the shape-separability question without a detector
benchmark. Per request, **no plots** in this cut — numeric metrics first.
**Consequences:** QE/RES (sub-DIS, low weight) unmodeled ⇒ conservative. ν_e
cross-flavor channels deferred. Absolute significance awaits the normalization layer.
**Update 2026-06-24:** the d²σ/dx dy source is now **GENIE tables** (G18_02a, isoscalar
CC, `resources/diffxsec` → `software/diffxsec/tables`), loaded by `TableDISModel` (the
analytic stand-in is now only a fallback). Real data corrected the τ production angle
upward (~93 mrad), giving S/B≈1.5e-3, KS(θ)≈0.08 — separation is modest, so it leans on
the 2D (r,θ)+energy shape. σ_τ now from the GENIE integral (cross-checked ~5% vs the
digitized total at 25–50 GeV).

## 2026-06-24 — alouette (TAUOLA) for τ→μ decay, pooled, fixed-helicity
**Context:** v0 `tau_decay.py` used the analytic massless unpolarized Michel spectrum.
The τ-spin↔μ correlation (sensitive to V−A and ultimately the production polarization)
is exactly the kind of off-axis kinematics we want modeled faithfully.
**Decision:** decay τ's with **alouette** (TAUOLA wrapper). Generate a **rest-frame
muon pool once** and sample+boost per event (not one alouette call per event). Default
τ polarization = **fixed helicity** −1 (τ⁻) / +1 (τ⁺); `polarization=0` recovers the
unpolarized check.
**Rationale:** TAUOLA gives the full decay (massive μ, radiative, spin correlation)
and the leptonic BR for free (keep muon-PID products). alouette is single-threaded
(~5e4 decays/s), so a pool keeps the per-event path vectorized. Fixed helicity is the
correct leading V−A behavior; energy-dependent CC polarization is a later refinement.
**Consequences:** adds an alouette dependency (Linux x86_64 wheel; no Apple-Silicon —
hence done on the cluster). Unpolarized alouette reproduces Michel v0; polarization on
shifts the observable forward — must carry the polarization assumption into sensitivity
claims. Validated against Michel limits (see [03-progress.md](03-progress.md)).

## 2026-06-18 — Lightweight rewrite instead of nuSQuIDS
**Context:** legacy `../tau_appearance/` study used nuSQuIDS + GENIE on the cluster.
**Decision:** this repo reimplements the chain analytically with a custom C
oscillator (`oschelper`), no heavy dependencies.
**Rationale:** portability, transparency, and tight control over the off-axis flux
and τ-decay modeling that nuSQuIDS doesn't provide out of the box.
**Consequences:** must validate our oscillation/flux against the legacy nuSQuIDS
results; we own the cross-section and decay modeling ourselves.

<!-- Next decision to record: agreed build order (see roadmap proposal). -->
