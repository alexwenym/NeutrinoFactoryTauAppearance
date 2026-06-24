# Physics goal & measurement strategy

## The question
A muon storage ring (neutrino factory) produces a clean, well-predicted neutrino
beam from muon decay:

- μ⁻ → e⁻ ν̄_e ν_μ  →  beam of **ν_μ** and **ν̄_e**
- μ⁺ → e⁺ ν_e ν̄_μ  →  beam of **ν̄_μ** and **ν_e**

We fire this beam through the Earth and look for **ν_τ appearance** in a far
detector. The **ν_e → ν_τ** channel carries the CP-violating phase **δ_CP**, so the
appearance rate vs δ_CP is the observable that probes CP violation.

## The measurement strategy (the heart of this project)
We do **not** try to reconstruct the τ directly. Instead we detect ν_τ appearance
through the **leptonic decay**

    τ → μ ν_μ ν_τ      (branching ratio ≈ 17.4%)

and measure the resulting **off-axis muons** — their angular displacement / "decay
angle" — precisely in a far detector.

**Why this is attractive:** a clean μ track is far easier to detect than a τ.
**The price:**
1. The ~17% leptonic branching ratio (rate hit).
2. Two unseen neutrinos in the decay smear the muon energy and angle, so the decay
   angle is only **statistically** reconstructable — not event-by-event.

## Signal vs background — the make-or-break question
The dominant background is **direct ν_μ CC interactions** producing muons. The whole
measurement hinges on:

> Can off-axis muon angular resolution statistically separate τ-daughter muons
> (which acquire transverse momentum from the τ decay) from direct-ν_μ CC muons
> (which point back along the beam)?

Nothing in the current code addresses this separation yet — see
[02-roadmap.md](02-roadmap.md).

## Physics ingredients needed end-to-end
1. **Beam flux** — angle-resolved ν_e flux from muon decay. ✅ (`flux.py`)
2. **Oscillation** — P(ν_e → ν_τ) through the Earth with matter effects, δ_CP
   dependence. 🚧 (`osc.py`; flavor-selection bug, see roadmap)
3. **ν_τ CC cross section** — production of τ in the detector. ✅ data, ⬜ not wired in
   (`cross_sections.py`)
4. **τ → μ decay kinematics** — BR + boosted leptonic (Michel) spectrum giving
   daughter-μ energy & angle. ⬜ **entirely unwritten — the core deliverable.**
5. **Geometry** — τ production-point distribution → off-axis μ displacement at the
   detector; angular binning of arriving muons. ⬜
6. **Background model** — direct ν_μ CC muon angular distribution for comparison. ⬜

## Reference papers
- Beam configs & muon numbers: arXiv:2407.12450 (Tables 1.1, 6.4)
- ν_τ / ν̄_τ cross sections: arXiv:2409.01258 (Fig. 4, digitized in `resources/`)
- Muon-decay flux: Geer, Phys. Rev. D 57, 6989 (1998)
- PMNS convention: PDG 2020 neutrino mixing review
