# Code map — current state

Snapshot of what each file does **today** (2026-06-18). Update when code changes.

## Pipeline overview
```
flux.py  ──►  osc.py / oschelper  ──►  cross_sections.py  ──►  rates.py
(beam)        (P νe→ντ, matter)        (ντ CC σ)               (event rate vs δ_CP)
                                                               └─ τ→μ decay: NOT YET
```

## `flux.py` — muon-decay neutrino flux ✅
Analytic lab-frame flux from muon decay (Geer-style), **fully angle-resolved**.
- `costhcm_v(Emuon, costhlab)` — lab→CM angle transform.
- `numu_flux`, `nue_flux` — fluence with boost factor `M = 1/(G(1−β cosθ))` and
  polarization terms. `nue_flux` flips muon polarization for μ⁻.
- Three beam normalizations from arXiv:2407.12450:
  - `*_accelerator` — straight section pointing at target (rate = ring rotation rate).
  - `*_dump` — beam dump (rate = rep rate).
  - `*_baseline` — fixed 1e15 muon decays (order-of-magnitude reference).
- Constants: `MU_MASS`, `MU_TAU`, `MU_PER_BUNCH`, `REPRATE`, `CIRCUMFERENCE`, straight-
  section lengths.
- **Reusable as-is** for off-axis work — the angular machinery is exactly what's needed.

## `osc.py` — oscillation probabilities 🚧
4-flavor oscillation (3 active + 1 sterile) with matter effects.
- Builds PMNS matrix `U` (PDG 2020 convention) + NSI-style matter potential `MP`.
- Walks the **AK135 Earth density profile** along the chord (`DensityProfile*AK135.txt`).
- Delegates per-step Schrödinger evolution to compiled C ext `oschelper.do_osc`.
- Knobs: `deltaCP`, sterile (`th14/24/34`, `dm41`), anomalous matter (`eps*`).
- 🐛 Returns amplitudes for **all 4 flavors**, shape `(4, N_E)`. Callers must select
  the τ component (index 2) for P(ν_e→ν_τ). See `rates.py` bug below.

## `oschelper/` — C/NumPy extension ✅ (built)
Gray Putnam's extension integrating the Schrödinger equation step-by-step along the
baseline. `do_osc(E, steps, density, U, MP, dm21, dm31, dm41, nuind, dL)`.
- Initial flavor set by `nuind` (0=e, 1=μ, 2=τ). Output `(4, N_E)` complex amplitudes.
- Built `.so` present under `oschelper/build/`. Rebuild: `cd oschelper && pip install -e .`

## `tau_decay.py` — τ→μ leptonic decay kinematics ✅ (NEW)
- `sample_michel_x` — Michel spectrum 2x²(3−2x); `decay_to_muon(E_tau, rng)` →
  lab-frame muon (E_μ, θ, φ) via isotropic rest-frame emission + boost.
- Constants `M_TAU`, `M_MU`, `BR_TAU_TO_MU`. v1: unpolarized, massless-μ.

## `simulate.py` — MC driver for the (position, angle) distribution ✅ (NEW)
- `muon_range(E_mu)` — CSDA range in standard rock (a+bE model).
- `simulate_tau_muons(...)` — wires `flux.numu_flux` × P(ν_μ→ν_τ) (`osc`, nuind=1,
  index 2) × `totalXS.sigma` → `tau_decay` → range-bounded geometry; returns
  per-muon (r, θ, energies) at the detector plane. `__main__` saves
  `tau_muons_detector_plane.pdf`.
- ⬜ Direct-ν_μ CC core (for the with/without-osc comparison) not yet added.
- ⬜ Absolute normalization (n_target, area, time) — weights are relative for now.

## `cross_sections.py` — ν_τ CC cross sections ✅
- `totalXS` class: interpolates digitized ν_τ / ν̄_τ CC σ (arXiv:2409.01258 Fig. 4).
- `sigma_nu`, `sigma_nubar`, `sigma(x, is_nubar)`; units E [GeV], σ [1e-38 cm²].
- Data in `resources/nutau_xs.txt`, `resources/nutaubar_xs.txt`.
- ✅ `import warnings` bug fixed (2026-06-22).
- ⬜ Not yet wired into `rates.py` (but now used by `simulate.py`).

## `rates.py` — event rates 🚧
`TauAppearanceRates`: combines flux × solid angle × P over a set of δ_CP values.
- Inputs: baseline, off-axis angle, detector radius, E_μ, P_μ, μ charge, N_μ, δ_CP grid.
- `detector_solid_angle = π r² / L²` — single on-axis acceptance only.
- 🐛 `rates.py:43` keeps all 4 amplitudes (`osc_prob` is `(4,N)`), not P(ν_e→ν_τ);
  must index the τ component (2).
- ⬜ No cross-section multiplication → not yet an actual event rate.
- ⬜ No τ→μ decay step → does not yet model the real observable.

## Supporting files
- `Sandbox.ipynb` — scratch/exploration notebook (large).
- `figures.mplstyle` — plot styling.
- `DensityProfile{Depth,Density}AK135.txt` — Earth density profile for matter effects.
- `resources/` — digitized cross-section data + readme.

## Legacy (`../tau_appearance/`)
FASRC cluster environment: full nuSQuIDS, SQuIDS, GENIE (Generator), WCG cross
sections, plus `analysis/` prototype scripts (nuSQuIDS-based P(ν_α→ν_τ) heatmaps,
simpler on-axis flux). This repo is the lightweight, dependency-free rewrite. Synced
via rsync (`../tau_appearance/sync_commands.txt`).
