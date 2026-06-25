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
- Rebuild per-machine (the `.C` source is C++, so force g++ for the link step):
  `LDSHARED="g++ -shared" CC=g++ pip install --no-build-isolation ./oschelper`.
  A plain build links with `gcc -shared` and import fails with
  `undefined symbol: _ZTVSt9basic_ios...` (missing libstdc++).

## `tau_decay.py` — τ→μ leptonic decay kinematics ✅ (alouette)
- `decay_to_muon(E_tau, rng, nusign=1, polarization=None)` → lab-frame muon
  (E_μ, θ, φ). Uses **alouette** (TAUOLA): taus decayed at rest, muon-PID products
  kept (realizes the leptonic BR), boosted along the beam. Full muon kinematics
  (massive μ, radiative corrections, τ-spin↔μ-direction correlation).
- **Rest-frame pool** (`POOL_SIZE=100k`, cached in `_POOLS` by charge+helicity):
  built once via `_build_pool`/`_get_pool`, then sampled+boosted per event — alouette
  is single-threaded (~5e4 decays/s) so the pool keeps the hot path vectorized.
- τ polarization = fixed helicity: default −1 for τ⁻ / +1 for τ⁺ (V−A);
  `polarization=0.0` → unpolarized (cross-checks Michel ⟨x⟩=0.70).
- `sample_michel_x` kept as the analytic massless reference. Constants `M_TAU`,
  `M_MU`, `BR_TAU_TO_MU`. v1 still assumes τ collinear with ν (E_τ=E_ν), set in
  `simulate.py`.

## `dis_kinematics.py` — DIS production kinematics + σ ✅ (GENIE tables)
- `sample_dis_lepton(E_nu, rng, channel)` → outgoing CC lepton `(E_lep, θ, φ)` (θ from
  the beam) by sampling (x, y) from d²σ/dx dy and exact massive-lepton kinematics
  (`Q²=2 M_N x y E`). `sigma(E, channel)` integrates the same differential (1e-38 cm²).
- **`TableDISModel`** (default) loads GENIE CC tables from `resources/diffxsec`
  (→ `software/diffxsec/tables`), **isoscalar (p+n)/2**, native 35×50×45 (E,x,y) grid
  (E 5–200 GeV log, x 1e-3–0.95 log, y 0.01–0.99 lin), units cm². `DISSampler` builds
  cell-width-weighted per-E-node 2D CDFs and samples within cells.
- `AnalyticDISModel` (LO valence+sea) is a fallback; `set_model(...)` swaps it.
- Channels `numu, numubar, nutau, nutaubar, nue, nuebar`. Named `dis_kinematics` (not
  `dis`) to avoid shadowing stdlib `dis`. (Tables are DIS-only ≥5 GeV; no QE/RES.)

## `simulate.py` — MC drivers for the (position, angle) distribution ✅
- `muon_range(E_mu)` — CSDA range in standard rock (a+bE model).
- `_project_to_detector(E_mu, θ, dmax, n_mc, rng)` — shared production-point + forward/
  range acceptance + `r=d·tanθ`. `_compose_directions(...)` — composes a parent
  production angle with a decay angle (unit-tested).
- `simulate_tau_muons(...)` — `numu_flux` × P(ν_μ→ν_τ) (`osc` index 2) × `dis.sigma`
  (GENIE ν_τ DIS) → **DIS τ production** (`dis_kinematics`, E_τ=(1−y)E_ν + angle) →
  `tau_decay` → compose → geometry. Returns per-muon (r, θ, E) + `w_integral` + `BR`.
- `simulate_numu_cc_muons(..., osc_on=True)` — direct-ν_μ CC core: P(ν_μ→ν_μ) (index 1,
  or 1 if `osc_on=False`) × `dis.sigma(numu)` → prompt DIS muon → geometry.
- ⬜ Absolute normalization (n_target, area, time) — weights relative; ratio is correct.

## `cross_sections.py` — ν_τ CC cross sections ✅
- `totalXS` class: interpolates digitized ν_τ / ν̄_τ CC σ (arXiv:2409.01258 Fig. 4).
- `sigma_nu`, `sigma_nubar`, `sigma(x, is_nubar)`; units E [GeV], σ [1e-38 cm²].
- Data in `resources/nutau_xs.txt`, `resources/nutaubar_xs.txt`.
- ✅ `import warnings` bug fixed (2026-06-22).
- Now used only as a **cross-check** of the GENIE ν_τ DIS σ (see `dis_kinematics`
  self-test); `simulate.py` takes σ from the GENIE tables.

## `rates.py` — event rates 🚧
`TauAppearanceRates`: combines flux × solid angle × P over a set of δ_CP values.
- Inputs: baseline, off-axis angle, detector radius, E_μ, P_μ, μ charge, N_μ, δ_CP grid.
- `detector_solid_angle = π r² / L²` — single on-axis acceptance only.
- 🐛 `rates.py:43` keeps all 4 amplitudes (`osc_prob` is `(4,N)`), not P(ν_e→ν_τ);
  must index the τ component (2).
- ⬜ No cross-section multiplication → not yet an actual event rate.
- ⬜ No τ→μ decay step → does not yet model the real observable.

## `detector.py` — idealized detection-plane ("LArTPC") layer ✅ (off by default)
- The sim output dict **is** the detection-plane readout (per-muon truth θ, r, φ, E_μ for
  every muon reaching the plane). `apply_detector(out, ...)` is a downstream,
  truth-preserving response: optional **footprint** (circular `radius` or rectangular
  `half_width`/`half_height` on the crossing position), **Gaussian smearing**
  (`sigma_theta` per-projection — fold TPC + rock-MCS in quadrature; `sigma_E_frac`), and
  `E_threshold`. All off by default → no-op passthrough. Adds `det_frac` (survival).
- `compare.run_comparison(..., detector=dict(...))` applies it to every sample; the
  channel rate (`_channel_weight`) folds in `det_frac`. Pheno-level only (no detector MC).

## `compare.py` — numeric core-vs-halo separability ✅ (no plots)
- `run_comparison(...)` runs halo + core(osc) + core(no-osc); reports S/B rate ratio
  (`W_halo/W_core`), θ-shape metrics (KS, tail fractions), and a per-exposure Asimov
  `q0` (absolute Z needs the deferred normalization K). Plots intentionally deferred.
- CLI for L–E scans: `python compare.py --emuon 50 --baseline 2000 [--nusign ±1
  --deltacp --n-mc --seed]` (baseline in km). Both knobs thread through flux/osc/σ/kin.

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
