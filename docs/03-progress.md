# Progress log

Newest first. One entry per meaningful chunk of work.

## 2026-06-24 (cluster) — multiple Coulomb scattering (MCS) in rock
- Added analytical **Highland MCS** to the muon propagation (`simulate._project_to_detector`,
  toggled by `mcs=` on the sims / `run_comparison` / `--mcs`; default off to preserve the
  truth-level baseline). Per-event, path-dependent: σ_θ integrated over the rock path
  d/cos θ **including ionization energy loss** (E_arr = E0 − dE/dx·L), so muons produced
  further upstream — and softer muons — scatter more. Smears **both** the arrival angle
  (RMS σ_θ per projection plane) **and**, correlated (ρ=√3/2), the landing position
  (RMS L·σ_θ/√3). Constants: X0=10.02 cm, dE/dx=5.3 MeV/cm (standard rock).
- Validated: σ_θ ≈ 12–35 mrad/plane for muons comfortably crossing the rock (matches the
  `docs/06` ~20 mrad estimate), diverging for near-stopping muons (floored). Effect on the
  result: halo/core median θ each broaden ~7–9 mrad (MCS adds in quadrature to the ~100
  mrad DIS angle); **separability preserved** (KS(θ) 0.077→0.084) — MCS does not wash out
  the halo. Full Moliere tails + range straggling remain a PROPOSAL-level upgrade.

## 2026-06-24 (cluster) — beam/baseline CLI + idealized detector-plane layer
- Added a CLI to `compare.py` (`--emuon`, `--baseline` [km], `--nusign`, `--deltacp`,
  `--n-mc`, `--seed`) for L–E scans; `Emuon`/`baseline` were already parameters threading
  through flux/osc/σ/kinematics (verified: 50 GeV/2000 km narrows angles ~1/E, shifts the
  osc point). No physics change.
- New `detector.py`: the simulation output already **is** a truth-level detection plane
  (per-muon θ, r, φ, E_μ). `apply_detector(out, ...)` is an off-by-default downstream
  "LArTPC" response — footprint (record-all / cut-later), Gaussian θ & E smearing
  (sigma_theta folds TPC + rock-MCS), energy threshold; adds `det_frac`. Wired as an
  optional `detector=` kwarg in `run_comparison` (`_channel_weight` folds in `det_frac`).
  Per the user: pheno only, plane formalized first, smearing/footprint left off for now.

## 2026-06-24 (cluster) — wired in GENIE differential cross sections
- **Replaced the analytic DIS stand-in with GENIE double-differential d²σ/dx dy tables**
  (tune G18_02a, E=5–200 GeV; produced by a separate agent via a custom `gDISDiffXSec`
  app). Symlinked `resources/diffxsec -> software/diffxsec/tables`. `dis_kinematics.py`
  gained `TableDISModel` (loads CC **isoscalar (p+n)/2** on the native 35×50×45 (E,x,y)
  grid; default model) and a **grid-based `DISSampler`** (cell-width-weighted CDFs);
  `AnalyticDISModel` kept as a fallback (`set_model`). σ(E) now integrates the tables.
- `simulate_tau_muons` now takes σ_τ from the GENIE tables (`dis.sigma`), dropping the
  digitized `totalXS` dependency (kept only as a self-test cross-check: GENIE DIS vs
  digitized total ν_τ CC agree ~5% at 25–50 GeV, ~30% near 10 GeV threshold).
- **Key correction from real data:** the **τ production angle is ~93 mrad** (the τ mass
  forces Q²≳m_τ²), not the stand-in's ~12 mrad — the τ is not ultra-forward. Updated
  numbers (Emuon=25, L=1300): S/B≈**1.5e-3**, halo median θ 116, core 102 mrad,
  **KS(θ)=0.076** (cores are broad from DIS p_T, so separation is modest — leans on the
  full 2D shape). ν_τ ⟨y⟩=0.37 (τ-mass-suppressed vs numu 0.51).

## 2026-06-24 (cluster) — DIS background core + DIS production upgrade
- **Built the direct-ν_μ CC "core"** background and **upgraded both channels to full DIS
  production** (signal τ was previously collinear). New `dis_kinematics.py`:
  samples (x, y) from a double-differential d²σ/dx dy model → outgoing lepton
  `(E_lep, θ, φ)` via exact massive-lepton kinematics (`Q²=2 M_N x y E`,
  `cosθ=(2EE'−Q²−m²)/(2E p')`); σ(E) by integrating the same differential. **Renamed
  from `dis.py`** (shadowed stdlib `dis`, broke `inspect`/numpy import).
- ⚠️ **The d²σ/dx dy splines (νμ,ν̄μ,ντ,ν̄τ,νe,ν̄e CC) are not yet delivered.** An
  `AnalyticDISModel` (LO valence+sea) stand-in drives the pipeline meanwhile; swap via
  `dis_kinematics.set_model(SplineDISModel(...))` (`load_xsec_splines` is a stub pending
  the file format).
- `simulate.py`: factored shared geometry into `_project_to_detector`; added
  `_compose_directions` (composes the τ production angle with the τ→μ decay kick, unit-
  tested); upgraded `simulate_tau_muons` to DIS τ production; added
  `simulate_numu_cc_muons` (prompt DIS muon, ν_μ survival `|amp[1]|²`, `osc_on` toggle
  for the brighter no-appearance core). Both now return `w_integral`.
- `compare.py` (new, **no plots** per request): numeric core-vs-halo separability —
  S/B rate ratio, θ-shape metrics (KS, tail fractions), and a per-exposure Asimov q0.
- **First numbers** (stand-in, Emuon=25, L=1300 km): S/B = W_halo/W_core ≈ **9.6e-4**;
  halo median θ ≈ 97 mrad (90% 261), core ≈ 79 mrad (90% 220); KS(θ)=0.10. DIS
  inelasticity lowers E_τ (~0.55 E_ν) so the decay cone widens (53→97 mrad). Stand-in
  ⟨x⟩≈0.27 is hard, so the core width is an upper estimate — real splines will sharpen it.
  Validation: ⟨y⟩ 0.48/0.33 (ν/ν̄); σ_μ/E=0.677 (PDG); σ_τ rises with E (threshold).

## 2026-06-24 (cluster)
- **Swapped the analytic Michel τ-decay in `tau_decay.py` for alouette** (TAUOLA
  wrapper, v1.0.1). Run taus at rest, keep muon-PID products (realizes the leptonic
  BR ≈ 0.177), boost a precomputed **rest-frame pool** (100k) along the beam — keeps
  the per-event path vectorized despite alouette being single-threaded (~5e4 decays/s).
  τ polarization = **fixed helicity** (−1 for τ⁻, +1 for τ⁺, V−A); `polarization=0`
  recovers the unpolarized spectrum. `decay_to_muon` gained optional `nusign`/
  `polarization` kwargs (interface preserved); `simulate.py` now passes `nusign` so
  antineutrino running gets τ⁺ polarization.
- **Verified:** unpolarized alouette ⟨x⟩=0.707 ≈ Michel 0.70; ⟨cosθ*⟩ = 0 / ∓0.11 for
  helicity 0 / ∓1; lab ⟨E_μ⟩/E_τ = 0.354 unpol → 0.405 with V−A. End-to-end `simulate.py`
  unpolarized reproduces Michel v0 (acceptance 0.264 vs 0.261, median θ 59 vs 60 mrad,
  median r 0.43 m); V−A shifts muons forward (acceptance 0.302, median θ 53 mrad).
- **Fixed the `oschelper` build on the cluster:** the installed `.so` failed import
  (`undefined symbol: _ZTVSt9basic_ios...`, a libstdc++ vtable) — the `.C` source is
  C++ but setuptools linked with `gcc -shared`. Rebuilt with `LDSHARED="g++ -shared"
  CC=g++ pip install --no-build-isolation ./oschelper`. `do_osc` now imports.
- Env note: cluster venv is **py3.8** (docs elsewhere said 3.9); alouette 1.0.1 was
  already installed.

## 2026-06-24 (laptop)
- Forked nickkamp1's repo to **alexwenym/NeutrinoFactoryTauAppearance**; rewired
  remotes (`origin`=fork, `upstream`=nickkamp1). Committed the v1 sim + docs (`b5008d9`)
  and pushed to the fork. Build artifacts gitignored.
- Decided to do the **alouette** (TAUOLA) τ-decay swap **on the Linux cluster**
  (alouette has no Apple-Silicon wheel; Linux x86_64 wheel exists).
- Wrote [07-cluster-handoff.md](07-cluster-handoff.md) so the cluster agent can pick up:
  env setup, the alouette API, the two open decisions (τ polarization; direct vs pool),
  and verification steps.

## 2026-06-22
- **First working build** of the ν_μ→ν_τ appearance muon distribution at the detector
  plane, built on the existing flux/osc/xs code. New files:
  - `tau_decay.py` — τ→μ leptonic decay kinematics (Michel sampling + boost).
    Self-tested: ⟨x⟩=0.700, ⟨E_μ⟩/E_τ=0.350, θ∝1/E_τ (all match analytics).
  - `simulate.py` — MC driver: flux × P(ν_μ→ν_τ) × σ → τ → μ → range-bounded
    geometry → (position r, angle θ) at the detector plane. Saves
    `tau_muons_detector_plane.pdf`.
  - Fixed the `import warnings` bug in `cross_sections.py` (blocked running).
  - Set up project `.venv` (numpy, matplotlib) and rebuilt `oschelper` for py3.9
    (the committed `.so` was a dead py3.12 build).
  - First result (E_μ=25 GeV, L=1300 km): production region 46 m, θ median ~60 mrad
    (90% < ~146 mrad), r median ~0.43 m (out to ~1.7 m). Off-axis angle ≫ MCS scale.

## 2026-06-18
- Refined sim design after user input: observable is **2D (position + angle)** at the
  detector plane; τ's produced **all along the beam**; comparison is with vs without
  ν_μ→ν_τ oscillation; **Milestone 1 = ν_μ→ν_τ** (technique), δ_CP/ν_e→ν_τ is
  Milestone 2. Muon range kept in v1 to bound the production region. Resolved the
  geometry open questions.
- Drafted [06-simulation-design.md](06-simulation-design.md): the MC chain, v1
  simplifications (pure kinematics, assume angular resolution), and the MCS/muon-range
  feasibility caveat. No code changed.
- Set up this `docs/` tree (physics, code map, roadmap, decisions, open questions).
  No code changed.
- Clarified the measurement strategy with the user: detect ν_τ appearance via the
  leptonic decay τ→μν_μν_τ by precisely measuring **off-axis muons** in a far
  detector (not direct τ reconstruction).
- Read the full codebase + legacy `../tau_appearance/`; identified two bugs and the
  major missing piece (τ→μ decay module). See [02-roadmap.md](02-roadmap.md).

## Pre-existing (from git history, before this session)
- `a722c0c` first commit
- `c37666f` cross section digitized
- `f48486d` cross section class
- `aa4182a` cross sections verification
- `1d0b3b4` merge PR #1
