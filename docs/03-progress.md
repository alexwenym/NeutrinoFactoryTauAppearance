# Progress log

Newest first. One entry per meaningful chunk of work.

## 2026-06-24
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
