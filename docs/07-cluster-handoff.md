# Cluster handoff — where the laptop agent left off

For the agent picking this up on the Linux cluster. Written 2026-06-24.

## ✅ DONE 2026-06-24 (on the cluster): the alouette swap
`tau_decay.py` now decays τ's with **alouette** (pooled, fixed-helicity V−A) and is
validated against the Michel limits — see [03-progress.md](03-progress.md) and
[04-decisions.md](04-decisions.md). The `oschelper` build was also fixed (libstdc++
link, `LDSHARED="g++ -shared"`). The rest of this doc is kept as the original handoff
record; **the live next steps are now the deferred physics in
[02-roadmap.md](02-roadmap.md)** (direct-ν_μ CC background core for the with/without-osc
comparison, absolute normalization, MCS realism).

## TL;DR (original)
A working v1 of the ν_μ→ν_τ appearance-muon simulation exists and runs (commit
`b5008d9` on `origin/main` = `alexwenym/NeutrinoFactoryTauAppearance`). The **next
task** is to replace the analytic Michel τ-decay in `tau_decay.py` with **alouette**
(a TAUOLA wrapper), which installs cleanly on Linux x86_64 but had no Apple-Silicon
wheel — which is why it was deferred to the cluster.

## What works right now
- `tau_decay.py` — analytic τ→μ leptonic (Michel) decay kinematics. Self-tested:
  ⟨x⟩=0.700, ⟨E_μ⟩/E_τ=0.350, opening angle ∝ 1/E_τ. **This is what alouette replaces.**
- `simulate.py` — MC driver: `flux.numu_flux` × P(ν_μ→ν_τ) (`osc.oscillate`, nuind=1,
  τ component index 2) × `totalXS.sigma` → τ → μ → range-bounded geometry →
  (transverse position r, arrival angle θ) at the detector plane. Saves
  `tau_muons_detector_plane.pdf`. Runs in <1 s for 1e6 events.
- `cross_sections.py` — `totalXS`; the `import warnings` bug is fixed.
- `osc.py` + `oschelper/` C extension — P(ν_α→ν_τ) with matter effects (AK135).

See [01-code-map.md](01-code-map.md) for the full module map and
[06-simulation-design.md](06-simulation-design.md) for the physics/MC design.

## Environment setup on the cluster (x86_64 Linux)
```bash
git clone https://github.com/alexwenym/NeutrinoFactoryTauAppearance.git
cd NeutrinoFactoryTauAppearance
python3 -m venv .venv && . .venv/bin/activate
pip install numpy matplotlib alouette          # alouette: Linux x86_64 wheel exists
LDSHARED="g++ -shared" CC=g++ pip install --no-build-isolation ./oschelper  # see note
```
**`oschelper` C++ link fix:** the `.C` source is C++; a plain build links with
`gcc -shared` and import fails with `undefined symbol: _ZTVSt9basic_ios...` (missing
libstdc++). Forcing `LDSHARED="g++ -shared"` fixes it. (Cluster venv is py3.8; alouette
1.0.1 was already present in the existing `.venv`.)
Notes / gotchas carried over from the laptop:
- `oschelper` MUST be rebuilt per-machine. The committed history had a dead build for
  another Python; the build dir is now gitignored. `--no-build-isolation` is required
  so `setup.py` sees numpy.
- `oschelper` falsely "imports" as a namespace package (the source dir) without
  `do_osc` if not actually built — verify with
  `python -c "import oschelper; print(hasattr(oschelper,'do_osc'))"` → must be True.
- Run via the venv's python. `osc.py` reads `DensityProfile*AK135.txt` from CWD, so run
  from the repo root.

## The next task: swap in alouette
Replace the analytic decay in `tau_decay.py` with alouette while keeping
`simulate.py`'s interface (`decay_to_muon(E_tau, rng) -> (E_mu, theta, phi)`) intact.

alouette API (v1.0.1):
```python
import alouette
p = alouette.decay(mode=0, pid=15, momentum=(px,py,pz), polarisation=(sx,sy,sz))
# p.pid : PDG ids of products, p.P : (N,4) four-momenta (px,py,pz,E) GeV, p.weight
```
- pid=15 is τ⁻ (use −15 for τ⁺ / antineutrino running). Select the μ (pid ±13) from
  the products. One τ per call; **single-threaded** (no multiprocessing inside).
- Decay in the τ **rest frame** (`momentum=(0,0,0)`) with the chosen polarisation, then
  boost by E_τ along the beam — mirrors the current code and decouples decay from
  kinematics. (Alternatively pass the lab momentum directly and let alouette boost.)

### Two open decisions (pick with the user, captured in [05-open-questions.md](05-open-questions.md))
1. **τ polarization** — the main reason to use alouette over Michel.
   - Recommended v1: **fixed helicity −1** for τ⁻ (V−A: ν gives left-handed τ⁻),
     +1 for τ⁺. Energy-dependent CC polarization (Hagiwara et al.) is a later refinement.
   - Sanity step: run alouette **unpolarised** first and confirm it reproduces the
     analytic Michel ⟨x⟩=0.7 before turning polarization on.
2. **Integration style** — direct `alouette.decay()` per event, vs precompute a
   rest-frame decay **pool** once and sample+boost it. Pool is faster (single-threaded
   calls are the bottleneck) and keeps runtime alouette-free; direct is simpler. On a
   cluster either is fine; decide based on event counts.

### Verification after the swap
- Unpolarised alouette vs analytic Michel: ⟨x⟩→0.70, ⟨E_μ⟩/E_τ→0.35 (cross-check).
- Re-run `python simulate.py`; the (r, θ) plot should be close to the Michel version
  unpolarised, and shift in a definite direction once polarization is on.
- Then proceed to the deferred physics in [02-roadmap.md](02-roadmap.md): the direct-ν_μ
  CC core (the with/without-oscillation comparison), absolute normalization, MCS realism.

## Git layout
- `origin` = `alexwenym/NeutrinoFactoryTauAppearance` (your fork) — push here.
- `upstream` = `nickkamp1/NeutrinoFactoryTauAppearance` — `git fetch upstream` for Nick's
  updates.
- Build artifacts (`.venv/`, `oschelper/build/`, `*.egg-info/`, the output PDF) are
  gitignored.
