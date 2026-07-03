# NeutrinoFactoryTauAppearance

Phenomenology study: can a **neutrino factory** probe ν_τ appearance (ultimately δ_CP via
ν_e → ν_τ) by detecting the **off-axis muons from τ → μ ν ν** — the "halo" — above the
direct-ν_μ CC muon "core", in an idealized far detector?

## The measurement
A muon storage ring produces an intense, well-predicted ν_μ (+ ν̄_e) beam. Fired through
the Earth, ν_μ → ν_τ appearance produces τ's in the rock near a far detector; ~17% decay
to a muon that inherits a transverse kick (τ production angle + decay). So the **signal is
an off-axis "halo"** of muons, versus the on-axis "core" from direct ν_μ CC. We ask how
well an idealized LArTPC — recording each muon's (position, angle, energy) at a detection
plane — can statistically isolate the halo.

## Pipeline
flux (muon-decay beam) → oscillation (matter effects, AK135 Earth) → **GENIE d²σ/dxdy DIS**
(cross section + kinematics) → **alouette/TAUOLA** τ decay → rock propagation (muon range,
**Highland MCS**, energy loss) → idealized detection plane → **absolute normalization** →
halo-detection Asimov significance over **(r, θ, E)**. All numpy/matplotlib + a small C
oscillator; no heavy detector Monte Carlo.

## Key result
At a nominal DUNE-like benchmark (25 GeV, 1300 km) the halo is buried under a ~1000×
brighter core: **Z ≈ 0.9σ** for 100 m² × 10 yr. But **the significance grows ~linearly
with baseline**: the halo rate stays ~constant (oscillation growth ∝ L² cancels the 1/L²
flux loss — at these energies we sit on the oscillation *rising slope*, the first maximum
being beyond Earth's diameter), while the core background falls as 1/L². So **S/B ∝ L²,
Z ∝ L**.

An (L, E) scan (accelerator beam, 100 m² × 10 yr, MCS + LArTPC):
- **Best in-grid: E ≈ 57 GeV, L = 5000 km → Z ≈ 4.6σ** — near discovery for a modest
  detector; Z ∝ √(area·years) scales from there.
- Energy has a broad optimum ~50–90 GeV; the **muon-energy dimension ~doubles** the
  significance (the τ→μ spectrum is softer, from energy lost to two unseen neutrinos).
- A realistic CONUS baseline (Fermilab→Seattle, ~2723 km chord) reaches ~2.5σ.
- ν_μ **disappearance** (~5% core dimming) is a separate, far stronger signature — but not
  the τ-appearance measurement this study targets.

**Bottom line:** τ appearance via the μ halo is not hopeless — it is a **long-baseline**
measurement, favored by the longest Earth chord and mid-range energy.

## Layout
- `flux.py`, `osc.py` + `oschelper/`, `dis_kinematics.py`, `tau_decay.py` — physics inputs.
- `simulate.py` — the MC (τ→μ halo + ν_μ core), rock propagation + MCS.
- `detector.py` — idealized detection-plane response (footprint / smearing / threshold).
- `normalization.py` — absolute through-going-muon rate (target = upstream rock).
- `compare.py` — significance + CLI · `scan.py` — (L, E) scan · `plots.py` — figures.
- `docs/` — physics, design, decisions, open questions, and the progress log.

## Running
```bash
# one working point
python compare.py --emuon 57 --baseline 5000 --mcs --sigma-theta-mrad 3 --e-threshold 1 --e-res 0.15
# the (L, E) scan + heatmap
python scan.py
# detector-plane distributions
python plots.py
```
Cross sections come from GENIE d²σ/dxdy tables symlinked at `resources/diffxsec`
(→ `software/diffxsec/tables`, ~0.9 GB, gitignored; regenerated separately).
