# NeutrinoFactoryTauAppearance

A phenomenology study asking: can a **neutrino factory** (muon storage ring) probe
**ν_τ appearance** — and ultimately **δ_CP** via ν_e → ν_τ — by detecting the **off-axis
muons from τ → μ ν ν** (the "halo") above the direct-ν_μ CC muon "core", in an idealized
far detector?

The whole chain is deliberately lightweight: `numpy` + `matplotlib`, a small compiled C
oscillator, `alouette` (TAUOLA) for τ decay, and precomputed GENIE cross-section tables.
No heavy detector Monte Carlo.

---

## 1. The physics idea

A muon storage ring produces an intense, well-predicted neutrino beam from muon decay
(μ⁻ → e⁻ ν̄_e ν_μ → a ν_μ + ν̄_e beam). Fired through the Earth to a far detector,
**ν_μ → ν_τ** oscillation makes τ's in the rock just upstream of the detector. About
**17%** of τ's decay to a muon, and that muon inherits a **transverse kick** from the τ
production angle plus the decay. So the appearance signal is an **off-axis "halo"** of
muons, sitting on top of the on-axis **"core"** from direct ν_μ CC interactions.

- **Signal (halo):** ν_μ → ν_τ → τ → μ. Off-axis, softer (energy lost to two unseen ν's).
- **Background (core):** ν_μ → ν_μ → μ (direct CC DIS). On-axis, harder, ~1000× more
  numerous.

The question is whether an idealized LArTPC — recording each through-going muon's
**(transverse position r, arrival angle θ, energy E)** at a detection plane — can
statistically isolate the halo. The dominant channel ν_μ → ν_τ is **Milestone 1**
(the technique); ν_e → ν_τ + δ_CP is the harder **Milestone 2** (not yet built).

---

## 2. The pipeline, stage by stage

Data flows: **beam flux → oscillation → cross section + DIS kinematics → τ decay →
rock propagation → detection plane → absolute normalization → significance**, scanned
over (L, E).

### (a) Beam flux — `flux.py`
Analytic **Geer-style** muon-decay neutrino flux, fully angle-resolved, with the forward
relativistic boost (γ² on-axis enhancement — the reason a neutrino factory beam is
intense). Parametrized by the stored-muon energy `Emuon`. Three absolute-intensity
variants, with constants from a **muon collider design** (arXiv:2407.12450, IMCC):
`accelerator` (a ring straight section points at the detector), `dump`, and `baseline`
(a fixed 10¹⁵ decays reference). The scan uses `accelerator` ≈ **3.6×10¹⁹ useful
decays/yr at 25 GeV**, scaling as 1/E_muon (higher-energy muons live longer → fewer
decays per straight-section pass). See the beam caveats in §4.

### (b) Oscillation — `osc.py` + `oschelper/`
4-flavor (3 active + 1 optional sterile) PMNS oscillation **with matter effects**. Walks
the actual chord through the Earth using the **AK135 density profile**, and integrates the
Schrödinger equation step-by-step via a compiled C extension (`oschelper.do_osc`, Gray
Putnam). Carries δ_CP and the full matter potential. The halo weight uses P(ν_μ→ν_τ)
(component index 2); the core uses P(ν_μ→ν_μ) (index 1).

### (c) Cross section + DIS kinematics — `dis_kinematics.py`
Both the cross section **and** the outgoing-lepton kinematics come from **GENIE
double-differential DIS tables** d²σ/dx dy (tune G18_02a, E = 5–200 GeV), symlinked at
`resources/diffxsec`. For each neutrino, `TableDISModel` samples (Bjorken-x, inelasticity
y) from the **isoscalar (p+n)/2** CC table on its native (E, x, y) grid, then exact
**massive-lepton kinematics** give the outgoing lepton energy and scattering angle:
E_lep = (1−y)E, Q² = 2 M_N x y E, cosθ = (2EE′−Q²−m²)/(2E p′). Integrating the same table
gives σ(E), so kinematics and rate are self-consistent. Tables are DIS-only (no QE/RES;
low weight at the ~12–24 GeV bulk).

### (d) τ decay — `tau_decay.py`
τ → μ via **alouette (a TAUOLA wrapper)**: τ's are decayed at rest into a cached rest-frame
muon pool, then boosted along the τ flight direction. Keeping only muon-PID products
realizes the ~17.4% leptonic branching ratio; the full kinematics (massive μ, radiative
corrections, τ-spin↔μ correlation) are included. τ polarization is **fixed helicity**
(−1 for τ⁻, +1 for τ⁺; V−A).

### (e) Signal + background MC — `simulate.py`
Two parallel Monte Carlos, both importance-sampling the neutrino energy ∝ flux × P × σ:
- `simulate_tau_muons` (**halo**): DIS τ production (E_τ, production angle) → alouette
  τ→μ decay → the muon direction composed from the τ production angle and the decay kick.
- `simulate_numu_cc_muons` (**core**): a prompt DIS muon; ν_μ survival weight (an `osc_on`
  toggle gives the no-appearance hypothesis).

### (f) Rock propagation + geometry — inside `simulate._project_to_detector`
τ's/μ's are produced **all along the beam in the rock** (uniform production point d,
bounded by the muon range ≈ **1.9 m/GeV** in standard rock). A muon is accepted if it is
forward and has the range to reach the plane; it lands at transverse offset r = d·tanθ.
Optional (`mcs=True`) **Highland multiple Coulomb scattering** in the rock over the path
d/cosθ: a per-event, energy-loss-integrated σ_θ (~15–35 mrad) that smears **both** the
arrival angle **and**, correlated, the landing position — so muons produced further
upstream scatter more. The output muon energy is the **arrival energy** at the plane
(production minus ionization loss, ~0.53 GeV/m), i.e. what a tracker measures.

### (g) Detection plane — `detector.py`
The MC output **is** the truth-level detection-plane readout (per-muon r, θ, φ, E).
`apply_detector` adds an idealized-LArTPC response — a finite **footprint**, Gaussian
**angular resolution** (TPC-intrinsic only; rock MCS is upstream) and **energy resolution**,
and an **energy threshold** — all off by default ("record all, cut later"). Pheno-level:
simple Gaussian smearing, no detector MC.

### (h) Absolute normalization — `normalization.py`
Neutrino-induced **through-going muons**: the target is the **upstream rock** (the detector
is the tracker, not the target). By transverse translational invariance the muon flux
crossing the plane is a local per-area quantity,
Φ_μ = ∫ Φ_ν P σ · n_rock · R_eff · BR, with R_eff = dmax·accept_frac from the MC; total
counts = Φ_μ × area × exposure. Beam intensity from `flux.n_muon_decays`.

### (i) Significance — `compare.py`
The **halo-detection Asimov significance**: H1 = core + halo vs H0 = core (**same** core,
so it isolates τ *appearance*, not ν_μ disappearance), computed over the 3D **(r, θ, E)**
histogram (empty core bins floored by MC statistics). Reports the incremental
Z(θ) / Z(r,θ) / Z(r,θ,E) — the muon **energy dimension ~doubles** Z — plus absolute
counts and the (separate, far stronger) ν_μ disappearance. Z ∝ √(area·years).

### (j) Scan + plots — `scan.py`, `plots.py`
`scan.py` maps Z(r,θ,E) over an (L, E) grid → the `scan_LE` contour heatmap.
`plots.py` renders the detector-plane distributions (2D x–y hit maps, and the θ / r
shapes and absolute spectra) at any working point.

---

## 3. Key results

- **Nominal DUNE-like point (25 GeV, 1300 km, 100 m²×10 yr):** the halo is buried under a
  ~1000× brighter core (S/B ≈ 1.5×10⁻³) → **Z ≈ 0.9σ**. Marginal.
- **Significance grows ~linearly with baseline (Z ∝ L).** The halo rate is ~L-independent
  (P(ν_μ→ν_τ) ∝ L² on the oscillation **rising slope** cancels the 1/L² flux loss — the
  first oscillation maximum is at L/E ≈ 491 km/GeV, i.e. **beyond Earth's diameter** for
  E ≥ 25 GeV, so it is never reached/passed), while the core ∝ 1/L² drops. Net **S/B ∝ L²**.
- **(L, E) scan optimum: E ≈ 59 GeV, L = 5000 km → Z ≈ 4.66σ.** Broad energy optimum
  ~50–90 GeV; contours nearly vertical (L is the dominant lever).
- **Muon energy roughly doubles** the significance (softer τ→μ spectrum survives rock loss).
- **Fermilab→Seattle** (2723 km chord) ≈ **2σ** — a realistic CONUS baseline.
- Z ∝ √(area·years), so a larger detector or longer exposure scales all of the above up.

**Bottom line:** τ appearance via the μ halo is not hopeless — it is a **long-baseline**
measurement, favored by the longest Earth chord and mid-range energy.

---

## 4. Assumptions & caveats

- **Beam.** Neutrinos from muon decays in a storage-ring/collider straight section (Geer
  spectrum — realistic). The absolute intensity uses **muon-collider design numbers**
  (1.8×10¹² μ/bunch, 5 Hz, 10 km ring) applied across 25–200 GeV with only 1/E lifetime
  scaling — **not a self-consistent per-energy facility design**, and the straight-section
  / pointing factors are order-of-magnitude. Since **Z ∝ √(N_μ·area·years)**, the **map
  shape (Z ∝ L) is robust** but the **absolute σ is good to a factor of a few**.
- **Hypothesis test.** Z is the halo-*appearance* significance with the core held fixed
  (perfectly known); flux/normalization **systematics/nuisances are not yet included**.
- **Cross sections.** GENIE DIS only (no QE/RES); low weight at the relevant energies.
- **MCS.** Analytical Highland (Gaussian core); full Molière tails + range straggling
  would be a PROPOSAL-level upgrade.
- **Flux on-axis** (beam angular divergence, ~few mrad, neglected); τ polarization
  fixed-helicity (energy-dependent CC polarization is a refinement).
- **ν_e / δ_CP channel not yet modeled** — the README's ultimate goal (Milestone 2).

See `docs/05-open-questions.md` and `docs/04-decisions.md` for the full list.

---

## 5. Repository layout

| Module | Role |
|--------|------|
| `flux.py` | muon-decay beam flux + intensity variants |
| `osc.py`, `oschelper/` | 4-flavor oscillation with matter effects (AK135 chord) |
| `dis_kinematics.py` | GENIE d²σ/dx dy: cross section + DIS lepton kinematics |
| `tau_decay.py` | τ → μ via alouette/TAUOLA (fixed-helicity) |
| `simulate.py` | halo + core MC, rock propagation, Highland MCS, arrival energy |
| `detector.py` | idealized detection-plane response (footprint/smearing/threshold) |
| `normalization.py` | absolute through-going-muon rate (target = upstream rock) |
| `compare.py` | halo-detection significance over (r, θ, E) + CLI |
| `scan.py` | (L, E) significance map |
| `plots.py` | detector-plane figures |
| `cross_sections.py` | digitized ν_τ CC σ (now a cross-check of the GENIE integral) |
| `docs/` | physics, design, decisions, open questions, progress log |

Most modules self-test / print sanity numbers when run directly (`python <module>.py`).

---

## 6. Running

```bash
# environment (Linux x86_64): venv + oschelper rebuilt with a C++ linker
python3 -m venv .venv && . .venv/bin/activate
pip install numpy matplotlib alouette
LDSHARED="g++ -shared" CC=g++ pip install --no-build-isolation ./oschelper

# a single working point (the optimum), realistic detector + MCS
python compare.py --emuon 59 --baseline 5000 --mcs --sigma-theta-mrad 3 --e-threshold 1 --e-res 0.15

# the (L, E) significance scan + heatmap
python scan.py

# detector-plane distributions at a chosen point
python plots.py
```

**Cross-section data:** GENIE d²σ/dx dy tables are symlinked at `resources/diffxsec`
(→ `software/diffxsec/tables`, ~0.9 GB, gitignored; produced separately — see the memory
note on the GENIE setup). Run from the repo root (the oscillator reads the AK135 profiles
and the DIS tables relative to CWD).
