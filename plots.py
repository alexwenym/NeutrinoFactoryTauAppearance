import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from simulate import simulate_tau_muons, simulate_numu_cc_muons
from detector import apply_detector
import normalization

# Visualize the muon hits on the idealized detection plane: the direct-nu_mu CC "core"
# (background) vs the tau-decay "halo" (signal). Shows why the halo is hard to isolate --
# the two overlap heavily. Realistic: MCS in rock + nominal LArTPC response.

_DET = dict(sigma_theta=0.003, sigma_E_frac=0.15, E_threshold=1.0)  # nominal LArTPC


def make_plots(Emuon=25.0, baseline=1300e3, n_mc=1_500_000, seed=0,
               out="detector_plane_hits"):
    halo = apply_detector(simulate_tau_muons(Emuon=Emuon, baseline=baseline,
                                             n_mc=n_mc, seed=seed, mcs=True), seed=seed + 5, **_DET)
    core = apply_detector(simulate_numu_cc_muons(Emuon=Emuon, baseline=baseline,
                                                 n_mc=n_mc, seed=seed + 1, mcs=True), seed=seed + 6, **_DET)

    def xy(o):
        return o["r"] * np.cos(o["phi"]), o["r"] * np.sin(o["phi"])

    cx, cy = xy(core)
    hx, hy = xy(halo)
    R = 2.0  # m, plot half-range for the hit maps
    xyb = np.linspace(-R, R, 91)

    # absolute counts for the "what the detector sees" panel (100 m^2 x 10 yr)
    area_m2, years = 100.0, 10.0
    N_core = normalization.counts(core, Emuon, "accelerator", area_m2, years)
    N_halo = normalization.counts(halo, Emuon, "accelerator", area_m2, years)

    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.5))
    thb = np.linspace(0, 350, 70)   # arrival angle [mrad]
    rb = np.linspace(0, R, 70)      # transverse offset [m]

    # --- (x, y) detector-plane hit maps, log color to show the spread not just the core ---
    for a, X, Y, ttl, cm in ((ax[0, 0], cx, cy, "direct $\\nu_\\mu$ CC core (background)", "viridis"),
                             (ax[0, 1], hx, hy, r"$\tau\to\mu$ halo (signal)", "magma")):
        h = a.hist2d(X, Y, bins=xyb, cmap=cm, cmin=1, norm=LogNorm())[3]
        fig.colorbar(h, ax=a, label="muons / bin (log)")
        a.set_xlabel("x at detector plane [m]")
        a.set_ylabel("y at detector plane [m]")
        a.set_title(ttl)
        a.set_aspect("equal")

    # --- normalized shapes: angle (discriminant) and radial offset ---
    for a, key, scale, xlabel, ttl in (
            (ax[1, 0], "theta", 1e3, r"arrival angle $\theta$ [mrad]",
             "angle shapes: halo broader (the discriminant)"),
            (ax[1, 1], "r", 1.0, "transverse offset $r$ [m]",
             "radial shapes: core wider (higher-E, from further upstream)")):
        bins = thb if key == "theta" else rb
        a.hist(core[key] * scale, bins=bins, density=True, histtype="step",
               color="C0", lw=1.8, label="core")
        a.hist(halo[key] * scale, bins=bins, density=True, histtype="step",
               color="C3", lw=1.8, label="halo")
        a.set_xlabel(xlabel)
        a.set_ylabel("normalized")
        a.set_title(ttl)
        a.legend()

    # --- absolute spectra: signal buried under the core (the ~0.5-sigma story) ---
    def abs_spectrum(a, key, scale, bins, xlabel, ttl):
        pc, _ = np.histogram(core[key] * scale, bins=bins, density=True)
        ph, _ = np.histogram(halo[key] * scale, bins=bins, density=True)
        w = np.diff(bins)
        ctr = 0.5 * (bins[1:] + bins[:-1])
        a.step(ctr, N_core * pc * w, color="C0", lw=1.8, where="mid", label="core")
        a.step(ctr, (N_core * pc + N_halo * ph) * w, color="C3", lw=1.3, where="mid",
               label="core + halo")
        a.set_yscale("log")
        a.set_xlabel(xlabel)
        a.set_ylabel("muons / bin  (100 m$^2\\times$10 yr)")
        a.set_title(ttl)
        a.legend()

    abs_spectrum(ax[0, 2], "theta", 1e3, thb, r"arrival angle $\theta$ [mrad]",
                 r"absolute angle: signal excess tiny (S/B$\approx$1.5$\times10^{-3}$)")
    abs_spectrum(ax[1, 2], "r", 1.0, rb, "transverse offset $r$ [m]",
                 "absolute radial: signal buried under the core")

    fig.suptitle(r"Detector-plane muon hits, $E_\mu=%.0f$ GeV, L=%.0f km "
                 r"(MCS + nominal LArTPC);  $N_{halo}\approx%.0f$, $N_{core}\approx%.0f$ "
                 r"(100 m$^2\times$10 yr)"
                 % (Emuon, baseline / 1e3, N_halo, N_core))
    fig.tight_layout()
    fig.savefig(out + ".png", dpi=130)
    fig.savefig(out + ".pdf")
    print("saved %s.png / .pdf  (core N=%d, halo N=%d)" % (out, cx.size, hx.size))


if __name__ == "__main__":
    make_plots()
