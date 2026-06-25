import numpy as np

from simulate import simulate_tau_muons, simulate_numu_cc_muons

# Numeric comparison of the tau-decay "halo" against the direct-nu_mu CC "core" at the
# detector plane (the make-or-break separability question, docs/00-physics.md). NO PLOTS
# yet -- this delivers the numbers (rate ratio, angular-shape separation, and a
# per-exposure Asimov significance) so they can be validated before any figure is made.
#
# Normalization is consistent-relative (docs decision): the core:halo rate ratio is fixed
# by w_integral * BR * acceptance per channel; the common factor (Nmuon, target density,
# detector area, exposure) cancels in the ratio and only multiplies the absolute event
# count. So shape metrics below are absolute-scale-independent, and the Asimov q0 is
# reported PER UNIT exposure K (absolute significance = sqrt(q0_perK * K), K deferred).


def _channel_weight(out):
    """Relative channel rate ~ integral(flux*P*sigma) * BR * geometric acceptance,
    times the detector footprint/threshold survival fraction if a detector was applied."""
    return out["w_integral"] * out["BR"] * out["accept_frac"] * out.get("det_frac", 1.0)


def _ks_1d(a, b):
    """Two-sample Kolmogorov-Smirnov distance (max |CDF_a - CDF_b|), pure numpy."""
    grid = np.sort(np.concatenate([a, b]))
    ca = np.searchsorted(np.sort(a), grid, side="right") / a.size
    cb = np.searchsorted(np.sort(b), grid, side="right") / b.size
    return float(np.max(np.abs(ca - cb)))


def _norm_hist2d(r, theta, r_edges, th_edges):
    """2D histogram of (theta, r) normalized to unit sum (a shape pdf over bins)."""
    H, _, _ = np.histogram2d(theta, r, bins=[th_edges, r_edges])
    s = H.sum()
    return H / s if s > 0 else H


def run_comparison(Emuon=25.0, baseline=1300e3, nusign=1, deltaCP=0.0,
                   n_mc=2_000_000, seed=0, n_rbins=40, n_thbins=40, verbose=True,
                   detector=None):
    """Run halo + core (osc) + core (no-osc) and report separability numbers.

    detector: optional dict of detector.apply_detector kwargs (footprint/smearing/
    threshold) applied to every sample. None -> truth-level detection plane (default).
    """
    halo = simulate_tau_muons(Emuon=Emuon, baseline=baseline, deltaCP=deltaCP,
                              nusign=nusign, n_mc=n_mc, seed=seed)
    core = simulate_numu_cc_muons(Emuon=Emuon, baseline=baseline, deltaCP=deltaCP,
                                  nusign=nusign, n_mc=n_mc, seed=seed + 1, osc_on=True)
    core0 = simulate_numu_cc_muons(Emuon=Emuon, baseline=baseline, deltaCP=deltaCP,
                                   nusign=nusign, n_mc=n_mc, seed=seed + 2, osc_on=False)

    if detector is not None:
        from detector import apply_detector
        halo = apply_detector(halo, seed=seed + 10, **detector)
        core = apply_detector(core, seed=seed + 11, **detector)
        core0 = apply_detector(core0, seed=seed + 12, **detector)

    w_halo = _channel_weight(halo)
    w_core = _channel_weight(core)
    w_core0 = _channel_weight(core0)
    sb_ratio = w_halo / w_core  # tau halo rate / direct-mu core rate

    # --- angular-shape separation (absolute-scale independent) -----------------------
    th_h, th_c = halo["theta"], core["theta"]
    ks_theta = _ks_1d(th_h, th_c)
    core_p90 = np.percentile(th_c, 90)
    halo_p90 = np.percentile(th_h, 90)
    frac_halo_beyond_core90 = float((th_h > core_p90).mean())
    frac_core_beyond_halo90 = float((th_c > halo_p90).mean())

    # --- per-exposure Asimov significance over the 2D (r, theta) plane ----------------
    # Common bin edges spanning both samples; combine with relative channel weights.
    r_all = np.concatenate([halo["r"], core["r"]])
    th_all = np.concatenate([th_h, th_c])
    r_edges = np.linspace(0.0, np.percentile(r_all, 99.5), n_rbins + 1)
    th_edges = np.linspace(0.0, np.percentile(th_all, 99.5), n_thbins + 1)

    p_halo = _norm_hist2d(halo["r"], th_h, r_edges, th_edges)
    p_core = _norm_hist2d(core["r"], th_c, r_edges, th_edges)
    p_core0 = _norm_hist2d(core0["r"], core0["theta"], r_edges, th_edges)

    # Expected bin contents per unit exposure K=1, relative weights:
    #   H1 (appearance): bright-ish core(osc) + halo;   H0 (no appearance): brighter core
    n1 = w_core * p_core + w_halo * p_halo
    n0 = w_core0 * p_core0
    mask = (n1 > 0) & (n0 > 0)
    # Asimov discovery of the appearance signal (data = H1, null = H0), per unit K.
    q0_perK = 2.0 * np.sum(n1[mask] * np.log(n1[mask] / n0[mask]) - (n1[mask] - n0[mask]))

    summary = {
        "sb_ratio": sb_ratio,
        "w_halo": w_halo, "w_core": w_core, "w_core0": w_core0,
        "halo_median_theta_mrad": float(np.median(th_h) * 1e3),
        "core_median_theta_mrad": float(np.median(th_c) * 1e3),
        "halo_p90_theta_mrad": float(halo_p90 * 1e3),
        "core_p90_theta_mrad": float(core_p90 * 1e3),
        "ks_theta": ks_theta,
        "frac_halo_beyond_core90": frac_halo_beyond_core90,
        "frac_core_beyond_halo90": frac_core_beyond_halo90,
        "q0_perK": float(q0_perK),
        "halo_accept": halo["accept_frac"], "core_accept": core["accept_frac"],
    }

    if verbose:
        print("== core vs halo @ Emuon=%.0f GeV, L=%.0f km, nusign=%+d ==" %
              (Emuon, baseline / 1e3, nusign))
        print("  rate:   W_halo/W_core (S/B) = %.3e   core osc/no-osc = %.3f" %
              (sb_ratio, w_core / w_core0))
        print("  theta:  halo median=%.1f (90%%=%.1f)   core median=%.1f (90%%=%.1f) mrad" %
              (summary["halo_median_theta_mrad"], summary["halo_p90_theta_mrad"],
               summary["core_median_theta_mrad"], summary["core_p90_theta_mrad"]))
        print("  shape:  KS(theta)=%.3f   halo beyond core-90%%=%.3f   core beyond halo-90%%=%.3f" %
              (ks_theta, frac_halo_beyond_core90, frac_core_beyond_halo90))
        print("  Asimov: q0 per unit exposure K = %.3e  ->  Z = sqrt(%.3e * K)" %
              (q0_perK, q0_perK))
        print("          (absolute Z needs K = Nmuon * n_target * area * exposure; deferred)")
    return summary


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Numeric core-vs-halo (tau-appearance) comparison at a given "
                    "beam energy and baseline.")
    p.add_argument("--emuon", type=float, default=25.0,
                   help="stored muon (beam) energy [GeV] (default 25)")
    p.add_argument("--baseline", type=float, default=1300.0,
                   help="baseline L [km] (default 1300)")
    p.add_argument("--nusign", type=int, default=1, choices=(1, -1),
                   help="+1 = neutrino (mu- beam), -1 = antineutrino (mu+ beam)")
    p.add_argument("--deltacp", type=float, default=0.0, help="delta_CP [rad]")
    p.add_argument("--n-mc", type=int, default=2_000_000, dest="n_mc",
                   help="Monte Carlo events per channel")
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()

    run_comparison(Emuon=a.emuon, baseline=a.baseline * 1e3, nusign=a.nusign,
                   deltaCP=a.deltacp, n_mc=a.n_mc, seed=a.seed)
