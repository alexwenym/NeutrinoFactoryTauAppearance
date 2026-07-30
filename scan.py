import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare import run_comparison

# Scan the halo-detection significance over beam energy (log 25-200 GeV, within the GENIE
# table range) and baseline. Everything else fixed: accelerator beam, nominal LArTPC, MCS,
# 100 m^2 x 10 yr. Z scales as sqrt(area*years), so the map shape is the physics; the
# absolute scale is for this nominal exposure. The beam intensity's ~1/E_muon fall-off is
# included (n_muon_decays), so the E axis is a fair comparison.

_DET = dict(sigma_theta=0.003, sigma_E_frac=0.15, E_threshold=1.0)

# Fermilab (Batavia, IL) -> Seattle, WA straight-line (beam chord) distance [km]
L_FNAL_SEATTLE = 2723.0
# Earth mean diameter [km]
L_EARTH_DIAM = 12742.0


def make_grids(nE=18, nL=22, E_lo=25.0, E_hi=200.0, L_lo=500.0, L_hi=12500.0):
    return np.logspace(np.log10(E_lo), np.log10(E_hi), nE), np.linspace(L_lo, L_hi, nL)


def run_scan(nE=18, nL=22, n_mc=800_000, seed=0, area_m2=100.0, years=10.0, nusign=1):
    E_GRID, L_GRID = make_grids(nE, nL)
    Z = np.zeros((E_GRID.size, L_GRID.size))       # Z(r, theta, E)
    Zrt = np.zeros_like(Z)                          # Z(r, theta)
    Nhalo = np.zeros_like(Z)
    label = "nu" if nusign == 1 else "nubar"
    for i, E in enumerate(E_GRID):
        for j, L in enumerate(L_GRID):
            s = run_comparison(Emuon=float(E), baseline=float(L) * 1e3, n_mc=n_mc,
                               seed=seed, mcs=True, detector=_DET, variant="accelerator",
                               det_area_m2=area_m2, years=years, verbose=False,
                               nusign=nusign)
            Z[i, j], Zrt[i, j], Nhalo[i, j] = s["Z"], s["Z_rtheta"], s["N_halo"]
        print("%s row E=%3.0f GeV done (Z range %.2f-%.2f)" % (
            label, E, Z[i].min(), Z[i].max()), flush=True)
    return dict(E=E_GRID, L=L_GRID, Z=Z, Zrt=Zrt, Nhalo=Nhalo,
                area_m2=area_m2, years=years, nusign=nusign)


def run_scan_channels(nE=18, nL=22, n_mc=800_000, seed=0, area_m2=100.0, years=10.0,
                      nusign=1):
    """Run (L,E) scan returning Z for numu-only, nue-only, and combined channels."""
    E_GRID, L_GRID = make_grids(nE, nL)
    Z_numu     = np.zeros((E_GRID.size, L_GRID.size))
    Z_nue      = np.zeros_like(Z_numu)
    Z_combined = np.zeros_like(Z_numu)
    Nhalo_numu = np.zeros_like(Z_numu)
    Nhalo_nue  = np.zeros_like(Z_numu)
    label = "nu" if nusign == 1 else "nubar"
    for i, E in enumerate(E_GRID):
        for j, L in enumerate(L_GRID):
            s = run_comparison(Emuon=float(E), baseline=float(L) * 1e3, n_mc=n_mc,
                               seed=seed, mcs=True, detector=_DET, variant="accelerator",
                               det_area_m2=area_m2, years=years, verbose=False,
                               nusign=nusign, include_nue=True)
            Z_numu[i, j]     = s["Z"]
            Z_nue[i, j]      = s["Z_nue"]
            Z_combined[i, j] = s["Z_combined"]
            Nhalo_numu[i, j] = s["N_halo"]
            Nhalo_nue[i, j]  = s["N_halo_nue"]
        print("%s channels row E=%3.0f GeV done  numu-Z=%.2f-%.2f  nue-Z=%.2f-%.2f  comb-Z=%.2f-%.2f" % (
            label, E, Z_numu[i].min(), Z_numu[i].max(),
            Z_nue[i].min(), Z_nue[i].max(),
            Z_combined[i].min(), Z_combined[i].max()), flush=True)
    return dict(E=E_GRID, L=L_GRID,
                Z_numu=Z_numu, Z_nue=Z_nue, Z_combined=Z_combined,
                Nhalo_numu=Nhalo_numu, Nhalo_nue=Nhalo_nue,
                area_m2=area_m2, years=years, nusign=nusign)


def plot_scan(res, out="scan_LE", key="Z"):
    Z, E, L = res[key], res["E"], res["L"]

    FS = 18   # base font size
    plt.rcParams.update({"font.size": FS, "axes.titlesize": FS,
                         "axes.labelsize": FS, "xtick.labelsize": FS - 2,
                         "ytick.labelsize": FS - 2, "legend.fontsize": FS - 2})

    fig, ax = plt.subplots(figsize=(14, 8))

    levels = np.linspace(0.0, np.ceil(Z.max() * 2) / 2, 41)
    cf = ax.contourf(L, E, Z, levels=levels, cmap="viridis", extend="max")
    sigma_levels = list(range(1, int(np.floor(Z.max())) + 1))
    cs = ax.contour(L, E, Z, levels=sigma_levels, colors="w", linewidths=1.4, alpha=0.8)
    ax.clabel(cs, fmt=r"%.0f$\sigma$", fontsize=FS - 2)
    ax.set_yscale("log")
    ax.set_yticks([25, 40, 60, 100, 150, 200])
    ax.set_yticklabels(["25", "40", "60", "100", "150", "200"])
    ax.set_xlabel("Baseline $L$ [km]")
    ax.set_ylabel(r"Beam energy $E_\mu$ [GeV]")
    ax.set_title(r"Halo-detection $Z(r,\theta,E)$  (accelerator, %.0f m$^2\times$%.0f yr, "
                 r"MCS + LArTPC)" % (res["area_m2"], res["years"]), pad=10)

    ax.axvline(L_FNAL_SEATTLE, color="C3", ls="--", lw=2.0)
    ax.text(L_FNAL_SEATTLE * 0.982, np.sqrt(E[0] * E[-1]),
            "FNAL$\\to$Seattle\n(%.0f km)" % L_FNAL_SEATTLE, color="C3",
            rotation=90, ha="right", va="center", fontsize=FS - 4)

    ax.axvline(L_EARTH_DIAM, color="w", ls=":", lw=2.0, alpha=0.85)
    ax.set_xlim(L[0], L[-1])

    ib, jb = np.unravel_index(np.argmax(Z), Z.shape)
    ax.plot(L[jb], E[ib], "*", color="w", ms=20, mec="k", mew=1.2)
    ax.text(L[jb], E[ib] * 0.87, r"  best %.1f$\sigma$" % Z[ib, jb], color="w",
            ha="center", va="top", fontsize=FS - 3)

    cb = fig.colorbar(cf, ax=ax, label=r"$Z$ [$\sigma$]")
    cb.ax.tick_params(labelsize=FS - 2)
    cb.set_label(r"$Z$ [$\sigma$]", fontsize=FS)
    fig.tight_layout()
    fig.savefig(out + ".png", dpi=140)
    fig.savefig(out + ".pdf")
    print("best: E=%.0f GeV, L=%.0f km -> Z=%.2f sigma; saved %s.png"
          % (E[ib], L[jb], Z[ib, jb], out))


def plot_scan_both(res_nu, res_nubar, out="scan_LE_charges", key="Z"):
    """Side-by-side heatmaps for nu (mu- beam) and nubar (mu+ beam)."""
    FS = 18
    plt.rcParams.update({"font.size": FS, "axes.titlesize": FS,
                         "axes.labelsize": FS, "xtick.labelsize": FS - 2,
                         "ytick.labelsize": FS - 2, "legend.fontsize": FS - 2})

    Znu,  E, L = res_nu[key],  res_nu["E"],  res_nu["L"]
    Znub, _, _ = res_nubar[key], res_nubar["E"], res_nubar["L"]

    vmax = np.ceil(max(Znu.max(), Znub.max()) * 2) / 2
    levels = np.linspace(0.0, vmax, 41)
    sigma_levels = list(range(1, int(np.floor(max(Znu.max(), Znub.max()))) + 1))

    fig, axes = plt.subplots(1, 2, figsize=(26, 8), sharey=True)
    titles = [r"$\mu^-$ beam: $\nu_\mu \to \nu_\tau$ appearance",
              r"$\mu^+$ beam: $\bar\nu_\mu \to \bar\nu_\tau$ appearance"]

    for ax, Z, title in zip(axes, [Znu, Znub], titles):
        cf = ax.contourf(L, E, Z, levels=levels, cmap="viridis", extend="max")
        cs = ax.contour(L, E, Z, levels=sigma_levels, colors="w", linewidths=1.2, alpha=0.8)
        ax.clabel(cs, fmt=r"%.0f$\sigma$", fontsize=FS - 3)
        ax.set_yscale("log")
        ax.set_yticks([25, 40, 60, 100, 150, 200])
        ax.set_yticklabels(["25", "40", "60", "100", "150", "200"])
        ax.set_xlabel("Baseline $L$ [km]")
        ax.set_title(title, pad=10)
        ax.set_xlim(L[0], L[-1])

        ax.axvline(L_FNAL_SEATTLE, color="C3", ls="--", lw=2.0)
        ax.text(L_FNAL_SEATTLE * 0.982, np.sqrt(E[0] * E[-1]),
                "FNAL$\\to$Seattle\n(%.0f km)" % L_FNAL_SEATTLE, color="C3",
                rotation=90, ha="right", va="center", fontsize=FS - 5)
        ax.axvline(L_EARTH_DIAM, color="w", ls=":", lw=1.8, alpha=0.8)

        ib, jb = np.unravel_index(np.argmax(Z), Z.shape)
        ax.plot(L[jb], E[ib], "*", color="w", ms=18, mec="k", mew=1.2)
        ax.text(L[jb], E[ib] * 0.87, r"  best %.1f$\sigma$" % Z[ib, jb],
                color="w", ha="center", va="top", fontsize=FS - 4)

        cb = fig.colorbar(cf, ax=ax, label=r"$Z$ [$\sigma$]")
        cb.ax.tick_params(labelsize=FS - 2)
        cb.set_label(r"$Z$ [$\sigma$]", fontsize=FS)

    axes[0].set_ylabel(r"Beam energy $E_\mu$ [GeV]")
    fig.suptitle(r"Halo-detection $Z(r,\theta,E)$  —  accelerator beam, "
                 r"%.0f m$^2\times$%.0f yr, MCS + LArTPC" % (
                     res_nu["area_m2"], res_nu["years"]), fontsize=FS + 1, y=1.01)
    fig.tight_layout()
    fig.savefig(out + ".png", dpi=140, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print("saved %s.png" % out)


def plot_scan_channels(res, out="scan_LE_channels"):
    """Three-panel (L,E) heatmap: nu_mu only | nu_e only | combined."""
    E, L = res["E"], res["L"]
    Z_sets = [
        (res["Z_numu"],     r"$\nu_\mu \to \nu_\tau$  (halo)",     "numu"),
        (res["Z_nue"],      r"$\bar\nu_e \to \bar\nu_\tau$  (nue)", "nue"),
        (res["Z_combined"], r"Combined  $\nu_\mu + \bar\nu_e$",     "comb"),
    ]
    if res.get("nusign", 1) == -1:
        Z_sets = [
            (res["Z_numu"],     r"$\bar\nu_\mu \to \bar\nu_\tau$  (halo)", "numu"),
            (res["Z_nue"],      r"$\nu_e \to \nu_\tau$  (nue)",            "nue"),
            (res["Z_combined"], r"Combined  $\bar\nu_\mu + \nu_e$",        "comb"),
        ]

    FS = 18
    plt.rcParams.update({"font.size": FS, "axes.titlesize": FS,
                         "axes.labelsize": FS, "xtick.labelsize": FS - 2,
                         "ytick.labelsize": FS - 2, "legend.fontsize": FS - 2})

    vmax = np.ceil(max(Z.max() for Z, _, _ in Z_sets) * 2) / 2
    levels = np.linspace(0.0, vmax, 41)

    fig, axes = plt.subplots(1, 3, figsize=(36, 8), sharey=True)
    for ax, (Z, title, _) in zip(axes, Z_sets):
        cf = ax.contourf(L, E, Z, levels=levels, cmap="viridis", extend="max")
        sigma_levels = list(range(1, int(np.floor(Z.max())) + 1))
        if sigma_levels:
            cs = ax.contour(L, E, Z, levels=sigma_levels, colors="w", linewidths=1.2, alpha=0.8)
            ax.clabel(cs, fmt=r"%.0f$\sigma$", fontsize=FS - 3)
        ax.set_yscale("log")
        ax.set_yticks([25, 40, 60, 100, 150, 200])
        ax.set_yticklabels(["25", "40", "60", "100", "150", "200"])
        ax.set_xlabel("Baseline $L$ [km]")
        ax.set_title(title, pad=10)
        ax.set_xlim(L[0], L[-1])
        ax.axvline(L_FNAL_SEATTLE, color="C3", ls="--", lw=2.0)
        ax.text(L_FNAL_SEATTLE * 0.982, np.sqrt(E[0] * E[-1]),
                "FNAL$\\to$Seattle", color="C3", rotation=90,
                ha="right", va="center", fontsize=FS - 5)
        ax.axvline(L_EARTH_DIAM, color="w", ls=":", lw=1.8, alpha=0.8)
        ib, jb = np.unravel_index(np.argmax(Z), Z.shape)
        ax.plot(L[jb], E[ib], "*", color="w", ms=18, mec="k", mew=1.2)
        ax.text(L[jb], E[ib] * 0.87, r"  best %.1f$\sigma$" % Z[ib, jb],
                color="w", ha="center", va="top", fontsize=FS - 4)
        cb = fig.colorbar(cf, ax=ax, label=r"$Z$ [$\sigma$]")
        cb.ax.tick_params(labelsize=FS - 2)
        cb.set_label(r"$Z$ [$\sigma$]", fontsize=FS)

    axes[0].set_ylabel(r"Beam energy $E_\mu$ [GeV]")
    beam_label = r"$\mu^-$ beam" if res.get("nusign", 1) == 1 else r"$\mu^+$ beam"
    fig.suptitle(r"Halo-detection $Z(r,\theta,E)$ by channel  —  %s, "
                 r"%.0f m$^2\times$%.0f yr, MCS + LArTPC" % (
                     beam_label, res["area_m2"], res["years"]),
                 fontsize=FS + 1, y=1.01)
    fig.tight_layout()
    fig.savefig(out + ".png", dpi=140, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print("saved %s.png" % out)
    for Z, _, tag in Z_sets:
        ib, jb = np.unravel_index(np.argmax(Z), Z.shape)
        print("  %s: best %.2f sigma at E=%.0f GeV, L=%.0f km" % (
            tag, Z[ib, jb], E[ib], L[jb]))


if __name__ == "__main__":
    import sys
    if "--channels" in sys.argv:
        print("=== channels scan (mu- beam, numu + nue + combined) ===")
        res = run_scan_channels(nusign=1)
        np.savez("scan_LE_channels.npz", **res)
        plot_scan_channels(res)
    elif "--charges" in sys.argv:
        print("=== nu scan (mu- beam) ===")
        res_nu = run_scan(nusign=1)
        np.savez("scan_LE_nu.npz", **res_nu)
        print("=== nubar scan (mu+ beam) ===")
        res_nubar = run_scan(nusign=-1)
        np.savez("scan_LE_nubar.npz", **res_nubar)
        plot_scan_both(res_nu, res_nubar)
        plot_scan(res_nu,   out="scan_LE_nu")
        plot_scan(res_nubar, out="scan_LE_nubar")
    else:
        res = run_scan()
        np.savez("scan_LE.npz", **res)
        plot_scan(res)
