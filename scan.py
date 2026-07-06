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


def make_grids(nE=18, nL=16, E_lo=25.0, E_hi=200.0, L_lo=500.0, L_hi=5000.0):
    return np.logspace(np.log10(E_lo), np.log10(E_hi), nE), np.linspace(L_lo, L_hi, nL)


def run_scan(nE=18, nL=16, n_mc=800_000, seed=0, area_m2=100.0, years=10.0):
    E_GRID, L_GRID = make_grids(nE, nL)
    Z = np.zeros((E_GRID.size, L_GRID.size))       # Z(r, theta, E)
    Zrt = np.zeros_like(Z)                          # Z(r, theta)
    Nhalo = np.zeros_like(Z)
    for i, E in enumerate(E_GRID):
        for j, L in enumerate(L_GRID):
            s = run_comparison(Emuon=float(E), baseline=float(L) * 1e3, n_mc=n_mc,
                               seed=seed, mcs=True, detector=_DET, variant="accelerator",
                               det_area_m2=area_m2, years=years, verbose=False)
            Z[i, j], Zrt[i, j], Nhalo[i, j] = s["Z"], s["Z_rtheta"], s["N_halo"]
        print("row E=%3.0f GeV done (Z range %.2f-%.2f)" % (E, Z[i].min(), Z[i].max()),
              flush=True)
    return dict(E=E_GRID, L=L_GRID, Z=Z, Zrt=Zrt, Nhalo=Nhalo,
                area_m2=area_m2, years=years)


def plot_scan(res, out="scan_LE", key="Z"):
    Z, E, L = res[key], res["E"], res["L"]
    fig, ax = plt.subplots(figsize=(9, 6.5))

    levels = np.linspace(0.0, np.ceil(Z.max() * 2) / 2, 41)
    cf = ax.contourf(L, E, Z, levels=levels, cmap="viridis", extend="max")
    cs = ax.contour(L, E, Z, levels=[1, 2, 3, 4, 5], colors="w", linewidths=0.8, alpha=0.7)
    ax.clabel(cs, fmt=r"%.0f$\sigma$", fontsize=8)
    ax.set_yscale("log")
    ax.set_yticks([25, 40, 60, 100, 150, 200])
    ax.set_yticklabels(["25", "40", "60", "100", "150", "200"])
    ax.set_xlabel("baseline L [km]")
    ax.set_ylabel(r"beam energy $E_\mu$ [GeV]")
    ax.set_title(r"halo-detection $Z(r,\theta,E)$  (accelerator, %.0f m$^2\times$%.0f yr, "
                 r"MCS + LArTPC)" % (res["area_m2"], res["years"]))

    ax.axvline(L_FNAL_SEATTLE, color="C3", ls="--", lw=1.8)
    ax.text(L_FNAL_SEATTLE * 0.985, np.sqrt(E[0] * E[-1]),
            "Fermilab$\\to$Seattle (%.0f km)" % L_FNAL_SEATTLE, color="C3",
            rotation=90, ha="right", va="center", fontsize=9)

    ib, jb = np.unravel_index(np.argmax(Z), Z.shape)
    ax.plot(L[jb], E[ib], "*", color="w", ms=16, mec="k")
    ax.text(L[jb], E[ib] * 0.9, r"  best %.1f$\sigma$" % Z[ib, jb], color="w",
            ha="center", va="top", fontsize=9)

    fig.colorbar(cf, ax=ax, label=r"$Z$ [$\sigma$]")
    fig.tight_layout()
    fig.savefig(out + ".png", dpi=140)
    fig.savefig(out + ".pdf")
    print("best: E=%.0f GeV, L=%.0f km -> Z=%.2f sigma; saved %s.png"
          % (E[ib], L[jb], Z[ib, jb], out))


if __name__ == "__main__":
    res = run_scan()
    np.savez("scan_LE.npz", **res)
    plot_scan(res)
