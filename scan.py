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

E_GRID = np.round(np.logspace(np.log10(25.0), np.log10(200.0), 6), 0)  # GeV
L_GRID = np.array([500, 1000, 1300, 2000, 3000, 5000], dtype=float)    # km
_DET = dict(sigma_theta=0.003, sigma_E_frac=0.15, E_threshold=1.0)

# Fermilab (Batavia, IL) -> Seattle, WA straight-line (beam chord) distance [km]
L_FNAL_SEATTLE = 2723.0


def _edges(v, log=False):
    """Cell edges bracketing grid-point centers (linear or log), for pcolormesh."""
    w = np.log(np.asarray(v, float)) if log else np.asarray(v, float)
    mid = 0.5 * (w[1:] + w[:-1])
    e = np.concatenate([[w[0] - (mid[0] - w[0])], mid, [w[-1] + (w[-1] - mid[-1])]])
    return np.exp(e) if log else e


def run_scan(n_mc=800_000, seed=0, area_m2=100.0, years=10.0):
    Z = np.zeros((E_GRID.size, L_GRID.size))       # Z(r, theta, E)
    Zrt = np.zeros_like(Z)                          # Z(r, theta)
    Nhalo = np.zeros_like(Z)
    for i, E in enumerate(E_GRID):
        for j, L in enumerate(L_GRID):
            s = run_comparison(Emuon=float(E), baseline=float(L) * 1e3, n_mc=n_mc,
                               seed=seed, mcs=True, detector=_DET, variant="accelerator",
                               det_area_m2=area_m2, years=years, verbose=False)
            Z[i, j], Zrt[i, j], Nhalo[i, j] = s["Z"], s["Z_rtheta"], s["N_halo"]
            print("E=%3.0f GeV  L=%4.0f km   Z(r,th,E)=%.2f   Z(r,th)=%.2f   N_halo=%.0f"
                  % (E, L, s["Z"], s["Z_rtheta"], s["N_halo"]), flush=True)
    return dict(E=E_GRID, L=L_GRID, Z=Z, Zrt=Zrt, Nhalo=Nhalo,
                area_m2=area_m2, years=years)


def plot_scan(res, out="scan_LE"):
    Z, E, L = res["Z"], res["E"], res["L"]
    fig, ax = plt.subplots(figsize=(9, 6.5))
    im = ax.pcolormesh(_edges(L), _edges(E, log=True), Z, cmap="viridis", shading="flat")
    ax.set_yscale("log")
    ax.set_xlabel("baseline L [km]")
    ax.set_ylabel(r"beam energy $E_\mu$ [GeV]")
    ax.set_yticks(E)
    ax.set_yticklabels([("%.0f" % y) for y in E])
    ax.set_title(r"halo-detection $Z(r,\theta,E)$  (accelerator, %.0f m$^2\times$%.0f yr, "
                 r"MCS + LArTPC)" % (res["area_m2"], res["years"]))
    for i, e in enumerate(E):
        for j, l in enumerate(L):
            ax.text(l, e, "%.2f" % Z[i, j], ha="center", va="center",
                    color="w" if Z[i, j] < 0.7 * Z.max() else "k", fontsize=9)
    # realistic CONUS baseline
    ax.axvline(L_FNAL_SEATTLE, color="C3", ls="--", lw=1.8)
    ax.text(L_FNAL_SEATTLE * 0.98, np.sqrt(E[0] * E[-1]),
            "Fermilab$\\to$Seattle (%.0f km)" % L_FNAL_SEATTLE, color="C3",
            rotation=90, ha="right", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, label=r"$Z$ [$\sigma$]")
    fig.tight_layout()
    fig.savefig(out + ".png", dpi=130)
    fig.savefig(out + ".pdf")
    ib, jb = np.unravel_index(np.argmax(Z), Z.shape)
    print("best: E=%.0f GeV, L=%.0f km -> Z=%.2f sigma (100 m^2 x 10 yr); saved %s.png"
          % (E[ib], L[jb], Z[ib, jb], out))


if __name__ == "__main__":
    res = run_scan()
    np.savez("scan_LE.npz", **res)
    plot_scan(res)
