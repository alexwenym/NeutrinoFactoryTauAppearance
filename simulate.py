import numpy as np

from flux import numu_flux
from cross_sections import totalXS
from osc import oscillate
from tau_decay import decay_to_muon, BR_TAU_TO_MU

# Monte Carlo of nu_mu -> nu_tau appearance muons reaching a far detector plane,
# built on the existing flux / oscillation / cross-section code.
# See docs/06-simulation-design.md for the physics and v1 assumptions.

_XS = totalXS("resources/nutau_xs.txt", "resources/nutaubar_xs.txt")


def muon_range(E_mu):
    """CSDA muon range in standard rock [m] for E_mu [GeV], dE/dx = a + b E."""
    E = np.asarray(E_mu, dtype=float) * 1e3  # MeV
    a, b, rho = 2.0, 4.0e-6, 2.65  # MeV/(g/cm^2), 1/(g/cm^2), g/cm^3
    X = (1.0 / b) * np.log1p((b / a) * E)  # g/cm^2
    return X / rho / 100.0  # m


def simulate_tau_muons(Emuon=25.0, baseline=1300e3, deltaCP=0.0, nusign=1,
                       Pmuon=0.0, n_mc=2_000_000, n_grid=200, seed=0):
    """
    Generate tau-decay muons arriving at the detector plane.

    Returns a dict of per-muon arrays (arrival angle theta [rad], transverse
    offset r [m] from the beam axis, energies, azimuth) for accepted muons, plus
    the production-region size dmax [m] and the geometric acceptance fraction.
    """
    rng = np.random.default_rng(seed)

    # Energy grid: P(nu_mu->nu_tau) at the full baseline, flux shape, and sigma.
    Egrid = np.linspace(0.2, Emuon, n_grid)
    amp = oscillate(Egrid, baseline, 1, nusign, deltaCP=deltaCP)  # (4, n_grid)
    P_grid = np.abs(amp[2]) ** 2                                  # nu_tau component
    flux_grid = numu_flux(Emuon, 1.0, Pmuon, Egrid, 1.0, baseline)  # relative (Nmu=1)
    xs_grid = _XS.sigma(Egrid, is_nubar=(nusign < 0))
    w_grid = flux_grid * P_grid * xs_grid  # relative tau production rate vs energy

    # Importance-sample nu energy proportional to the production rate.
    cdf = np.cumsum(w_grid)
    cdf /= cdf[-1]
    E_nu = np.interp(rng.uniform(0.0, 1.0, n_mc), cdf, Egrid)

    # v0: tau collinear with nu, E_tau = E_nu. Decay to a muon.
    E_tau = E_nu
    E_mu, theta, phi = decay_to_muon(E_tau, rng)

    # Production point: uniform along the beam, bounded by the max muon range.
    dmax = float(muon_range(Emuon))
    d = rng.uniform(0.0, dmax, n_mc)

    # Geometry + range acceptance: muon must travel its path length d/cos(theta).
    forward = np.cos(theta) > 0.0
    path = np.where(forward, d / np.where(forward, np.cos(theta), 1.0), np.inf)
    accept = forward & (muon_range(E_mu) >= path)

    r = d * np.tan(theta)  # transverse offset at the detector plane [m]

    return {
        "E_nu": E_nu[accept],
        "E_mu": E_mu[accept],
        "theta": theta[accept],
        "r": r[accept],
        "phi": phi[accept],
        "dmax": dmax,
        "accept_frac": float(accept.mean()),
        "BR": BR_TAU_TO_MU,
    }


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("figures.mplstyle")

    out = simulate_tau_muons(Emuon=25.0, baseline=1300e3)
    theta_mrad = out["theta"] * 1e3  # rad -> mrad
    print("accepted muons: %d  (acceptance %.3f, production region %.1f m)"
          % (out["theta"].size, out["accept_frac"], out["dmax"]))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].hist(theta_mrad, bins=80, histtype="step", color="C0")
    axes[0].set_xlabel(r"arrival angle $\theta$ from beam axis [mrad]")
    axes[0].set_ylabel("muons / bin")
    axes[1].hist(out["r"], bins=80, histtype="step", color="C1")
    axes[1].set_xlabel(r"transverse offset $r$ at detector plane [m]")
    axes[1].set_ylabel("muons / bin")
    fig.suptitle(r"$\nu_\mu\to\nu_\tau$ appearance muons ($E_\mu=25$ GeV, L=1300 km)")
    fig.tight_layout()
    fig.savefig("tau_muons_detector_plane.pdf")
    print("saved tau_muons_detector_plane.pdf")
