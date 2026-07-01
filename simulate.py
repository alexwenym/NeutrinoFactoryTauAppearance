import numpy as np

from flux import numu_flux
from osc import oscillate
from tau_decay import decay_to_muon, BR_TAU_TO_MU
import dis_kinematics as dis

# Monte Carlo of nu_mu -> nu_tau appearance muons reaching a far detector plane,
# built on the existing flux / oscillation code plus GENIE DIS kinematics+cross sections
# (dis_kinematics). See docs/06-simulation-design.md for the physics and v1 assumptions.


def muon_range(E_mu):
    """CSDA muon range in standard rock [m] for E_mu [GeV], dE/dx = a + b E."""
    E = np.asarray(E_mu, dtype=float) * 1e3  # MeV
    a, b, rho = 2.0, 4.0e-6, 2.65  # MeV/(g/cm^2), 1/(g/cm^2), g/cm^3
    X = (1.0 / b) * np.log1p((b / a) * E)  # g/cm^2
    return X / rho / 100.0  # m


def _compose_directions(theta_prod, phi_prod, theta_rel, phi_rel):
    """Compose a decay direction (theta_rel, phi_rel measured w.r.t. the parent's flight
    direction) with the parent's production direction (theta_prod, phi_prod w.r.t. the
    beam +z), returning the daughter's (theta, phi) w.r.t. the beam.

    Builds an orthonormal frame (e1, e2, parent_hat) and expresses the daughter unit
    vector in the lab. For a parent along +z this reduces to (theta_rel, phi_rel).
    """
    tx = np.sin(theta_prod) * np.cos(phi_prod)
    ty = np.sin(theta_prod) * np.sin(phi_prod)
    tz = np.cos(theta_prod)

    # e1 = normalize(z_hat x parent) = (-ty, tx, 0)/|...|; fall back to x_hat if parallel
    e1x, e1y = -ty, tx
    n1 = np.hypot(e1x, e1y)
    small = n1 < 1e-12
    inv = np.where(small, 1.0, 1.0 / np.where(small, 1.0, n1))
    e1x = np.where(small, 1.0, e1x * inv)
    e1y = np.where(small, 0.0, e1y * inv)
    e1z = np.zeros_like(tz)
    # e2 = parent x e1
    e2x = ty * e1z - tz * e1y
    e2y = tz * e1x - tx * e1z
    e2z = tx * e1y - ty * e1x

    sR, cR = np.sin(theta_rel), np.cos(theta_rel)
    a, b = sR * np.cos(phi_rel), sR * np.sin(phi_rel)
    mx = a * e1x + b * e2x + cR * tx
    my = a * e1y + b * e2y + cR * ty
    mz = a * e1z + b * e2z + cR * tz

    theta = np.arccos(np.clip(mz, -1.0, 1.0))
    phi = np.arctan2(my, mx)
    return theta, phi


# multiple Coulomb scattering in standard rock (Highland/Moliere, Gaussian core)
_X0_ROCK_CM = 10.02     # radiation length of standard rock [cm] (26.54 g/cm^2 / 2.65)
_DEDX_ROCK = 5.3e-3     # ionization loss [GeV/cm] (2.0 MeV/(g/cm^2) * 2.65 g/cm^3)


def _mcs_sigma_theta_plane(E_mu, path):
    """Per-projection-plane RMS MCS angle [rad] (Highland), integrated over the rock path
    including ionization energy loss along it. E_mu [GeV] is the production energy, path
    [m] the slant path to the plane.

        sigma^2 = (13.6 MeV)^2 / (X0 dE/dx) * (1/E_arr - 1/E0) * [1 + 0.038 ln(L/X0)]^2

    with E_arr = E0 - dE/dx * L the muon energy at the plane (floored: the Gaussian-core
    approximation breaks down for near-stopping muons, which scatter enormously anyway).
    A full treatment (Moliere tails, stochastic loss/range straggling) is PROPOSAL's job;
    this is the pheno-level analytical approximation.
    """
    L = np.asarray(path, dtype=float) * 100.0          # cm
    E0 = np.asarray(E_mu, dtype=float)                  # GeV
    E_arr = np.maximum(E0 - _DEDX_ROCK * L, 0.5)        # GeV at the plane (floored)
    inv_diff = np.maximum(1.0 / E_arr - 1.0 / E0, 0.0)  # 1/GeV
    var = (0.0136 ** 2) / (_X0_ROCK_CM * _DEDX_ROCK) * inv_diff
    logcorr = np.maximum(1.0 + 0.038 * np.log(np.maximum(L / _X0_ROCK_CM, 1e-3)), 0.0)
    return np.sqrt(var) * logcorr


def _project_to_detector(E_mu, theta, phi, dmax, n_mc, rng, mcs=False):
    """Shared geometry: sample a production point d ~ U(0, dmax) along the beam, require
    the muon to be forward and to have enough range to reach the plane, and project to
    the transverse landing position (arrival angle/azimuth returned too).

    With mcs=True, add Highland multiple-Coulomb scattering in the rock over the path
    d/cos(theta): a correlated angular deflection (RMS sigma_theta per projection plane)
    and lateral displacement (RMS L*sigma_theta/sqrt(3), correlation rho=sqrt(3)/2). Both
    grow with the path length, so muons produced further upstream scatter more.

    Returns (accept, r, theta_out, phi_out, d).
    """
    d = rng.uniform(0.0, dmax, n_mc)
    forward = np.cos(theta) > 0.0
    path = np.where(forward, d / np.where(forward, np.cos(theta), 1.0), np.inf)
    accept = forward & (muon_range(E_mu) >= path)

    # straight-line arrival: transverse angle components and landing position
    ang_x = theta * np.cos(phi)
    ang_y = theta * np.sin(phi)
    pos_x = d * np.tan(theta) * np.cos(phi)
    pos_y = d * np.tan(theta) * np.sin(phi)

    if mcs:
        L_safe = np.where(accept, path, 0.0)             # no MCS for rejected muons
        s_theta = _mcs_sigma_theta_plane(E_mu, L_safe)   # per projection plane [rad]
        s_pos = L_safe * s_theta / np.sqrt(3.0)          # lateral displacement RMS [m]
        rho = np.sqrt(3.0) / 2.0
        z1x, z2x = rng.normal(size=n_mc), rng.normal(size=n_mc)
        z1y, z2y = rng.normal(size=n_mc), rng.normal(size=n_mc)
        ang_x = ang_x + s_theta * z1x
        ang_y = ang_y + s_theta * z1y
        anti = np.sqrt(1.0 - rho ** 2)
        pos_x = pos_x + s_pos * (rho * z1x + anti * z2x)
        pos_y = pos_y + s_pos * (rho * z1y + anti * z2y)

    theta_out = np.hypot(ang_x, ang_y)
    phi_out = np.arctan2(ang_y, ang_x)
    r = np.hypot(pos_x, pos_y)
    return accept, r, theta_out, phi_out, d


def simulate_tau_muons(Emuon=25.0, baseline=1300e3, deltaCP=0.0, nusign=1,
                       Pmuon=0.0, n_mc=2_000_000, n_grid=200, seed=0, mcs=False):
    """
    Generate tau-decay muons arriving at the detector plane.

    Returns a dict of per-muon arrays (arrival angle theta [rad], transverse
    offset r [m] from the beam axis, energies, azimuth) for accepted muons, plus
    the production-region size dmax [m] and the geometric acceptance fraction.
    """
    rng = np.random.default_rng(seed)

    # Energy grid: P(nu_mu->nu_tau) at the full baseline, flux shape, and sigma.
    channel_tau = "nutaubar" if nusign < 0 else "nutau"
    Egrid = np.linspace(0.2, Emuon, n_grid)
    amp = oscillate(Egrid, baseline, 1, nusign, deltaCP=deltaCP)  # (4, n_grid)
    P_grid = np.abs(amp[2]) ** 2                                  # nu_tau component
    flux_grid = numu_flux(Emuon, 1.0, Pmuon, Egrid, 1.0, baseline)  # relative (Nmu=1)
    xs_grid = dis.sigma(Egrid, channel_tau)  # nu_tau CC (GENIE DIS)
    w_grid = flux_grid * P_grid * xs_grid  # relative tau production rate vs energy
    w_integral = float(np.trapz(w_grid, Egrid))

    # Importance-sample nu energy proportional to the production rate.
    cdf = np.cumsum(w_grid)
    cdf /= cdf[-1]
    E_nu = np.interp(rng.uniform(0.0, 1.0, n_mc), cdf, Egrid)

    # tau production via DIS: E_tau = (1-y) E_nu and a production angle off the beam.
    E_tau, theta_tau, phi_tau = dis.sample_dis_lepton(E_nu, rng, channel_tau)

    # tau -> mu decay (muon angle relative to the tau direction), then compose with the
    # tau production direction to get the muon angle relative to the beam.
    E_mu, theta_rel, phi_rel = decay_to_muon(E_tau, rng, nusign=nusign)
    theta, phi = _compose_directions(theta_tau, phi_tau, theta_rel, phi_rel)

    dmax = float(muon_range(Emuon))
    accept, r, theta, phi, d = _project_to_detector(E_mu, theta, phi, dmax, n_mc, rng,
                                                    mcs=mcs)

    return {
        "E_nu": E_nu[accept],
        "E_mu": E_mu[accept],
        "theta": theta[accept],
        "r": r[accept],
        "phi": phi[accept],
        "dmax": dmax,
        "accept_frac": float(accept.mean()),
        "w_integral": w_integral,
        "BR": BR_TAU_TO_MU,
    }


def simulate_numu_cc_muons(Emuon=25.0, baseline=1300e3, deltaCP=0.0, nusign=1,
                           Pmuon=0.0, n_mc=2_000_000, n_grid=200, seed=0,
                           osc_on=True, mcs=False):
    """
    Generate direct nu_mu CC ("core") muons arriving at the detector plane -- the
    background against which the tau-decay halo is compared.

    Same machinery as simulate_tau_muons but: weight by nu_mu survival P(nu_mu->nu_mu)
    (component index 1) instead of nu_tau appearance, use the nu_mu CC cross section,
    and produce a prompt DIS muon (no tau, no decay, no branching ratio). With
    osc_on=False the survival is forced to 1 -- the "no nu_mu->nu_tau oscillation"
    hypothesis (a brighter core, no halo).

    Returns the same dict structure as simulate_tau_muons (BR is 1 here).
    """
    rng = np.random.default_rng(seed)

    Egrid = np.linspace(0.2, Emuon, n_grid)
    if osc_on:
        amp = oscillate(Egrid, baseline, 1, nusign, deltaCP=deltaCP)  # (4, n_grid)
        P_grid = np.abs(amp[1]) ** 2                                  # nu_mu survival
    else:
        P_grid = np.ones_like(Egrid)                                  # no disappearance
    flux_grid = numu_flux(Emuon, 1.0, Pmuon, Egrid, 1.0, baseline)  # relative (Nmu=1)
    channel_mu = "numubar" if nusign < 0 else "numu"
    xs_grid = dis.sigma(Egrid, channel_mu)  # nu_mu CC (from the DIS model/splines)
    w_grid = flux_grid * P_grid * xs_grid
    w_integral = float(np.trapz(w_grid, Egrid))

    cdf = np.cumsum(w_grid)
    cdf /= cdf[-1]
    E_nu = np.interp(rng.uniform(0.0, 1.0, n_mc), cdf, Egrid)

    # prompt muon directly from DIS (angle already w.r.t. the beam).
    E_mu, theta, phi = dis.sample_dis_lepton(E_nu, rng, channel_mu)

    dmax = float(muon_range(Emuon))
    accept, r, theta, phi, d = _project_to_detector(E_mu, theta, phi, dmax, n_mc, rng,
                                                    mcs=mcs)

    return {
        "E_nu": E_nu[accept],
        "E_mu": E_mu[accept],
        "theta": theta[accept],
        "r": r[accept],
        "phi": phi[accept],
        "dmax": dmax,
        "accept_frac": float(accept.mean()),
        "w_integral": w_integral,
        "BR": 1.0,
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
