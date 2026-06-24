import numpy as np

# tau -> mu nu_mu nu_tau leptonic decay kinematics.
#
# v1 approximations (see docs/06-simulation-design.md):
#  - unpolarized tau
#  - massless muon Michel spectrum  dN/dx = 2 x^2 (3 - 2x),  x = 2 E*_mu / m_tau
#  - muon direction isotropic in the tau rest frame
#  - muon treated as massless in the boost (m_mu / m_tau ~ 6%, sub-leading)

M_TAU = 1.77686  # GeV
M_MU = 0.105658  # GeV
BR_TAU_TO_MU = 0.1739  # tau -> mu nu nu branching ratio


def sample_michel_x(size, rng):
    """Sample x = 2 E*_mu / m_tau from the Michel spectrum 2 x^2 (3 - 2x) on [0, 1]."""
    # rejection sampling; pdf peaks at x = 1 with value 2
    x = np.empty(size)
    n_filled = 0
    while n_filled < size:
        n = size - n_filled
        xt = rng.uniform(0.0, 1.0, n)
        yt = rng.uniform(0.0, 2.0, n)
        accept = yt < 2.0 * xt**2 * (3.0 - 2.0 * xt)
        n_acc = np.count_nonzero(accept)
        x[n_filled:n_filled + n_acc] = xt[accept]
        n_filled += n_acc
    return x


def decay_to_muon(E_tau, rng):
    """
    Decay an array of taus (each travelling along +z with energy E_tau [GeV])
    to a muon, returning the lab-frame muon (E_mu [GeV], theta [rad], phi [rad]).

    theta is the muon angle relative to the tau (= beam) direction; phi is the
    azimuth about it.
    """
    E_tau = np.asarray(E_tau, dtype=float)
    n = E_tau.size

    # muon in the tau rest frame (massless)
    x = sample_michel_x(n, rng)
    E_star = x * M_TAU / 2.0
    cos_star = rng.uniform(-1.0, 1.0, n)
    sin_star = np.sqrt(1.0 - cos_star**2)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)

    pz_star = E_star * cos_star
    pT = E_star * sin_star

    # boost along +z into the lab
    gamma = E_tau / M_TAU
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    E_mu = gamma * (E_star + beta * pz_star)
    pz = gamma * (pz_star + beta * E_star)

    theta = np.arctan2(pT, pz)
    return E_mu, theta, phi
