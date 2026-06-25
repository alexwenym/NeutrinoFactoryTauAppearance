import numpy as np
import alouette

# tau -> mu nu_mu nu_tau leptonic decay kinematics, via alouette (a TAUOLA wrapper).
#
# alouette gives the full TAUOLA decay (all channels), so we run taus at rest and
# keep only the decays that contain a muon -- this reproduces the leptonic branching
# ratio (~17.4%) and the full muon kinematics (massive muon, radiative corrections,
# and the tau-spin -> muon-direction correlation) rather than the analytic massless
# Michel spectrum used in v0.
#
# Design (see docs/06-simulation-design.md):
#  - Decays are generated in the tau REST FRAME once into a pool, then sampled and
#    boosted along +z (the beam = tau direction). alouette.decay() is single-threaded
#    (~5e4 decays/s), so a precomputed pool keeps the per-event hot path vectorized
#    and alouette-free. The rest-frame spectrum is independent of E_tau, so one pool
#    serves all tau energies.
#  - tau POLARIZATION is fixed-helicity: -1 for tau^- (V-A: a nu makes a left-handed
#    tau^-), +1 for tau^+. The polarization axis is the beam (+z); decays are
#    azimuthally symmetric about it, so we redraw the muon azimuth per event.
#  - v1 still assumes the tau is collinear with the neutrino (E_tau = E_nu, no tau
#    production p_T); that lives in simulate.py, not here.
#
# Energy-dependent CC polarization (Hagiwara et al.) is a later refinement; set
# polarization=0.0 to recover the unpolarized spectrum (cross-checks Michel <x>=0.70).

M_TAU = 1.77686  # GeV
M_MU = 0.105658  # GeV
BR_TAU_TO_MU = 0.1739  # tau -> mu nu nu branching ratio

# Number of rest-frame tau->mu decays precomputed per (charge, polarization) pool.
# Sampled with replacement and boosted, so it need not match the MC event count.
POOL_SIZE = 100_000

# Cache of rest-frame muon pools, keyed by (nusign, helicity, pool_size).
# Each entry is (pz_star, pT_star, E_star) arrays [GeV] in the tau rest frame.
_POOLS = {}


def sample_michel_x(size, rng):
    """Sample x = 2 E*_mu / m_tau from the massless Michel spectrum 2 x^2 (3 - 2x).

    Kept for cross-checking the unpolarized alouette pool (both give <x> = 0.70).
    """
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


def _build_pool(nusign, helicity, pool_size, rng):
    """Generate `pool_size` rest-frame tau->mu decays with alouette and return the
    muon (pz_star, pT_star, E_star) along/about the polarization axis (+z)."""
    tau_pid = 15 if nusign > 0 else -15
    mu_pid = 13 if nusign > 0 else -13
    polarisation = (0.0, 0.0, float(helicity))

    # alouette has its own PRNG; seed it from the caller's rng for reproducibility.
    alouette.random.set(int(rng.integers(1, 2**31 - 1)))

    pz_star = np.empty(pool_size)
    pT_star = np.empty(pool_size)
    E_star = np.empty(pool_size)
    i = 0
    while i < pool_size:
        prod = alouette.decay(mode=0, pid=tau_pid, momentum=(0.0, 0.0, 0.0),
                              polarisation=polarisation)
        ids = np.asarray(prod.pid)
        sel = np.where(ids == mu_pid)[0]
        if sel.size == 0:
            continue  # non-leptonic decay; skip (realizes the leptonic BR)
        px, py, pz, E = np.asarray(prod.P)[sel[0]]
        pz_star[i] = pz
        pT_star[i] = np.hypot(px, py)
        E_star[i] = E
        i += 1
    return pz_star, pT_star, E_star


def _get_pool(nusign, helicity, pool_size, rng):
    key = (1 if nusign > 0 else -1, round(float(helicity), 6), pool_size)
    pool = _POOLS.get(key)
    if pool is None:
        pool = _build_pool(nusign, helicity, pool_size, rng)
        _POOLS[key] = pool
    return pool


def decay_to_muon(E_tau, rng, nusign=1, polarization=None, pool_size=POOL_SIZE):
    """
    Decay an array of taus (each travelling along +z with energy E_tau [GeV]) to a
    muon, returning the lab-frame muon (E_mu [GeV], theta [rad], phi [rad]).

    theta is the muon angle relative to the tau (= beam) direction; phi is the
    azimuth about it.

    nusign       : +1 for neutrino running (tau^-), -1 for antineutrino (tau^+).
    polarization : tau helicity along the beam. None -> V-A default (-1 for tau^-,
                   +1 for tau^+). Pass 0.0 for the unpolarized cross-check.
    """
    E_tau = np.asarray(E_tau, dtype=float)
    n = E_tau.size

    if polarization is None:
        helicity = -1.0 if nusign > 0 else 1.0
    else:
        helicity = float(polarization)

    pz_star, pT_star, E_star = _get_pool(nusign, helicity, pool_size, rng)

    idx = rng.integers(0, pz_star.size, n)
    pzs = pz_star[idx]
    pTs = pT_star[idx]
    Es = E_star[idx]
    phi = rng.uniform(0.0, 2.0 * np.pi, n)  # azimuthal symmetry about the beam axis

    # boost along +z into the lab (pT is unchanged by a z-boost)
    gamma = E_tau / M_TAU
    beta = np.sqrt(1.0 - 1.0 / gamma**2)
    E_mu = gamma * (Es + beta * pzs)
    pz = gamma * (pzs + beta * Es)

    theta = np.arctan2(pTs, pz)
    return E_mu, theta, phi


if __name__ == "__main__":
    rng = np.random.default_rng(0)

    # Unpolarized alouette vs analytic Michel: rest-frame <x> should match ~0.70.
    pz0, pT0, E0 = _get_pool(1, 0.0, 20_000, rng)
    x_alouette = 2.0 * E0 / M_TAU
    x_michel = sample_michel_x(200_000, rng)
    print("rest-frame <x> = 2 E*_mu / m_tau:")
    print("  alouette unpolarized : %.4f" % x_alouette.mean())
    print("  analytic Michel      : %.4f  (massless reference)" % x_michel.mean())

    # Polarization asymmetry: <cos theta*> of the muon about the beam axis.
    for hel in (0.0, -1.0, +1.0):
        pz, pT, E = _get_pool(1, hel, 20_000, rng)
        cstar = pz / np.sqrt(pz**2 + pT**2)
        print("  helicity %+.0f : <cos theta*> = %+.4f" % (hel, cstar.mean()))

    # End-to-end: boosted muon energy fraction should track the rest-frame <x>.
    E_tau = np.full(200_000, 25.0)
    E_mu, theta, phi = decay_to_muon(E_tau, rng)
    print("lab (E_tau=25 GeV, tau^- V-A): <E_mu>/E_tau = %.4f, median theta = %.1f mrad"
          % ((E_mu / E_tau).mean(), np.median(theta) * 1e3))
