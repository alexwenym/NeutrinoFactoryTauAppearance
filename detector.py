import numpy as np

# Idealized detection-plane ("LArTPC") model -- a downstream, truth-preserving layer.
#
# The simulations (simulate_tau_muons / simulate_numu_cc_muons) record, for every muon
# that reaches the detector plane, its truth-level arrival angle theta, transverse
# crossing position r (and azimuth phi), energy E_mu. THAT dict already *is* the
# detection-plane readout. This module applies an optional, idealized detector response
# on top of it -- a finite footprint, Gaussian angular/energy smearing, and an energy
# threshold -- WITHOUT touching the simulation. Everything is off by default, so
# apply_detector(out) returns the truth-level muons unchanged ("record all, cut later").
#
# This is a phenomenology-level model on purpose: simple Gaussian smearing, no detector
# Monte Carlo. The angle is the signal observable, so sigma_theta is the key knob. Here
# sigma_theta is the TPC's intrinsic angular resolution ONLY; multiple Coulomb scattering
# in the rock is a separate propagation effect applied upstream in simulate (mcs=True).


def _smear_angle(theta, phi, sigma_theta, rng):
    """Gaussian-smear the muon direction by sigma_theta applied to each transverse
    projected angle (the usual way TPC resolution / MCS are quoted). Returns (theta, phi).
    """
    n = theta.size
    tx = theta * np.cos(phi) + rng.normal(0.0, sigma_theta, n)
    ty = theta * np.sin(phi) + rng.normal(0.0, sigma_theta, n)
    return np.hypot(tx, ty), np.arctan2(ty, tx)


def apply_detector(out, rng=None, seed=0,
                   radius=None, half_width=None, half_height=None,
                   sigma_theta=0.0, sigma_E_frac=0.0, E_threshold=0.0):
    """
    Apply an idealized detection-plane response to a simulation output dict and return a
    NEW dict (input is not mutated). All effects are off by default -> truth passthrough.

    Footprint (truth crossing position r is used; off if all None):
      radius        -- keep muons with r <= radius [m] (circular footprint).
      half_width,
      half_height   -- keep |x|<=half_width and/or |y|<=half_height [m] (rectangular),
                       where x=r*cos(phi), y=r*sin(phi).
    Smearing (Gaussian; off if 0):
      sigma_theta   -- TPC per-projection angular resolution [rad] (intrinsic only; rock
                       MCS is applied upstream in simulate). Applied to the measured theta.
      sigma_E_frac  -- fractional muon energy resolution (E -> E*(1+N(0,sigma_E_frac))).
    Threshold (off if 0):
      E_threshold   -- drop muons with measured E_mu < E_threshold [GeV].

    The returned dict keeps the same keys; the scalars (dmax, accept_frac, w_integral,
    BR) are passed through and a 'det_frac' (fraction surviving the footprint/threshold)
    is added.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    theta = np.asarray(out["theta"], dtype=float)
    phi = np.asarray(out["phi"], dtype=float)
    r = np.asarray(out["r"], dtype=float)
    E_mu = np.asarray(out["E_mu"], dtype=float)
    n0 = theta.size

    # measured angle (smeared)
    if sigma_theta > 0.0:
        theta, phi = _smear_angle(theta, phi, sigma_theta, rng)

    # measured energy (smeared), clipped to be positive
    if sigma_E_frac > 0.0:
        E_mu = np.maximum(E_mu * (1.0 + rng.normal(0.0, sigma_E_frac, E_mu.size)), 0.0)

    # footprint on the (truth) crossing position + energy threshold
    keep = np.ones(n0, dtype=bool)
    if radius is not None:
        keep &= r <= radius
    if half_width is not None:
        keep &= np.abs(r * np.cos(phi)) <= half_width
    if half_height is not None:
        keep &= np.abs(r * np.sin(phi)) <= half_height
    if E_threshold > 0.0:
        keep &= E_mu >= E_threshold

    result = dict(out)  # copy scalars (dmax, accept_frac, w_integral, BR, ...)
    for key in ("E_nu", "E_mu", "theta", "r", "phi"):
        if key in out:
            arr = E_mu if key == "E_mu" else (
                theta if key == "theta" else (
                    phi if key == "phi" else np.asarray(out[key])))
            result[key] = arr[keep]
    result["det_frac"] = float(keep.mean()) if n0 else 0.0
    return result


if __name__ == "__main__":
    # smoke test: truth passthrough is a no-op; footprint + smearing reduce/spread.
    from simulate import simulate_tau_muons

    out = simulate_tau_muons(Emuon=25.0, baseline=1300e3, n_mc=200_000)
    truth = apply_detector(out)
    assert truth["theta"].size == out["theta"].size, "default must be a no-op"
    print("truth passthrough: N=%d  median theta=%.1f mrad"
          % (truth["theta"].size, np.median(truth["theta"]) * 1e3))

    det = apply_detector(out, radius=2.0, sigma_theta=0.01, E_threshold=2.0, seed=1)
    print("with R=2 m, sigma_theta=10 mrad, E>2 GeV: det_frac=%.3f  N=%d  median theta=%.1f mrad"
          % (det["det_frac"], det["theta"].size, np.median(det["theta"]) * 1e3))
