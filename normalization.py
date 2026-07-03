import numpy as np

import flux

# Absolute normalization for neutrino-induced through-going muons.
#
# The taus/muons are produced by nu CC interactions in the ROCK upstream of the detector
# (the detector is the tracker, not the target), all along the last ~muon-range of rock.
# The beam is ~uniform transversely over the muon lateral spread (a far detector sees a
# km-scale beam spot), so by translational invariance the muon flux crossing the plane is
# a local per-area quantity, independent of detector size:
#
#   Phi_mu [muons/cm^2/yr] = int dE  Phi_nu(E) P(E) sigma(E) * n_rock * R_eff(E) * BR
#
# where R_eff is the effective longitudinal reach. The MC samples E ~ Phi_nu*P*sigma and
# production point d ~ U(0, dmax), so int Phi_nu P sigma * R_eff = dmax * w_integral *
# accept_frac (accept_frac is the rate-weighted mean acceptance over d and angle). The
# MC's w_integral uses Nmuon=1 (flux per m^2) and sigma in 1e-38 cm^2, so:
#
#   Phi_mu = (Nmuon/1e4) * 1e-38 * w_integral * n_rock * (dmax*100) * accept_frac
#            * det_frac * BR
#
# Detector transverse size does NOT enter Phi_mu (it's per cm^2); the total detected
# count is Phi_mu * area * years -- so "infinitely large / cut later" just means we report
# the per-area flux and scale by whatever area is chosen.

N_A = 6.02214e23
RHO_ROCK = 2.65                 # g/cm^3, standard rock
N_ROCK = RHO_ROCK * N_A         # nucleons/cm^3 (~1 g/mol per nucleon) ~ 1.60e24


def muon_flux(out, Emuon, variant="accelerator"):
    """Absolute muon flux crossing the detector plane [muons / cm^2 / yr] for one
    channel's simulation output dict (see module docstring for the derivation)."""
    Nmuon = flux.n_muon_decays(Emuon, variant)
    dmax_cm = out["dmax"] * 100.0
    det_frac = out.get("det_frac", 1.0)
    return ((Nmuon / 1e4) * 1e-38 * out["w_integral"] * N_ROCK * dmax_cm
            * out["accept_frac"] * det_frac * out["BR"])


def counts(out, Emuon, variant="accelerator", area_m2=100.0, years=10.0):
    """Expected detected muons for a detector transverse area [m^2] and exposure [yr]."""
    return muon_flux(out, Emuon, variant) * (area_m2 * 1e4) * years


if __name__ == "__main__":
    # sanity: beam intensity and a core-channel muon flux magnitude
    for E in (25.0, 50.0):
        for v in ("accelerator", "dump", "baseline"):
            print("Nmuon(E=%2.0f, %-11s) = %.3e /yr" % (E, v, flux.n_muon_decays(E, v)))
    from simulate import simulate_numu_cc_muons, simulate_tau_muons
    core = simulate_numu_cc_muons(Emuon=25.0, baseline=1300e3, n_mc=400_000, mcs=True)
    halo = simulate_tau_muons(Emuon=25.0, baseline=1300e3, n_mc=400_000, mcs=True)
    for name, o in (("core", core), ("halo", halo)):
        phi = muon_flux(o, 25.0, "accelerator")
        n = counts(o, 25.0, "accelerator", area_m2=100.0, years=10.0)
        print("%s: Phi_mu = %.3e muons/cm^2/yr   N(100 m^2, 10 yr) = %.3e" % (name, phi, n))
