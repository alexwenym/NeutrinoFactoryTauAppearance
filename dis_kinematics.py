import os

import numpy as np

# Deep-inelastic-scattering (DIS) production kinematics for neutrino CC interactions.
#
# Given a neutrino energy E_nu, this samples the outgoing charged lepton (E_lep, theta,
# phi) -- theta measured from the neutrino (= beam) direction -- from the double-
# differential cross section d^2 sigma / dx dy (x = Bjorken-x, y = inelasticity). The
# muon scattering angle (the make-or-break observable separating the direct-nu_mu core
# from the tau-decay halo, docs/00-physics.md) comes from the momentum transfer
#
#     E_lep = (1 - y) E_nu
#     Q^2   = 2 M_N x y E_nu
#     cos(theta) = (2 E_nu E_lep - Q^2 - m_lep^2) / (2 E_nu p_lep)      (exact, massive)
#
# CROSS SECTIONS: sigma(E, channel) is the integral of the same d^2sigma/dx dy, so
# kinematics and total rate are self-consistent (returned in 1e-38 cm^2, matching totalXS).
#
# DATA: GENIE double-differential DIS tables (tune G18_02a, E=5-200 GeV) produced with a
# custom gDISDiffXSec app, symlinked at resources/diffxsec (-> software/diffxsec/tables).
# 24 files = 6 flavors x {p,n} x {CC,NC}; we use the CC isoscalar average (p+n)/2.
# TableDISModel is the default; AnalyticDISModel (LO valence+sea) remains as a fallback.

M_N = 0.9383      # GeV, isoscalar nucleon
M_MU = 0.105658   # GeV
M_TAU = 1.77686   # GeV
M_E = 0.000511    # GeV

CHANNELS = ("numu", "numubar", "nutau", "nutaubar", "nue", "nuebar")
_LEPTON_MASS = {"numu": M_MU, "numubar": M_MU, "nutau": M_TAU,
                "nutaubar": M_TAU, "nue": M_E, "nuebar": M_E}
_IS_NUBAR = {"numu": False, "numubar": True, "nutau": False,
             "nutaubar": True, "nue": False, "nuebar": True}

DEFAULT_TABLE_DIR = "resources/diffxsec"


def _edges_from_centers(c):
    """Cell edges bracketing grid-point centers c (for cell-width weighting/sampling)."""
    c = np.asarray(c, dtype=float)
    mid = 0.5 * (c[1:] + c[:-1])
    first = c[0] - (mid[0] - c[0])
    last = c[-1] + (c[-1] - mid[-1])
    return np.concatenate(([first], mid, [last]))


# ---- cross-section models: provide grid_and_values(channel) ---------------------------

class TableDISModel:
    """Loads GENIE CC double-differential DIS tables (isoscalar (p+n)/2) on their native
    (E, x, y) grid. Tables are ASCII: columns E[GeV], x, y, d2sigma_dxdy[cm^2], ordered
    y-inner / x-middle / E-outer."""

    def __init__(self, table_dir=DEFAULT_TABLE_DIR, current="ccdis"):
        self.table_dir = table_dir
        self.current = current
        self._grid = None
        self._values = {}

    def _path(self, channel, target):
        return os.path.join(self.table_dir,
                            "%s_%s_%s.dat" % (channel, target, self.current))

    def _load_grid(self):
        if self._grid is not None:
            return
        ref_path = self._path("numu", "p")
        if not os.path.exists(ref_path):
            raise FileNotFoundError(
                "GENIE DIS tables not found at %r. Expected the symlink "
                "resources/diffxsec -> software/diffxsec/tables (see dis_kinematics docs)."
                % self.table_dir)
        ref = np.loadtxt(ref_path)
        nE = np.unique(ref[:, 0]).size
        nx = np.unique(ref[:, 1]).size
        ny = np.unique(ref[:, 2]).size
        if nE * nx * ny != ref.shape[0]:
            raise ValueError("table %s is not a full regular grid" % ref_path)
        self.nE, self.nx, self.ny = nE, nx, ny
        self.e_nodes = ref[:, 0].reshape(nE, nx, ny)[:, 0, 0].copy()
        self.xc = ref[:, 1].reshape(nE, nx, ny)[0, :, 0].copy()
        self.yc = ref[:, 2].reshape(nE, nx, ny)[0, 0, :].copy()
        self._grid = (self.e_nodes, self.xc, self.yc)

    def grid_and_values(self, channel):
        if channel not in CHANNELS:
            raise ValueError("unknown channel %r" % channel)
        self._load_grid()
        if channel not in self._values:
            fp = np.loadtxt(self._path(channel, "p"), usecols=3)
            fn = np.loadtxt(self._path(channel, "n"), usecols=3)
            iso = 0.5 * (fp + fn)  # isoscalar target
            self._values[channel] = iso.reshape(self.nE, self.nx, self.ny)
        return self.e_nodes, self.xc, self.yc, self._values[channel]


class AnalyticDISModel:
    """Leading-order DIS stand-in (valence+sea, d^2sigma/dx dy ~ x[q + qbar(1-y)^2]).
    Fallback only; the GENIE tables (TableDISModel) are the default."""

    def __init__(self, sea_frac=0.15, nx=120, ny=120, e_nodes=None):
        self.sea_frac = sea_frac
        self.nx, self.ny = nx, ny
        self.e_nodes = np.linspace(5.0, 60.0, 40) if e_nodes is None else np.asarray(e_nodes)
        self.xc = np.linspace(1e-3, 1.0, nx, endpoint=False) + 0.5 / nx
        self.yc = (np.arange(ny) + 0.5) / ny

    def grid_and_values(self, channel):
        m = _LEPTON_MASS[channel]
        nub = _IS_NUBAR[channel]
        X = self.xc[:, None]
        Y = self.yc[None, :]
        xq = np.sqrt(X) * (1.0 - X) ** 3 + self.sea_frac * (1.0 - X) ** 7
        xqbar = self.sea_frac * (1.0 - X) ** 7
        shape = (xqbar + xq * (1.0 - Y) ** 2) if nub else (xq + xqbar * (1.0 - Y) ** 2)
        vals = np.empty((self.e_nodes.size, self.nx, self.ny))
        for k, E in enumerate(self.e_nodes):
            vals[k] = np.where((1.0 - Y) * E > m, shape, 0.0)
        return self.e_nodes, self.xc, self.yc, vals


# numu CC slope (PDG, 1e-38 cm^2/GeV) -- only used to anchor the analytic fallback's scale
_SIGMA_SLOPE_NUMU = 0.677


class DISSampler:
    """Builds per-energy-node 2D (x, y) CDFs from a gridded DIS model for one channel and
    samples (E_lep, theta, phi) for arrays of E_nu; also provides sigma(E)."""

    def __init__(self, model, channel, sigma_scale=None):
        self.channel = channel
        self.mass = _LEPTON_MASS[channel]
        e_nodes, xc, yc, vals = model.grid_and_values(channel)
        self.e_nodes = np.asarray(e_nodes, dtype=float)
        self.nE = self.e_nodes.size
        self.xedges = _edges_from_centers(xc)
        self.yedges = _edges_from_centers(yc)
        self.xedges[0] = max(self.xedges[0], 1e-6)
        xw = np.diff(self.xedges)
        yw = np.diff(self.yedges)
        self.nx, self.ny = xc.size, yc.size

        # cell weight = d2sigma/dxdy * dx * dy ; per-node integral and flattened CDF
        W = np.clip(vals, 0.0, None) * xw[None, :, None] * yw[None, None, :]
        self.integral = W.reshape(self.nE, -1).sum(axis=1)          # cm^2 per node
        flat = W.reshape(self.nE, self.nx * self.ny)
        c = np.cumsum(flat, axis=1)
        tot = c[:, -1:].copy()
        tot[tot == 0] = 1.0
        self.cdf = c / tot

        # geometric midpoints for nearest-E-node digitizing (E grid is log-spaced)
        self.e_edges = np.concatenate(
            ([-np.inf], np.sqrt(self.e_nodes[1:] * self.e_nodes[:-1]), [np.inf]))

        # sigma units: tables are cm^2 -> report in 1e-38 cm^2 (matches totalXS).
        # Analytic fallback carries no absolute scale, so anchor it to the PDG numu slope.
        self.sigma_scale = sigma_scale

    def sigma(self, E):
        """Total CC cross section [1e-38 cm^2]."""
        E = np.asarray(E, dtype=float)
        I = np.interp(E, self.e_nodes, self.integral,
                      left=self.integral[0], right=self.integral[-1])
        if self.sigma_scale is not None:        # analytic fallback: scale * E * I
            return self.sigma_scale * E * I
        return I / 1e-38                         # tables: integral already in cm^2

    def sample(self, E_nu, rng):
        """Return lab-frame lepton (E_lep [GeV], theta [rad], phi [rad]) per E_nu,
        theta measured from the neutrino (= beam) direction."""
        E_nu = np.asarray(E_nu, dtype=float)
        n = E_nu.size
        node = np.clip(np.searchsorted(self.e_edges, E_nu) - 1, 0, self.nE - 1)

        x = np.empty(n)
        y = np.empty(n)
        u = rng.uniform(size=n)
        jx = rng.uniform(size=n)
        jy = rng.uniform(size=n)
        ncell = self.nx * self.ny
        for k in range(self.nE):
            m = node == k
            if not np.any(m):
                continue
            flat = np.clip(np.searchsorted(self.cdf[k], u[m]), 0, ncell - 1)
            ix = flat // self.ny
            iy = flat % self.ny
            x[m] = self.xedges[ix] + jx[m] * (self.xedges[ix + 1] - self.xedges[ix])
            y[m] = self.yedges[iy] + jy[m] * (self.yedges[iy + 1] - self.yedges[iy])

        E_lep = np.maximum((1.0 - y) * E_nu, self.mass)
        Q2 = 2.0 * M_N * x * y * E_nu
        p_lep = np.sqrt(np.maximum(E_lep**2 - self.mass**2, 0.0))
        with np.errstate(divide="ignore", invalid="ignore"):
            cos_t = (2.0 * E_nu * E_lep - Q2 - self.mass**2) / (2.0 * E_nu * p_lep)
        cos_t = np.where(p_lep > 0.0, cos_t, 1.0)
        theta = np.arccos(np.clip(cos_t, -1.0, 1.0))
        phi = rng.uniform(0.0, 2.0 * np.pi, n)
        return E_lep, theta, phi


# ---- module-level convenience (cached samplers per channel) ---------------------------

_MODEL = None
_SAMPLERS = {}


def _default_model():
    global _MODEL
    if _MODEL is None:
        _MODEL = TableDISModel()
    return _MODEL


def set_model(model):
    """Swap the DIS model (e.g. AnalyticDISModel() fallback) and clear cached samplers."""
    global _MODEL
    _MODEL = model
    _SAMPLERS.clear()


def get_sampler(channel):
    s = _SAMPLERS.get(channel)
    if s is None:
        model = _default_model()
        scale = None
        if isinstance(model, AnalyticDISModel):
            # anchor analytic fallback so sigma_numu = PDG slope * E
            tmp = DISSampler(model, "numu")
            I_ref = float(np.interp(1e3, tmp.e_nodes, tmp.integral,
                                    left=tmp.integral[0], right=tmp.integral[-1]))
            scale = _SIGMA_SLOPE_NUMU / (1e3 * I_ref)
        s = DISSampler(model, channel, sigma_scale=scale)
        _SAMPLERS[channel] = s
    return s


def sample_dis_lepton(E_nu, rng, channel):
    """Sample the outgoing CC lepton (E_lep, theta, phi) for `channel` at energies E_nu."""
    return get_sampler(channel).sample(E_nu, rng)


def sigma(E, channel):
    """Total CC cross section [1e-38 cm^2] for `channel`."""
    return get_sampler(channel).sigma(E)


if __name__ == "__main__":
    rng = np.random.default_rng(0)

    print("== model:", type(_default_model()).__name__,
          "  table_dir:", getattr(_default_model(), "table_dir", "-"))

    print("== mean inelasticity <y> and median theta (E_nu = 20 GeV) ==")
    E = np.full(400_000, 20.0)
    for ch in ("numu", "numubar", "nutau", "nutaubar"):
        El, th, ph = sample_dis_lepton(E, rng, ch)
        print("  %-9s <y>=%.3f  median theta=%5.1f mrad  <E_lep>/E=%.3f"
              % (ch, (1.0 - El / E).mean(), np.median(th) * 1e3, (El / E).mean()))

    print("== sigma(E)/E  [1e-38 cm^2 / GeV] ==")
    for ch in ("numu", "numubar", "nutau"):
        row = "  %-7s" % ch
        for Etest in (10.0, 25.0, 50.0):
            row += "  E=%4.0f: %.3f" % (Etest, float(sigma(Etest, ch)) / Etest)
        print(row)

    # cross-check tau DIS sigma against the digitized TOTAL nu_tau CC (they should agree
    # at high E where DIS dominates and differ near threshold where QE/RES matter).
    try:
        from cross_sections import totalXS
        xs = totalXS("resources/nutau_xs.txt", "resources/nutaubar_xs.txt")
        print("== nu_tau: GENIE DIS sigma vs digitized TOTAL [1e-38 cm^2] ==")
        for Etest in (10.0, 25.0, 50.0):
            print("  E=%4.0f  DIS=%.3f  total(digitized)=%.3f"
                  % (Etest, float(sigma(Etest, "nutau")), float(xs.sigma(Etest))))
    except Exception as e:
        print("  (skipped digitized cross-check:", e, ")")

    print("== theta -> 0 as y -> 0 (numu) ==")
    El, th, ph = sample_dis_lepton(np.full(500_000, 20.0), rng, "numu")
    yv = 1.0 - El / 20.0
    print("  median theta  y<0.2: %.1f mrad   y>0.6: %.1f mrad"
          % (np.median(th[yv < 0.2]) * 1e3, np.median(th[yv > 0.6]) * 1e3))
