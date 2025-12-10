import numpy as np

# XS are digitized from https://arxiv.org/pdf/2409.01258

class totalXS:

    def __init__(self, nu_file, nubar_file, delimiter=','):

        self.delimiter = delimiter

        # Load nu XS
        self._x_nu, self._y_nu = self._read_xs_file(nu_file)
        # Load nubar XS: use same file if not given
        self._x_nubar, self._y_nubar = self._read_xs_file(nubar_file)
        
        # store min/max for fast checking
        self._min_x_nu = self._x_nu[0]
        self._max_x_nu = self._x_nu[-1]

        self._min_x_nubar = self._x_nubar[0]
        self._max_x_nubar = self._x_nubar[-1]

    def _read_xs_file(self, filename):
        """
        Read a cross-section file with a possible header E [GeV], sigma [1e-38 cm2]
        and return sorted arrays.
        """
        xs = []
        ys = []
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Skip header or commented lines
                if line.lower().startswith('x') or line.startswith('#'):
                    continue
                parts = [p.strip() for p in line.split(self.delimiter)]
                if len(parts) < 2:
                    continue
                x_val = float(parts[0])
                y_val = float(parts[1])
                xs.append(x_val)
                ys.append(y_val)

        x_arr = np.array(xs)
        y_arr = np.array(ys)

        # Ensure ascending order in x
        order = np.argsort(x_arr)
        return x_arr[order], y_arr[order]
    
    def _check_bounds(self, x, min_x, max_x, label=""):
        """
        Warn if x is outside [min_x, max_x].
        """
        x = np.asarray(x)
        if np.any(x < min_x) or np.any(x > max_x):
            warnings.warn(
                f"[totalXS] Energy outside interpolation range for {label}: "
                f"valid range = [{min_x}, {max_x}] GeV. "
                f"Values outside will be clamped to the nearest boundary.",
                RuntimeWarning
            )
            
    def get_range(self):
        return [[self._min_x_nu, self._max_x_nu],[self._min_x_nubar, self._max_x_nubar]]
    
    def sigma_nu(self, x):
        """
        inputs E[GeV] and return sigma [1e-38 cm2]
        """
        self._check_bounds(x, self._min_x_nu, self._max_x_nu, label="nutau")
        
        x = np.asarray(x)
        
        return np.interp(x, self._x_nu, self._y_nu)

    def sigma_nubar(self, x):
        """
        inputs E[GeV] and return sigma [1e-38 cm2]
        """
        self._check_bounds(x, self._min_x_nubar, self._max_x_nubar, label="nutaubar")
        
        x = np.asarray(x)
        
        return np.interp(x, self._x_nubar, self._y_nubar)

    def sigma(self, x, is_nubar=False):

        if is_nubar:
            return self.sigma_nubar(x)
        else:
            return self.sigma_nu(x)