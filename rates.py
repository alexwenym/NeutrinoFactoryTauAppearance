from osc import *
from flux import *

class TauAppearanceRates:


    def __init__(self,
                 baseline, off_axis_angle, detector_radius, # detector parameters
                 Emuon, Pmuon, mucharge=-1, Nmuon = 1e22, # muon beam parameters
                 deltaCPvals=[0,np.pi/2,np.pi,3*np.pi/2], # oscillation parameters
                 N_e_points = 100,
                 ):

        self.baseline = baseline
        self.off_axis_angle = off_axis_angle
        self.detector_solid_angle = detector_radius**2 / (self.baseline**2) * np.pi

        self.Emuon = Emuon
        self.Nmuon = Nmuon
        self.Pmuon = Pmuon
        self.mucharge = mucharge
        if self.mucharge == -1:
            self.nusign = 1 # nue from mu- decay
        elif self.mucharge == 1:
            self.nusign = -1 # nuebar from mu+ decay
        else:
            print("invalid muon charge %s"%str(self.mucharge))


        self.enurange = np.linspace(0, Emuon, N_e_points)

        self.deltaCPvals = deltaCPvals

        self._compute_flux()
        self._compute_osc()

    def _compute_flux(self):
        self.flux = nue_flux(self.Emuon, self.Nmuon, self.Pmuon, self.enurange, np.cos(self.off_axis_angle), self.baseline, mucharge=self.mucharge) # per GeV per sr

    def _compute_osc(self):
        self.osc_prob = {}
        for deltaCP in self.deltaCPvals:
           self.osc_prob[deltaCP] =  np.abs(oscillate(self.enurange,self.baseline,0,self.nusign,deltaCP=deltaCP))**2

    def nus_thru_detctor(self,osc=True):
        if osc:
            return {deltaCP:self.flux * self.detector_solid_angle * self.osc_prob[deltaCP] for deltaCP in self.deltaCPvals}
        else:
            return self.flux * self.detector_solid_angle # per GeV