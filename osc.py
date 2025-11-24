try:
    import oschelper
except:
    raise ImportError("You must build the oschelper library first!")

import numpy as np

# Neutrino oscillation parameters
th12 = 33.44 *np.pi/180
th23 = 49.20 *np.pi/180
th13 =  8.57 *np.pi/180

dm21 = 7.42e-5
dm31 = 2.517e-3

# matter effect
Ye = 0.5
R_earth = 6371.0 # km

# PREM model
DPdepth = np.genfromtxt("DensityProfileDepthAK135.txt")
DPdensity = np.genfromtxt("DensityProfileDensityAK135.txt")

deltaCP0 = 0.

def oscillate(E, baseline, nuind, nusign, dL=1e2, # m
              rhoset=None, rhoshift=0, rhoscale=1,
              # three-neutrino oscillations
              th12=th12, th23=th23, th13=th13, dm21=dm21, dm31=dm31, deltaCP=deltaCP0,
              # sterile neutrino
              th14=0, th24=0, th34=0, dm41=0,
              # anomalous matter potential
              epsee=0, epsemu=0, epsetau=0, epsmumu=0, epsmutau=0, epstautau=0):

    if isinstance(E, float) or isinstance(E, int):
        E = np.array([E]).astype(float)

    c12 = np.cos(th12)
    c23 = np.cos(th23)
    c13 = np.cos(th13)
    c14 = np.cos(th14)
    c24 = np.cos(th24)
    c34 = np.cos(th34)

    s12 = np.sin(th12)
    s23 = np.sin(th23)
    s13 = np.sin(th13)
    s14 = np.sin(th14)
    s24 = np.sin(th24)
    s34 = np.sin(th34)

    # # Vacuum oscillations -- from https://arxiv.org/pdf/hep-ph/0209097
    # U = np.array([
	# [c12*c13*c14,	c13*c14*s12,	c14*s13,	s14],
	# [-c23*c24*s12 -c12*c24*s13*s23 -c12*c13*s14*s24,	c12*c23*c24 -c24*s12*s13*s23 -c13*s12*s14*s24,	c13*c24*s23 -s13*s14*s24,	c14*s24],
	# [-c12*c23*c34*s13 +c34*s12*s23 -c12*c13*c24*s14*s34 +c23*s12*s24*s34 +c12*s13*s23*s24*s34,	-c23*c34*s12*s13 -c12*c34*s23 -c13*c24*s12*s14*s34 -c12*c23*s24*s34 +s12*s13*s23*s24*s34,	c13*c23*c34 -c24*s13*s14*s34 -c13*s23*s24*s34,	c14*c24*s34],
	# [-c12*c13*c24*c34*s14 +c23*c34*s12*s24 +c12*c34*s13*s23*s24 +c12*c23*s13*s34 -s12*s23*s34,	-c13*c24*c34*s12*s14 -c12*c23*c34*s24 +c34*s12*s13*s23*s24 +c23*s12*s13*s34 +c12*s23*s34,	-c24*c34*s13*s14 -c13*c34*s23*s24 -c13*c23*s34,	c14*c24*c34]
    # ]).astype(np.complex128)

    # U += np.array([
    #     [0, 0, s13*np.exp(-1.0j*deltaCP), 0],
    #     [-c12*s23*s13*np.exp(1.0j*deltaCP), -s12*s23*s13*np.exp(1.0j*deltaCP), 0, 0],
    #     [-c12*c23*s13*np.exp(1.0j*deltaCP), -s12*c23*s13*np.exp(1.0j*deltaCP), 0, 0],
    #     [0, 0, 0, 0]
    # ])

    # 3 neutrino vacuum oscillations -- from https://pdg.lbl.gov/2020/reviews/rpp2020-rev-neutrino-mixing.pdf
    U = np.array([
    [c12*c13,	c13*s12,	s13*np.exp(-1.0j*deltaCP), 0],
    [-c23*s12 -c12*s13*s23*np.exp(1.0j*deltaCP),	c12*c23 -s12*s13*s23*np.exp(1.0j*deltaCP),	c13*s23, 0],
    [-c12*c23*s13*np.exp(1.0j*deltaCP) + s12*s23, -c23*s12*s13*np.exp(1.0j*deltaCP) -c12*s23,	c13*c23, 0],
    [0, 0, 0, 1]
    ]).astype(np.complex128)

    # Matter Potential
    MP = nusign*np.array([
        [1-0.5+epsee, epsemu, epsetau, 0],
        [np.conj(epsemu), -0.5+epsmumu, epsmutau, 0],
        [np.conj(epsetau), np.conj(epsmutau), -0.5+epstautau, 0],
        [0, 0, 0, 0]
    ])

    if nusign < 0:
        U = np.conj(U)

    steps = np.linspace(0, baseline, int(baseline/dL))

    B = steps[1:]/1e3
    cosalpha =  B[-1] / (2*R_earth)

    depths = R_earth - np.sqrt(B**2 + R_earth**2 - 2*B*R_earth*cosalpha)

    density = np.interp(depths, DPdepth, DPdensity, left=DPdensity[0], right=DPdensity[-1])
    density += rhoshift
    density *= rhoscale

    return oschelper.do_osc(E, steps, density, U, MP, dm21, dm31, dm41, nuind, dL)
