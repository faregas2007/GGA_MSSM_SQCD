import numpy as np
from scipy.integrate import odeint

def RGE(api, t, beta0, beta1, beta2, beta3):
	"""
	RG-equation: (d api)/(d log(mu^2)) = api*beta(api), t = logmu2
	"""
	return api*(-beta0*api - beta1*api**2 - beta2*api**3 - beta3*api**4)


def inibeta(nf, nloopin):
	"""
	initial b_i beta function
	SM model ?
	"""
	z3 = 1.20205690031595942853997

	beta0 = (33 - 2*nf)/12.0
	beta1 = (102 - (38*nf)/3.0)/16.0
	beta2 = (2857/2.0 - (5033*nf)/18.0 + (325*nf**2)/54.0)/64.0
	beta3 = (149753/6.0 + (1093*nf**3)/729.0 + 3564*z3 + nf**2*(50065/162.0 + (6472*z3)/81.0)
		- nf*(1078361/162.0 + (6508*z3)/27.0))/256.0
	nloop = nloopin

	if (nloop > 4):
		print("5-loop beta function unknow. Using 4-loop instead")
		nloop = 4
	if (nloop < 4):
		beta3 = 0.0
		if (nloop < 3):
			beta2 = 0.0
			if (nloop < 2):
				beta1 = 0.0
				if (nloop < 1):
					beta0 = 0.0

	return beta0, beta1, beta2, beta3

def runalphas(api0, mu0, mu, nf, nloop, size):
	"""
	purpose: computes the value of api(mu) from api(mu0)
	method: RG-equation from Adaptive Runge-kutta method

	api = alphas/pi
	api0 : api(mu0)
	nf : number of flavor
	nloop: number of loops
	size: running points
	apiout: api(mu)
	"""

	if (nloop == 0):
		return api0


	# integration bounds (logmu2 is the integration variable)
	t0 = 0.0
	tf = 2.0*np.log(mu/mu0)
	#apif = api0
	t = np.linspace(t0, tf, size)

	# intitalize beta function
	beta0, beta1, beta2, beta3 = inibeta(nf, nloop)

	# chek if input are reasonable
	dlam = mu0*np.exp(-1.0/(2.0*beta0*api0))
	if (mu  < dlam):
		print(dlam+" "+mu+" "+mu0+" "+api0*np.pi)

	# integrate RG-equation
	sol = odeint(RGE, api0, t, args=(beta0, beta1, beta2, beta3))
	return sol[-1]

"""
api0 = 1.1720000e-1/np.pi
mu0 = np.float(9.1876000e+01)
mu = 1000.0
nf = 5
nloop =3
size = 100
runalphas(api0, mu0, mu, nf, nloop, size)
"""

def runmass(mass0, api0, apif, nf, nloop):
	"""
	evaluate the running of the MS-bar quark mass
	m(mu) = m(mu0)*exp(\int_a0^af dx gamma(x)/x/beta(x))
	in terms of alphas.

	Inpput:
	mass0 : m(mu0)
	api0 : alphas(mu0)/pi
	apif : alphas(mu)/pi
	nf : number of flavor
	nloop : order of calculation (nloop=1..4)

	Output:
	massout: m(muf)
	"""

	if (nloop == 0):
		return mass0

	gamma0, gamma1, gamma2, gamma3 = inigamma(nf, nloop)
	beta0, beta1, beta2, beta3 = inibeta(nf, nloop)

	bb1 = beta1/beta0
	bb2 = beta2/beta0
	bb3 = beta3/beta0

	cc0 = gamma0/beta0
	cc1 = gamma1/beta0
	cc2 = gamma2/beta0
	cc3 = gamma3/beta0

	cfunc1 = 1.0
	cfunc2 = cc1 - bb1*cc0
	cfunc3 = 1/2.0*((cc1 - bb1*cc0)**2 + cc2 - bb1*cc1 + bb1**2*cc0 - bb2*cc0)
	cfunc4 = (1/6.0*(cc1 - bb1*cc0)**3 + 1/2.0*(cc1 - bb1*cc0)*(cc2 -bb1*cc1 + bb1**2*cc0-bb2*cc0)
		+1/3.0*(cc3 - bb1*cc2 + bb1**2*cc1 - bb2*cc1 - bb1**3*cc0 + 2*bb1*bb2*cc0 - bb3*cc0))

	if (nloop < 4):
		cfunc4 = 0.0
		if (nloop < 3):
			cfunc3 = 0.0
			if (nloop < 2):
				cfunc2 = 0.0
				if (nloop < 1):
					cfunc1 = 0.0

	cfuncmu0 = cfunc1 + cfunc2*api0 + cfunc3*api0**2 + cfunc4*api0**3
	cfuncmuf = cfunc1 + cfunc2*apif + cfunc3*apif**2 + cfunc4*apif**3

	massout = mass0*(apif/api0)**cc0*cfuncmuf/cfuncmu0
	return massout

def inigamma(nfin, nloopin):
	"""
	initial gammma functions
	"""
	z3 = 1.2020569031595942853997
	z5 = 1.0369277551433699263

	nf = nfin

	gamma0 = 1.0
	gamma1 = (67.33333333333333 - (20*nf)/9.0)/16.0
	gamma2 = (1249.0 - (140*nf**2)/81.0 + 2*nf*(-20.59259259259259 - 48*z3) +(8*nf*(-46 + 48*z3))/9.0)/64.0
	gamma3 = (28413.919753086420 + (135680*z3)/27.0 + nf**3*(-1.3662551440329218 + (64*z3)/27.0)
		+ nf**2*(21.57201646090535 - (16*np.pi**4)/27.0 + (800*z3)/9.0) + 8800*z5
		+ nf*(-3397.1481481481483 + (88*np.pi**4)/9.0 - (34192*z3)/9.0 + (18400*z5)/9.0))/256.0
	nloop = nloopin

	if (nloop > 4):
		print ("5-loop gamma function unknown. Using 4-loop instead")
		nloop =4

	if (nloop < 4):
		gamma3 = 0.0
		if (nloop < 3):
			gamma2 = 0.0
			if (nloop < 2):
				gamma1 = 0.0
				if (nloop < 1):
					gamma0 = 0.0
	return gamma0, gamma1, gamma2, gamma3
