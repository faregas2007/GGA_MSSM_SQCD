import numpy as np
import const
from getdata import *

# C0 function
def Ifunc(a, b, c):
	return (a*b*np.log(a/b) + b*c*np.log(b/c) + a*c*np.log(c/a))/((a-b)*(b-c)*(a-c))

def ftau(x):
	if x>= 1:
		ftau = (np.arcsin(1.0/np.sqrt(x)))**2
	else:
		ftau = -(np.log((1 + np.sqrt(1-x))/(1-np.sqrt(1-x))) - complex(0.0,1.0)*np.pi)**2/4.0
	return ftau

# quark contribution
def Aqtau(mq, mhiggs, lpesudo):
	tau = (2*mq/mhiggs)**2
	if not lpesudo:
		Aqtau = -(3.0/2.0)*tau*(1 + (1-tau)*ftau(tau))
	else:
		Aqtau = tau*ftau(tau)
	return Aqtau

# squark contributions
def Asqtau(msq, mhiggs, lpesudo):
	tau = (2*msq/mhiggs)**2
	if not lpesudo:
		Asqtau = -3.0/4.0 * tau * (1+ tau*ftau(tau))
	else:
		Asqtau = 0.0
	return Asqtau

# squared amplitude
def AMPQ(mt, mb, mG, lpesudo, tanbresum, tantresum, alphas, norder):
	"""
	CqQCD: NLO QCD contributions from quark part contributions --> todo list.
	CqSQCD: NLO SQCD contributions from quark part contributions
	mG: will use Gluino mass, fixed parameters in Micheal input
	mu: will use muSUSY scale, fixed parameters in Micheal input
	"""

	CqQCD = 0.0

	AMPQ = []
	AMPQ0 = []
	AMPQe = []

	CtQCD = 0.0
	CbQCD = 0.0

	# get Cdata from Cdatatable
	params = getdata('summaryfile.csv')

	reCt = params.reCt
	reCterr = params.reCterr
	imCt = params.imCt
	imCterr = params.imCterr
	reCb = params.reCb
	reCberr = params.reCberr
	imCb = params.imCb
	imCberr = params.imCberr
	mu = params.mu
	At = params.At
	Ab = params.Ab

	mA = params.mA
	tanb = params.tanb

	mst1 = params.mst1
	mst2 = params.mst2
	msb1 = params.msb1
	msb2 = params.msb2
	alphaqst = params.alphaqst
	alphaqsb = params.alphaqsb

	gb = tanb
	gt = 1.0/tanb

	CF = 4.0/3.0
	for idx in range(len(mA)):
		temp = 0.0
		tempe = 0.0
		CtSQCD = complex(0.0, 0.0)
		CbSQCD = complex(0.0, 0.0)
		deltab = (CF/2.0)*(alphaqsb[idx]/np.pi)*mG*mu[idx]*tanb[idx]*Ifunc(msb1[idx]**2, msb2[idx]**2, mG**2)
		deltat = (CF/2.0)*(alphaqst[idx]/np.pi)*mG*mu[idx]*(1/tanb[idx])*Ifunc(mst1[idx]**2, mst2[idx]**2, mG**2)

		if (norder == 0):
			CtSQCD = complex(0.0,0.0)
			CtSQCDe = complex(0.0,0.0)
			CbSQCD = complex(0.0,0.0)
			CbSQCDe = complex(0.0,0.0)
			CbLE = 0.0
			CtLE = 0.0
		elif (norder == 1):
			CtSQCD = complex(reCt[idx], imCt[idx])
			CtSQCDe = complex(reCterr[idx], imCterr[idx])
			CbSQSCD = complex(reCb[idx], imCb[idx])
			CbSQCDe = complex(reCberr[idx], imCberr[idx])

			CbLE = deltab*(1.0 + 1.0/tanb[idx]**2)*np.pi/alphaqsb[idx]
			CtLE = deltat*(1.0 + tanb[idx]**2)*np.pi/alphaqst[idx]

		A_bv = Aqtau(mb, mA[idx], lpesudo)
		A_tv = Aqtau(mt, mA[idx], lpesudo)

		gbeff = gb[idx]*(1.0 - deltab/tanb[idx]**2)/(1 + deltab)
		gteff = gt[idx]*(1.0 - deltat*tanb[idx]**2)/(1 + deltat)

		if (tantresum == 1 and tanbresum == 1):
			Atq = A_tv*gteff
			Abq = A_bv*gbeff
		elif (tantresum == 0 and tanbresum == 0):
			Atq = A_tv*gt[idx]
			Abq = A_bv*gb[idx]
			CtLE = 0.0
			CbLE = 0.0
		elif (tantresum == 1 and tanbresum == 0):
			Atq = A_tv*gteff
			Abq = A_bv*gb[idx]
			CbLE = 0.0
		elif (tantresum == 0 and tanbresum == 1):
			Atq = A_tv*gt[idx]
			Abq = A_bv*gbeff
			CtLE = 0.0

		temp = np.abs(Atq + Abq)**2
		AMPQ0.append(temp)
		temp = temp + (Atq + Abq).conjugate()*(gt[idx]*A_tv*(CtSQCD-CtLE) + gb[idx]*A_bv*(CbSQCD - CbLE))*(alphas[idx]/np.pi)
		AMPQ.append(temp)

		AMPQeb = 2.0*(Atq + Abq).conjugate()*(gb[idx]*A_bv*CbSQCDe)*(alphas[idx]/np.pi)
		AMPQet = 2.0*(Atq + Abq).conjugate()*(gt[idx]*A_tv*CtSQCDe)*(alphas[idx]/np.pi)
		tempe = np.sqrt(np.abs(AMPQeb)**2 + np.abs(AMPQet)**2)
		AMPQe.append(tempe)
	
	return AMPQ0, AMPQ, AMPQe
