import numpy as np
from pylooptools.ltools import *
import pandas as pd
from slha_input import *
from utils import *
from const import *
from integrals import *
from renscheme import *
from reals import *

# QCD corrections part
B1SMA_pd = pd.read_csv('/B1SMA_pd.csv')

def interpolate(df, values, colname1, colname2):
    df = df.iloc[(df[colname1]-values).abs().argsort()[:2]]
    interpolate =  (df[:1][colname2].values - df[1:2][colname2].values)/(df[:1][colname1].values - df[1:2][colname1].values) + df[1:2][colname2].values
    return interpolate[0]

def ImF0AB1A(tau):
    if (tau <= 1.00000001):
        # tau <= 1
        ImF0AB1A = 0.0
    elif (tau >= 10.0):
        # tau >= 10
        ImF0AB1A = 1.45470520 * (tau+0.336397620)*(tau +20.6155650)*(tau**2 -0.283058680 * tau+0.152948090) -0.33173440*(tau-1.5058080)*(tau-0.46593250)*(tau**2 +0.39337320 * tau+0.099801230)* (np.log(tau**2))**2 -2.03308280 *(tau-4.89360940)* (tau-0.127529040)*(tau**2+0.180641150*tau+0.0309529080) * np.log(tau**2) +(0.054541539124820 * tau**4-0.0054541539124820 * tau**2-0.0054541539124820 *tau-0.0051132692929520)*(np.log(tau**2))**3
        ImF0AB1A = ImF0AB1A/(tau**5)
    else:
        # 1 < tau < 10
        # interpolate result
        ImF0AB1A =  interpolate(B1SMA_pd, tau, 'tauimval', 'imvalues')

    return ImF0AB1A

def ReF0AB1A(tau):
    # check if tau was too close to 1, divergence, set reF0AB1A = 0
    # tau >= 20.0
    if (tau >= 20.0):
        ReF0AB1A = -20.0587960*(tau-1.26452900)*(tau+0.376730300)*(tau**2-0.107011770 * tau+0.179938730)+0.0175990560*(tau-1.50580820)*(tau-0.46593250)*(tau**2+0.393373240* tau+0.099801230)*(np.log(tau**2))**3 +0.504482190* (tau-1.68174090)*(tau+0.292852660)*(tau**2-0.163460670 *tau+0.1167291530)*(np.log(tau**2))**2-1.621089300* (tau+0.541737520)*(tau+2.309399280)*(tau**2-0.326851450*tau+0.330847650)*np.log(tau**2)+(-0.00217013888888890*tau**4+0.000217013888888890*tau**2+0.000217013888888890*tau+0.000203450520833330)*(np.log(tau**2))**4
        ReF0AB1A = ReF0AB1A/(tau**5)
    else:
        # tau near threshold
        if(tau > 0.999 and tau < 1.001): 
            ReF0AB1A = 0.0
        # 0.0 <tau < 20.0
        else:
            ReF0AB1A =  interpolate(B1SMA_pd, tau, 'taureval', 'revalues')  
    return ReF0AB1A

def B1SMA(x):
    B1SMA = 3/2.0*complex(ReF0AB1A(x),ImF0AB1A(x))
    return B1SMA

def B2SMA(tau):
    if (tau <= 1.0):
        dfdtau = complex(np.arcsin(np.sqrt(tau)/np.sqrt(tau*(1-tau))))
    else:
        dfdtau = -complex(np.log((1+np.sqrt(1-1/tau)/(1-np.sqrt(1-1/tau)))), -np.pi)/(2.0*np.sqrt(tau*(tau-1.0)))

    ftau = -Integral3(0.0, tau, 1/4.0)*tau/2.0
    B2SMA = 3.0*(ftau - 2.0*dfdtau*tau)/tau
    return B2SMA

def B3SMA(tau):
    # correction dgb in SM
    ftau = -Integral3(0.0, tau,  1/4.0)*tau/2.0
    B3SMA = 3.0*ftau/tau
    return B3SMA

def c_virt():
    tauc = mh2/mc2/4.0
    taub = mh2/mb2/4.0
    taut = mh2/mt2/4.0

    sum1 = gc*B1SMA(tauc) + gb*(B1SMA(taub) - B2SMA(taub)*dmb[0] - B3SMA(taub)*dgb) + gt*B1SMA(taut)
    sum2 = AMPLO(mh2)

    # susy contribution, SM parameters are fixed for QCD 
    # and SUSY parameters could be changed with these fixed SM parameters.
    # stored SQCD parameters as an array. 
    # csusy, error  = sumsusy('summaryfile.csv')

    c_virt = np.real(sum1/sum2)
    #return c_virt. np.real(csusy), error
    return c_virt
"""
def sumsusy(filename, norder):
    AMPQ, AMPQe = 0.0,0.0
    params = getdata(filename)

    retCt = params.reCt
    reCterr params.reCterr
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

    mst12 = mst1*mst1
    mst22 = mst2*mst2
    msb12 = msb1*msb1
    msb22 = msb2*msb2
    tanb2 = tanb*tanb

    alphaqst = params.alphaqst
    alphaqsb = params.alphaqsb

    gb = tanb
    gt = 1.0/tanb

    CF = 4.0/3.0
    deltab = (CF/2.0)*(alphaqsb/np.pi)*mG*mu*tanb*Ifun(msb12, msb22, mG**2)
    deltat = (CF/2.0)*(alphaqst/np.pi)*mG*mu*(1/tanb)*Ifunc(mst12, mst22, mG**2)

    if(norder == 0):
        Ct = complex(0.0,0.0)
        Ct_e = complex(0.0,0.0)
        Cb = complex(0.0,0.0)
        Cb_e = complex(0.0,0.0)
        CbLE = 0.0
        CtLE = 0.0

    elif(norder == 1):
        Ct = complex(reCt, imCt)
        Ct_e = complex(reCterr, imCterr)
        Cb = complex(reCb, imCb)
        Cb_e = complex(reCberr, imCberr)
        CbLE = deltab*(1.0 - 1.0/tanb**2)*np.pi/alphaqsb
        CtLE = deltat*(1.0 - tanb2)*np.pi/alphaqst

    A_bv = AMPLO(mb, mA)/3.0
    A_tb = AMPLO(mt, mA)/3.0

	gbeff = gb*(1.0 - deltab/tanb**2)/(1 + deltab)
	gteff = gt*(1.0 - deltat*tanb**2)/(1 + deltat)

    # need a better coding for tanbresum, tantresum
    if (tantresum == 1 and tanbresum == 1):
			Atq = A_tv*gteff
			Abq = A_bv*gbeff
		elif (tantresum == 0 and tanbresum == 0):
			Atq = A_tv*gt
			Abq = A_bv*gb
			CtLE = 0.0
			CbLE = 0.0
		elif (tantresum == 1 and tanbresum == 0):
			Atq = A_tv*gteff
			Abq = A_bv*gb
			CbLE = 0.0
		elif (tantresum == 0 and tanbresum == 1):
			Atq = A_tv*gt
			Abq = A_bv*gbeff
			CtLE = 0.0

	c_susy = np.abs(Atq+Abq)**2 + (Atq + Abq).conjugate()*(gt*A_tv*(Ct-CtLE) + gb*A_bv*(Cb - CbLE))*(alphas/np.pi)

	csusy_eb = 2.0*(Atq + Abq).conjugate()*(gb*A_bv*Cb_e)*(alphas/np.pi)
	csusy_et = 2.0*(Atq + Abq).conjugate()*(gt*A_tv*Ct_e)*(alphas/np.pi)
	error = np.sqrt(np.abs(csusy_eb)**2 + np.abs(csusy_et)**2)
	
	return c_susy, error
"""