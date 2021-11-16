# virtuals.py

import numpy as np
import pandas as pd
from argparse import Namespace
from pylooptools.ltools import *

from consts import *
from init import *
from utils import *
from config import *

from xs.integrals import *
from xs.reals import *
from renschemes.renschemes import *

class virtual(einital): 
    def __init__(self, params_fp:str=Path(model_dir), *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.params_fp = params_fp
        self.params = Namespace(**load_dict(Path(self.params_fp, 'models.json')))
        self.cqcd = pd.read_csv(Path(data_dir, "CQCD.csv"))

        self.mG = self.params.mgl
        self.norder = self.params.norder
        self.mb = self.params.mb
        self.mt = self.params.mt

    def interpolate(self, df, values, colname1, colname2):
        df = df.iloc[(df[colname1]-values).abs().argsort()[:2]]
        interpolate =  (df[:1][colname2].values - df[1:2][colname2].values)/(df[:1][colname1].values - df[1:2][colname1].values) + df[1:2][colname2].values
        return interpolate[0]

    def ImF0AB1A(self, tau):
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
            ImF0AB1A =  self.interpolate(self.cqcd, tau, 'tauimval', 'imvalues')

        return ImF0AB1A

    
    def ReF0AB1A(self, tau):
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
                ReF0AB1A = self.interpolate(self.cqcd, tau, 'taureval', 'revalues')  
        return ReF0AB1A

    def B1SMA(self, x):
        print(f'x: {x}')
        if(np.array(x).shape == ()):
            B1SMA = 3/2.0*complex(self.ReF0AB1A(x),self.ImF0AB1A(x))
        else:
            B1SMA = [3/2.0*complex(self.ReF0AB1A(ix), self.ImF0AB1A(ix)) for ix in x]
        return np.array(B1SMA)
    
    
    def B2SMA(self, tau):
        if(np.array(tau).shape == ()):
            if (tau <= 1.0):
                dfdtau = complex(np.arcsin(np.sqrt(tau)/np.sqrt(tau*(1-tau))))
            else:
                dfdtau = -complex(np.log((1+np.sqrt(1-1/tau)/(1-np.sqrt(1-1/tau)))), -Pi)/(2.0*np.sqrt(tau*(tau-1.0)))

            ftau = -Integral3(0.0, tau, 1/4.0)*tau/2.0
            B2SMA = 3.0*(ftau - 2.0*dfdtau*tau)/tau
        else:
            temp = tau[ids == tau[tau<=1.0]]
            if(temp.size != 0):
                dfdtau = [complex(np.arcsin(np.sqrt(tau)/np.sqrt(itemp*(1-tau)))) for tau in temp]
            temp = tau[ids == tau[tau>1]]
            if(temp.size != 0):
                 dfdtau = [-complex(np.log((1+np.sqrt(1-1/tau)/(1-np.sqrt(1-1/tau)))), -Pi)/(2.0*np.sqrt(tau*(tau-1.0))) for tau in temp]
        
            B2SMA = 3.0*(ftau - 2.0*np.array(dfdtau)*tau)/tau
        return B2SMA
    
    def B3SMA(self, tau):
    # correction dgb in SM
        if(np.array(tau).shape == ()):
            ftau = -Integral3(0.0, tau,  1/4.0)*tau/2.0
            B3SMA = 3.0*ftau/tau
        else:
            ftau = -Integral3_batch(0.0, tau, 1/4.0)*tau/2.0
        
        return B3SMA
    
    def c_virts(self):
        dmb = self.params.dmb
        dgb = self.params.dgb
        
        mc  = self.params.mcos
        mb =  self.params.mb
        mt =  self.params.mt

        mc2 = mc*mc  
        mb2 = mb*mb
        mt2 = mt*mt
        mh2 = self.params.mh2

        gc = self.params.gc
        gb = self.params.gb
        gt = self.params.gt

        tauc = mh2/mc2/4.0
        taub = mh2/mb2/4.0
        taut = mh2/mt2/4.0

        sum1 = gc*self.B1SMA(tauc) + gb*(self.B1SMA(taub) - self.B2SMA(taub)*dmb[0] - self.B3SMA(taub)*dgb) + gt*self.B1SMA(taut)
        sum2 = real_amps().AMPLO(mh2)

        c_virt = np.real(sum1/sum2)
        return c_virt 

    def c_virts_batch(self):
        dmb = self.params.dmb
        dgb = self.params.dgb
        
        mc  = self.params.mcos
        mb =  self.params.mb
        mt =  self.params.mt

        mc2 = mc*mc  
        mb2 = mb*mb
        mt2 = mt*mt
        mh2 = self.params.mh2

        gc = self.params.gc
        gb = self.params.gb
        gt = self.params.gt

        tauc = mh2/mc2/4.0
        taub = mh2/mb2/4.0
        taut = mh2/mt2/4.0

        sum1 = gc*self.B1SMA(tauc) + gb*(self.B1SMA(taub) - self.B2SMA(taub)*dmb[0] - self.B3SMA(taub)*dgb) + gt*self.B1SMA(taut)
        sum2 = real_amps().AMPLO_batch(mh2)

        c_virt = np.real(sum1/sum2)
        return c_virt
    
    def c_virts_sqcd(self, alphas):
        mG = self.mG
        norder = self.norder
        mb = self.mb
        mt = self.mt

        # fixed tanbresum and tanbresum0 for now
        tantresum = 0
        tanbresum0 = 1

        AMPQ, AMPQe = 0.0, 0.0
        params = getdata(Path(data_dir, 'CSQCD.csv'))

        reCt = params.reCt
        reCterr = params.reCterr
        imCt = params.imCt
        imCterr = params.imCterr
        reCb = params.reCb
        reCberr = params.reCberr
        imCb = params.imCb
        imCberr = params.imCberr

        mu = params.mu[10]
        At = params.At
        Ab = params.Ab

        mA = params.mA
        tanb = params.tanb

        mst1 = params.mst1
        mst2 = params.mst2
        msb1 = params.msb1
        msb2 = params.msb2

        gb = tanb
        gt = 1/tanb

        msb12 = msb1*msb1
        msb22 = msb2*msb2
        mst12 = mst1*mst1
        mst22 = mst2*mst2
        tanb2 = tanb*tanb

        #print(f'reCt, imCt: {reCt, imCt}')
        alphaqst = params.alphaqst[10]
        alphaqsb = alphaqst

        CF = 4.0/3.0
        c_susy = []
        error = []
        for idx  in range(len(mA)):
            deltab = (CF/2.0)*(alphaqsb/Pi)*mG*mu*tanb[idx]*ifunc(msb1[idx]**2, msb2[idx]**2, mG*mG)
            deltat = (CF/2.0)*(alphaqst/Pi)*mG*mu*(1/tanb[idx])*ifunc(mst1[idx]**2, mst2[idx]**2, mG*mG)
            
            if(norder == 0):
                Ct = complex(0.0,0.0)
                Ct_e = complex(0.0,0.0)
                Cb = complex(0.0,0.0)
                Cb_e = complex(0.0,0.0)
                CbLE = 0.0
                CtLE = 0.0
            elif(norder==1):
                Ct = complex(reCt[idx], imCt[idx])
                Ct_e = complex(reCterr[idx], imCterr[idx])
                Cb = complex(reCb[idx], imCb[idx])
                Cb_e = complex(reCberr[idx], imCberr[idx])

                CbLE = deltab*(1.0 - 1/tanb[idx]**2)*Pi/alphaqsb
                CtLE = deltat*(1.0 -tanb[idx]**2)*Pi/alphaqst

            # full Leading order terms without factor
            A_bv = Aqtau(mb, mA[idx])
            A_tv = Aqtau(mt, mA[idx])

            # have changd into AMPLO of Sushi
            # let make it work first
            #A_bv = AMPLO(mb, mA[idx])/3.0
            #A_tb = AMPLO(mt, mA[idx])/3.0

            gbeff = gb[idx]*(1.0 - deltab/tanb[idx]**2)/(1 + deltab)
            gteff = gt[idx]*(1.0 - deltat*tanb[idx]**2)/(1 + deltat)
            if(tantresum == 1 and tanbresum0 == 1):
                Atq = A_tv*gteff
                Abq = A_bv*gbeff
            elif(tantresum == 0 and tanbresum0 == 0):
                Atq = A_tv*gt[idx]
                Abq = A_bv*gb[idx]
                CtLE = 0.0
                CbLE = 0.0
            elif(tantresum == 1 and tanbresum0 == 0):
                Atq = A_tv*gteff
                Abq = A_bv*gb[idx]
                CbLE = 0.0
            elif(tantresum == 0 and tanbresum0 == 1):
                Atq = A_tv*gt[idx]
                Abq = A_bv*gbeff
                CtLE = 0.0

            #print(f'Atq, Abq: {Atq, Abq}')
            c_susy.append(np.abs(Atq+Abq)**2 + (Atq + Abq).conjugate()*(gt[idx]*A_tv*(Ct-CtLE) + gb[idx]*A_bv*(Cb-CbLE))*(alphas/Pi))
            csusy_eb = 2.0*(Atq + Abq).conjugate()*(gb[idx]*A_bv*Cb_e)*(alphas/Pi) 
            csusy_et = 2.0*(Atq + Abq).conjugate()*(gt[idx]*A_tv*Ct_e)*(alphas/Pi)
            error.append(np.sqrt(np.abs(csusy_eb)**2 + np.abs(csusy_et)**2))

        return np.array(c_susy), np.array(error)

def ftau(x):
    if(x>=1):
        ftau = (np.arcsin(1.0/np.sqrt(x)))**2
    else:
        ftau = -(np.log((1.0 + np.sqrt(1-x))/(1 - np.sqrt(1-x))) - complex(0.0,1.0)*Pi)**2/4.0
    return ftau

def Aqtau(mq, mhiggs):
    tau = (2*mq/mhiggs)
    Aqtau = tau*ftau(tau)
    return Aqtau