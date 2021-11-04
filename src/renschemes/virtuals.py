# virtuals.py

import numpy as np
import pandas as pd
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
        self.params_fp = params_fp
        self.params = Namespace(**load_dict(Path(self.params_fp, 'models.json')))
        self.cqcd = pd.read_csv(Path(data_dir, "CQCD.csv"))

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
        B1SMA = 3/2.0*complex(self.ReF0AB1A(x),self.ImF0AB1A(x))
        return B1SMA
    
    
    def B2SMA(self, tau):
        if (tau <= 1.0):
            dfdtau = complex(np.arcsin(np.sqrt(tau)/np.sqrt(tau*(1-tau))))
        else:
            dfdtau = -complex(np.log((1+np.sqrt(1-1/tau)/(1-np.sqrt(1-1/tau)))), -Pi)/(2.0*np.sqrt(tau*(tau-1.0)))

        ftau = -Integral3(0.0, tau, 1/4.0)*tau/2.0
        B2SMA = 3.0*(ftau - 2.0*dfdtau*tau)/tau
        return B2SMA
    
    def B3SMA(self, tau):
    # correction dgb in SM
        ftau = -Integral3(0.0, tau,  1/4.0)*tau/2.0
        B3SMA = 3.0*ftau/tau
        return B3SMA
    
    def c_virt(self):
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