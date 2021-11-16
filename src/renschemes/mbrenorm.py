# mbrenorm.py
# renormalization of bootm quark mass
# MS, DR and deltab resummation
# not yet included 2-loop electroweak effects into deltab

import sys

sys.path.insert(0, '/users/tp/dat/aggrenorm2/src/')

from typing import Dict, List

from init import *
from consts import *
from utils import *
from interface import *

class mbrenorm(einital):

    def dmbnontanb(self, 
        msb12: float=1.0,
        msb22: float=1.5,
        mu: float=einital().initial()['muRggh'],
        mb: float= einital().initial()['mbos'],
        mgl: float=inputs().M3,
        smflag: bool=True, 
        drflag: bool=False,   
        Ab: float=inputs().Ab)->float:
        
        mgl2 = mgl*mgl
        mb2 = mb**2
        lb = np.log(mu**2/mb2)
        lb1 = np.log(mu**2/msb12)
        lb2 = np.log(mu**2/msb22)
        lg = np.log(mu**2/mgl2)

        dmbnontanb = (- 2.0*(1.0+lg)*mgl**2+ msb12 + lb1*msb12 + msb22 + lb2*msb22 +(mb2 +mgl2- msb12 - 4.0*mb2*mgl*Ab/(msb12-msb22))*b0fin(mb2,mgl,np.sqrt(msb12),mu) + (mb2 +mgl2- msb22 + 4*mb2*mgl*Ab/(msb12-msb22))*b0fin(mb2,mgl,np.sqrt(msb22),mu))/(6.0*mb)

        if(smflag == True):
            if(drflag == True):
                dmbnontanb = dmbnontanb - 2.0*(5.0 + 3.0*lb)*mb2/(6.0*mb)
            else:
                dmbnontanb = dmbnontanb - 2.0*(4.0 + 3.0*lb)*mb2/(6.0*mb)
        
        return dmbnontanb
        
    def dmbtanb(self, 
        mb: float= einital().initial()['mbos'],
        msb12: float=1.0, 
        msb22: float=1.5,
        mu: float=einital().initial()['muRggh'],
        mgl: float=inputs().M3,
        beta: float=np.arctan(inputs().tanb),
        muSUSY: float=inputs().mu,
        )->float:

        tanb = self.tanb
        mb2 = mb**2
        
        dmbtanb = (4*mb2*mgl*muSUSY*tanb)/(msb12-msb22)*b0fin(mb2,mgl,np.sqrt(msb12),mu) - 4*mb2*mgl*muSUSY*tanb/(msb12-msb22)*b0fin(mb2,mgl,np.sqrt(msb22),mu)/(6.0*mb)
        return dmbtanb
    
    def mbMSDRtrans(self, 
        mb: float= einital().initial()['mbos'],
        norder: float=inputs().norder, 
        mu_in: float=einital().initial()['muRggh'])->float:

        alphas = SUSHI_alphas(inputs().alphasmz, mu_in, inputs().mZ)
        mbMSDRtrans = mb*(1.0 - 1.0/3.0 * alphas/Pi)
        if(norder == 1):
            mbMSDRtrans = mbMSDRtrans - mb*(29.0/72.0 * (alphas/Pi)**2)
        return mbMSDRtrans

    def mbMSbarmuX(self, 
        mbmb: float=inputs().mbmb,
        mu_in: float=einital().initial()['muRggh']
        ):
        alphas1 = SUSHI_alphas(inputs().alphasmz, inputs().mbmb, inputs().mZ)/Pi
        alphas2 = SUSHI_alphas(inputs().alphasmz, mu_in, inputs().mZ)/Pi
        return runmass(mbmb, alphas1, alphas2, nf, 4)

    # to modify with 2-loop electroweak effects in deltab term
    # for now, use mu_in and alphasren as input
    # latter use 1 of them. 
    def deltab(self, 
        msb12: float=1.0, 
        msb22: float=1.5, 
        mgl: float=inputs().M3,
        Ab: float=inputs().Ab,
        muSUSY: float=inputs().mu, 
        mu_in: float=einital().initial()['muRggh'],
        alphasren: float=SUSHI_alphas(inputs().alphasmz,einital().initial()['muRggh'],inputs().mZ))->Dict:
        delmb = ifunc(msb12, msb22, mgl*mgl)/(2.0*Pi)
        delAb = -cf*alphasren*mgl*Ab*delmb
        delmb = cf*alphasren*mgl*muSUSY*self.tanb*delmb
        return {'delmb': delmb, 'delAb': delAb}
