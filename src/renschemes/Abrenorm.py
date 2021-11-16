# Abrenorm.py
# renormalization of Ab, finite part only from
# HM: Harlander, Mantek
# DS: Degrassi, Slavich

import sys

from typing import Dict, List

from init import *
from consts import *
from utils import *
from interface import *

class Abrenorm(einital):

    # all parameters could vary through LO and NLO values. 
    def dAbosfinmbHM(self,
        msb12: float=1.0,
        msb22: float=1.5,
        deltamb: float=0.0,
        mu_in: float=einital().initial()['muRggh'],
        beta: float=np.arctan(inputs().tanb),
        mb: float= einital().initial()['mbos'],
        mgl: float=inputs().M3,
        Ab: float=inputs().Ab, 
        muSUSY: float=inputs().mu,
        )->float:
    
        tanb = self.tanb
        mgl2 = mgl*mgl

        lb1 = np.log(mu_in**2/msb12)
        lb2 = np.log(mu_in**2/msb22)
        mb2 = mb**2

        dAbosfinmbHM = -((4*lb1*mb*msb12 + (3*deltamb+ 8*mb)*(msb12 - msb22) 
            - 4*lb2*mb*msb22)*(Ab + muSUSY/tanb) 
            + 2*mb*(Ab*(mb2 +mgl2- msb12) + mgl*(-msb12 + msb22) 
            + ((mb2 +mgl2- msb12)*muSUSY)/tanb)*b0fin(msb12,mb,mgl,mu_in) 
            - 2*mb*(Ab*(mb2 +mgl2- msb22) + mgl*(msb12 - msb22) 
            + ((mb2 +mgl2- msb22)*muSUSY)/tanb)*b0fin(msb22,mb,mgl,mu_in))/(3.0*mb*(msb12 - msb22))

        return dAbosfinmbHM

    def dAbosfintbHM(self,
        msb12:float=1.0,
        msb22:float=1.5,
        sthetab:float=0.0,
        cthetab:float=0.0,
        dmsb1_1:float=0.0,
        dmsb2_1:float=0.0,
        deltab:float=0.0,
        mu_in: float=einital().initial()['muRggh'],
        beta: float=np.arctan(inputs().tanb),
        mb: float= einital().initial()['mbos'],
        mgl: float=inputs().M3,
        Ab: float=inputs().Ab,
        muSUSY: float=inputs().mu)->float:

        tanb = self.tanb

        dmsb12 = 2.0*dmsb1_1*msb12
        dmsb22 = 2.0*dmsb2_1*msb22
        lb1 = np.log(mu_in**2/msb12)
        lb2 = np.log(mu_in**2/msb22)
        mb2 = mb**2
        mu_in2 = mu_in**2
        mgl2 = mgl**2

        s2thetab = 2*sthetab*cthetab
        c2thetab = cthetab**2-sthetab**2

        dAbosfintbHM =  ((6*c2thetab*deltab*(msb12 - msb22) +
            (3*dmsb12 - 3*dmsb22 + 8*msb12 + 4*lb1*msb12 - 8*msb22 -
            4*lb2*msb22)*s2thetab)*(Ab + muSUSY/tanb) -
            (2*s2thetab*b0fin(msb22,mb,mgl,mu_in)*((mb2 +mgl2- msb22)*muSUSY*np.cos(beta) +
            (Ab*(mb2 +mgl2- msb22) + mgl*(msb12 - msb22))*np.sin(beta)))/np.sin(beta) 
            + (2*s2thetab*b0fin(msb12,mb,mgl,mu_in)*((mb2 +mgl2- msb12)*muSUSY*np.cos(beta) +
            (Ab*(mb2 +mgl2- msb12) + mgl*(-msb12 + msb22))*np.sin(beta)))/np.sin(beta))/(6*Ab*mb - 3*(msb12 - msb22)*s2thetab + (6*mb*muSUSY)/tanb)
        
        return dAbosfintbHM