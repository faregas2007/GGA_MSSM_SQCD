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

    def fsushi(self,
        m1:float=0.0, 
        m2:float=1.0,
        mb:float=einital().initial()['mbos'],
        mgl:float=inputs().M3,
        Ab:float=inputs().Ab,
        mu_in:float=einital().initial()['muRggh'],
        beta: float=np.arctan(inputs().tanb),
        muSUSY:float=inputs().mu
        ):

        tanb = np.tan(beta)
        lb1 = np.log(mu_in**2/m1**2)
        vert = -(2/3.0)*(-mgl/(Ab + muSUSY*(1/tanb))*b0fin(m1**2, mb, mgl,mu_in))
        waver = (2/3.0)*(m1**2/(m1**2 - m2**2))*(4 + 2*lb1 - (1 - mgl**2/m1**2 - mb**2/m1**2)*b0fin(m1**2, mb, mgl,mu_in))
        return (vert + waver)

    def dAbos(self,
        msb12:float=1.0,
        msb22:float=1.5,
        sthetab:float=0.0,
        cthetab:float=0.0,
        dmsb1_1:float=0.0,
        dmsb2_1:float=0.0,
        dthetab:float=0.0,
        mu_in: float=einital().initial()['muRggh'],
        beta: float=np.arctan(inputs().tanb),
        mb: float= einital().initial()['mbos'],
        mgl: float=inputs().M3,
        Ab: float=inputs().Ab,
        muSUSY: float=inputs().mu)->float:

        msb1 = np.sqrt(msb12)
        msb2 = np.sqrt(msb22)
        tanb = self.tanb
        sinb = np.sin(beta)
        cosb = np.cos(beta)

        dmsb12 = 2.0*dmsb1_1*msb1
        dmsb22 = 2.0*dmsb2_1*msb2
        lb1 = np.log(mu_in**2/msb12)
        lb2 = np.log(mu_in**2/msb22)
        mb2 = mb**2
        mu_in2 = mu_in**2
        mgl2 = mgl**2

        s2thetab = 2*sthetab*cthetab
        c2thetab = cthetab**2-sthetab**2
        f12 = self.fsushi(m1=msb1, 
                m2=msb2,
                mb=mb,
                mgl=mgl,
                Ab=Ab,
                mu_in=mu_in,
                muSUSY=muSUSY)

        f21 = self.fsushi(m1=msb2, 
                m2=msb1,
                mb=mb,
                mgl=mgl,
                Ab=Ab,
                mu_in=mu_in,
                muSUSY=muSUSY)

        dAbOS = (Ab + muSUSY*(cosb/sinb))*(-(msb12 - msb22)*s2thetab*(f12 + f21) + (dmsb12 - dmsb22)*s2thetab + 2*dthetab*(msb12 - msb22)*c2thetab)/(2*(Ab*mb - (msb12 - msb22)*sinb*cosb) + mb*muSUSY*cosb/sinb)
        return dAbOS

    def dAbOSdep(self,
        msb12:float=1.0,
        msb22:float=1.5,
        sthetab:float=0.0,
        cthetab:float=0.0,
        dmsb1_1:float=0.0,
        dmsb2_1:float=0.0,
        dthetab:float=0.0,
        mu_in: float=einital().initial()['muRggh'],
        beta: float=np.arctan(inputs().tanb),
        mb: float= einital().initial()['mbos'],
        mgl: float=inputs().M3,
        Ab: float=inputs().Ab,
        muSUSY: float=inputs().mu)->float:

        msb1 = np.sqrt(msb12)
        msb2 = np.sqrt(msb22)
        tanb = self.tanb
        sinb = np.sin(beta)
        cosb = np.cos(beta)

        dmsb12 = 2.0*dmsb1_1*msb1
        dmsb22 = 2.0*dmsb2_1*msb2
        lb1 = np.log(mu_in**2/msb12)
        lb2 = np.log(mu_in**2/msb22)
        mb2 = mb**2
        mu_in2 = mu_in**2
        mgl2 = mgl**2

        s2thetab = 2*sthetab*cthetab
        c2thetab = cthetab**2-sthetab**2
        f12 = self.fsushi(m1=msb1, 
                m2=msb2,
                mb=mb,
                mgl=mgl,
                Ab=Ab,
                mu_in=mu_in,
                muSUSY=muSUSY)

        f21 = self.fsushi(m1=msb2, 
                m2=msb1,
                mb=mb,
                mgl=mgl,
                Ab=Ab,
                mu_in=mu_in,
                muSUSY=muSUSY)

        dAbOS = - (Ab - muSUSY*tanb)/(muSUSY*sinb*cosb)*(Ab + muSUSY/tanb)*((c2thetab/s2thetab)*dthetab + (dmsb1_1 - dmsb2_1)/(msb1**2 - msb2**2) + (f21 + f12))

        return dAbOS
