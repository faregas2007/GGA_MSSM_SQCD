# sbottomrenorm.py
# renormalization in bottom/sbottom sector
# from sushi

import sys

from typing import Dict, List

from init import *
from consts import *
from utils import *

from renschemes.Abrenorm import *

class sbottomrenorm(einital):

    def sbottomrenormHM(self,  
        mbsb: float=1.0, 
        msb12: float=1.0, 
        msb22: float=1.5, 
        sthetab: float=0.0, 
        cthetab: float=0.0, 
        mu_in: float=einital().initial()['muRggh'],
        mgl: float=inputs().MQ3,
        mb: float= einital().initial()['mbos'])->Dict:
        """
        The inputs could varied, LO-values or corrected parameters.
        Keep them seperate from inptus and einitial parameters.
        """
        s2thetab = 2.0*sthetab*cthetab
        c2thetab = cthetab**2 - sthetab**2
        msbot1 = np.sqrt(msb12)
        msbot2 = np.sqrt(msb22)
        mb2 = mb*mb

        dmsb11 = 0.0
        dmsb12 = 0.0
        dmsb21 = 0.0
        dmsb22 = 0.0
        dthetab1 = 0.0
        dthetab2 = 0.0
        dhbhb1 = 0.0
        dhbhb2 = 0.0
        dhbhb3 = 0.0
 
        dmsb11 = self.dmsb1osfin(mb2=mb2,mgl=mgl,msb12=msb12,msb22=msb22,sthetab=sthetab,cthetab=cthetab,mu_in=mu_in)/msbot1
        dmsb21 = self.dmsb2osfin(mb2=mb2,mgl=mgl,msb12=msb12,msb22=msb22,sthetab=sthetab,cthetab=cthetab,mu_in=mu_in)/msbot2
        dthetab1 = self.dthetabosfin(mb2=mb2,mgl=mgl,msb12=msb12,msb22=msb22,s2thetab=s2thetab,c2thetab=c2thetab,mu_in=mu_in)
        dhbhb2 = self.dmbfin(True,True, mb=mbsb,mgl=mgl,msb12=msb12,msb22=msb22,s2thetab=s2thetab,mu_in=mu_in)

        return {'dmsb11':dmsb11, 'dmsb12':dmsb12, 'dmsb21':dmsb21, 'dmsb22':dmsb22, 'dthetab1':dthetab1,'dthetab2':dthetab2, 'dhbhb1':dhbhb1, 'dhbhb2':dhbhb2, 'dhbhb3':dhbhb3}

    def sbottomrenormDS(self, 
        thetab: float=0.0,
        msb12: float=1.0,
        msb22: float=1.5,
        mu: float=einital().initial()['muRggh'],
        mb: float= einital().initial()['mbos'],
        mgl: float=inputs().MQ3)->Dict:

        cthetab2 = np.cos(2.0*thetab)
        sthetab2 = np.sin(2.0*thetab)
        mb2 = mb*mb
        dhbhb1, dhbhb2, dhbhb3 = 0.0,0.0,0.0
        dmsb11, dmsb12, dmsb22 = 0.0,0.0,0.0
        dthetab1, dthetab2 = 0.0, 0.0

        dmsb11 = msb12 * 1/3.0 *(3.0 * np.log(msb12/mu**2) -3.0 - cthetab2**2 * (np.log(msb12/mu**2) - 1.0) -sthetab2**2 * msb22/msb12 * (np.log(msb22/mu**2) - 1.0) -6.0 * mgl**2/msb12 -2.0 * (1.0 - 2.0 * mgl**2/msb12)*np.log(mgl**2/mu**2) -2.0 * (1.0 - mgl**2/msb12)**2 * np.log(np.abs(1.0 - msb12/mgl**2))  )
        dmsb12 = (-4*sthetab2*mgl*(-mgl**2*np.log(np.abs(1.0 - msb12/mgl**2)) + msb12*(-2 + np.log(np.abs(mgl**2 - msb12)/mu**2)))*np.sqrt(mb2))/(3.0*msb12)
        dmsb21 = msb22 * 1/3.0 * (3.0 * np.log(msb22/mu**2) - 3.0 - cthetab2**2 * (np.log(msb22/mu**2) - 1.0) - sthetab2**2 * msb12/msb22 * (np.log(msb12/mu**2) - 1.0) - 6.0 * mgl**2/msb22 - 2.0 * (1.0 - 2.0 * mgl**2/msb22)*np.log(mgl**2/mu**2)- 2.0 * (1.0 - mgl**2/msb22)**2 * np.log(np.abs(1.0 - msb22/mgl**2))  )
        dmsb22 = (4*sthetab2*mgl*(-mgl**2*np.log(np.abs(1.0 - msb22/mgl**2)) + msb22*(-2 + np.log(np.abs(mgl**2 - msb22)/mu**2)))*np.sqrt(mb2))/(3.0*msb22)
        dthetab1 = sthetab2 * 1/3.0 *(-2.0 * cthetab2**2 + 2.0 * cthetab2**2 /(msb12 - msb22)* (msb12 * np.log(msb12/mu**2) - msb22 * np.log(msb22/mu**2)))
        dthetab2 = (4*mgl*cthetab2**2*(mgl**2*msb22*np.log(np.abs(1.0 - msb12/mgl**2)) + msb12*(mgl**2*np.log(np.abs(1.0 - msb22/mgl**2)) + msb22*(4 - np.log(np.abs((mgl**2 - msb12)*(mgl**2 - msb22))/mu**4))))*np.sqrt(mb2))/(3.0*msb12*(msb12 - msb22)*msb22)

        x1 = msb12/mgl**2
        x2 = msb22/mgl**2

        #dhbhb_pole from (32) and (29)
        dhbhb1 = (sthetab2*mgl*((x1*np.log(x1))/(-1 + x1) - (x2*np.log(x2))/(-1 + x2)))/(3.*np.sqrt(mb2))
        dhbhb2 = np.log(mb2/mu**2) + (-5 + (-3 + x1)/(4.*(-1 + x1)) + (-3 + x2)/(4.*(-1 + x2)) - np.log(mgl**2/mu**2) - ((-2 + x1)*x1*np.log(x1))/(2.*(-1 + x1)**2) - ((-2 + x2)*x2*np.log(x2))/(2.*(-1 + x2)**2))/3.0
        dhbhb3 = sthetab2 * (((-3 + x1)*x1)/(2.*(-1 + x1)**2) - ((-3 + x2)*x2)/(2.*(-1 + x2)**2) + (x1*np.log(x1))/(-1 + x1)**3 - (x2*np.log(x2))/(-1 + x2)**3 ) * np.sqrt(mb2)/(3.0*mgl)

        dhbhb1 = dhbhb1*mb
        dhbhb2 = dhbhb2*mb
        dhbhb3 = dhbhb3*mb
        dmsb11 = dmsb11/2.0/msb12
        dmsb12 = dmsb12/2.0/msb12
        dmsb21 = dmsb21/2.0/msb22
        dmsb22 = dmsb22/2.0/msb22
        dthetab1 = dthetab1/2.0/np.cos(2*thetab)
        dthetab2 = dthetab2/2.0/np.cos(2*thetab)

        return {'dmsb11':dmsb11, 'dmsb12':dmsb12,'dmsb21':dmsb21, 
            'dmsb22':dmsb22, 'dthetab1':dthetab1, 'dthetab2':dthetab2,
            'dhbhb1':dhbhb1, 'dhbhb2':dhbhb2, 'dhbhb3':dhbhb3}

    def dmsb1osfin(self,
        msb12: float=1.0,
        msb22: float=1.5,
        sthetab: float=0.0,
        cthetab: float=0.0,
        mu_in: float=1000.0,
        mb2: float= einital().initial()['mbos']*einital().initial()['mbos'],
        mgl: float=inputs().M3)->float:
        
        lb1 = np.log(mu_in**2/msb12)
        lb2 = np.log(mu_in**2/msb22)
        lg = np.log(mu_in**2/mgl**2)
        
        mgl2 = mgl*mgl
        if (mb2 != 0.0):
            lb = np.log(mu_in**2/mb2)
        else:
            lb = 0.0     
        

        s2thetab = 2*sthetab*cthetab     
        dmsb1osfin =  -((1 + lb)*mb2 + (1 + lg)*mgl**2 + (3 + lb1)*msb12 + 2*cthetab**2*((1 + lb1)*msb12 - (1 + lb2)*msb22)*sthetab**2 + (mb2 +mgl2- msb12 - 2*np.sqrt(mb2)*mgl*s2thetab)*b0fin(msb12,np.sqrt(mb2),mgl,mu_in))/(3.0*np.sqrt(msb12))      
        
        return dmsb1osfin        

    def dmsb2osfin(self,
        msb12:float=1.0,
        msb22:float=1.5,
        sthetab: float=0.0,
        cthetab: float=0.0,
        mu_in: float=1000.0,
        mb2: float= einital().initial()['mbos']*einital().initial()['mbos'],
        mgl: float=inputs().M3)->float:
        mgl2 = mgl*mgl

        if (mb2 != 0.0):
            lb = np.log(mu_in**2/mb2)
        else:
            # lb always accompany with mb -> limit mb log(mb/mu)-> 0
            lb = 0.     
        lb1 = np.log(mu_in**2/msb12)
        lb2 = np.log(mu_in**2/msb22)
        lg = np.log(mu_in**2/mgl**2)     
        s2thetab = 2*sthetab*cthetab     
        dmsb2osfin = -((1 + lb)*mb2 + (1 + lg)*mgl**2 + (3 + lb2)*msb22 + 2*cthetab**2*(-((1 + lb1)*msb12) + (1 + lb2)*msb22)*sthetab**2 + (mb2 +mgl2- msb22 + 2*np.sqrt(mb2)*mgl*s2thetab)*b0fin(msb22,np.sqrt(mb2),mgl,mu_in))/(3.0*np.sqrt(msb22))
        return dmsb2osfin

    def dthetabosfin(self,
        msb12: float=1.0, 
        msb22: float=1.5,
        s2thetab: float=0.0,
        c2thetab: float=0.0,
        mu_in: float=1000.0,
        mb2:float= einital().initial()['mbos']*einital().initial()['mbos'],
        mgl: float=inputs().M3)->float:
        
        lb1 = np.log(mu_in**2/msb12)
        lb2 = np.log(mu_in**2/msb22)

        mgl2 = mgl**2
        mu_in2 = mu_in**2

        dthetabosfin = (c2thetab*((-((1 + lb1)*msb12)
                + (1 + lb2)*msb22)*s2thetab
                + 2.0*np.sqrt(mb2)*mgl*(b0fin(msb12, np.sqrt(mb2), mgl, mu_in)
                + b0fin(msb22, np.sqrt(mb2), mgl, mu_in))))/(3.0*(msb12 - msb22))

        return dthetabosfin
    
    # should belongs to mbrenorm --> todo.
    def dmbfin(self,
        smflag: bool=True, 
        drflag: bool=False, 
        mb: float= einital().initial()['mbos'],
        mgl:float=inputs().M3,
        msb12:float=1.0,
        msb22:float=1.5,
        s2thetab:float=0.0, 
        mu_in: float=1000.0)->float:

        mb2 = mb*mb
        if (mb2 != 0.0):
            lb = np.log(mu_in**2/mb2)
        else:
            lb = 0.0

        lb1 = np.log(mu_in**2/msb12)
        lb2 = np.log(mu_in**2/msb22)
        lg = np.log(mu_in**2/mgl**2)

        mgl2 = mgl**2
        mu2 = mu_in**2
        
        msb1 = np.sqrt(msb12)
        msb2 = np.sqrt(msb22)

        if (mb2 != 0.0):
            dmbfin = (-2.0*(1.0 + lg)*mgl**2
                + msb12 + lb1*msb12 + msb22 + lb2*msb22
                + (mb2 +mgl2- msb12 - 2.0*mb*mgl*s2thetab)*b0fin(mb2,mgl,msb1,mu_in)
                + (mb2 +mgl2- msb22
                + 2*mb*mgl*s2thetab)*b0fin(mb2,mgl,msb2,mu_in))/(6.0*mb)
        else:
            dmbfin = 0.0

        if (smflag == True):
            if (drflag == True):
                if (mb2 != 0.0):
                    dmbfin = dmbfin - 2.0*(5.0 + 3.0*lb)*mb2/(6.0*mb)
            else:
                if (mb2 != 0.0):
                    dmbfin = dmbfin - 2.0*(4.0 + 3.0*lb)*mb2/(6.0*mb)

        return dmbfin
    
    # mb in sbottom sector
    def mb_sbottom(self,
        mbsb_var1: float=4.0,
        mbsb_var2: float=0.0, 
        msb12: float=1.0, 
        msb22: float=1.5, 
        dmsb11: float=0.0, 
        dmsb21: float=0.0,
        s2thetab: float=0.0,
        c2thetab: float=0.0,
        dthetabos1: float=0.0,
        mb: float= einital().initial()['mbos'],
        mgl: float=inputs().M3, 
        beta: float=np.arctan(inputs().tanb), 
        muSUSY: float=inputs().mu,
        muD: float=einital().initial()['muRggh'],
        prefac: float=SUSHI_alphas(inputs().alphasmz, 1000.0, inputs().mZ)/Pi)->float:
        
        mb2 = mb*mb
        mgl2 = mgl*mgl
        lb1 = np.log(muD**2/msb12)
        lb2 = np.log(muD**2/msb22)
        thetab = np.arcsin(s2thetab)/2.0
        dAbosfintbHM = Abrenorm().dAbosfintbHM(
            msb12=msb12,
            msb22=msb22,
            sthetab=np.sin(thetab),
            cthetab=np.cos(thetab),
            dmsb1_1=dmsb11,
            dmsb2_1=dmsb21,
            deltab=0.0,
            mu_in=muD,
            beta=beta,
            mb=mb,
            mgl=mgl,
            Ab=self.Ab, 
            muSUSY=muSUSY
        )
        if(np.abs(s2thetab) > 10**-8):
            mbsb_var1 = mbsb_var2 - (2.0*mb/(msb12-msb22)*(dmsb11*msb12 - dmsb21*msb22)
                - (c2thetab*(msb12-msb22)*np.cos(beta)*np.sin(beta))/muSUSY*dthetabos1
                -2.0*mb2/(msb12-msb22)/s2thetab*dAbosfintbHM)*prefac
        else:
            # limit for sin(2*thetab)->0
            mbsb_var1 = mbsb_var2 - (2.0*mb/(msb12-msb22)*(dmsb11*msb12 - dmsb21*msb22)
                +(mb*((-6*msb12*dmsb11
                + 6*msb22*dmsb21
                - 8*msb12 - 4*lb1*msb12 + 8*msb22
                + 4*lb2*msb22)*muSUSY
                + 2*b0fin(msb22, mb, mgl, muD)*((mb2 +mgl2- msb22)*muSUSY
                + mgl*(msb12 - msb22)*np.cos(beta)*np.sin(beta))
                + b0fin(msb12, mb, mgl, muD)*(-2*(mb2 +mgl2- msb12)*muSUSY
                + mgl*(msb12 - msb22)*np.sin(2*beta))))/(3.0*(msb12-msb22)*muSUSY))*prefac
            if(thetab < 1.0):
            # thetab -> 0
                mbsb_var1 = mbsb_var1 + ((msb12-msb22)*np.cos(beta)*np.sin(beta))/muSUSY*dthetabos1*prefac
            else:
            # thetab-> pi/2
                mbsb_var1 = mbsb_var1 - ((msb12-msb22)*(np.sin(beta)*np.cos(beta)))/muSUSY*dthetabos1*prefac
        
        return mbsb_var1

    def get_json(self, obj: str='sbottomrenormHM')->Dict:
        if(obj=='sbottomrenormDS'):
            return json.dumps(self.sbottomrenormDS(), cls=NumpyEncoder, indent=4)
        elif(obj=='sbottomrenormHM'):
            return json.dumps(self.sbottomrenormHM(), cls=NumpyEncoder, indent=4)
    