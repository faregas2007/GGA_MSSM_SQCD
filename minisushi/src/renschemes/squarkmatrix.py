# squarkmatrix.py
# Find egienvalues of squark masses and mixing angles

import sys
import numpy as np
#sys.path.insert(0, '/users/tp/dat/aggrenorm2/src/common/')

from typing import Dict, List


from init import *
from consts import *
from utils import *
from interface import *

#from src import init, consts, utils

class squarkmatrix(inputs):

    # (M3SU, M3SQ, M3SD) -> (MQ3, MTR, MBR)
    def calcsquarkmass(self, 
        M3SU: float = inputs().MQ3, 
        M3SQ: float = inputs().MTR, 
        M3SD: float = inputs().MBR, 
        mw: float= inputs().mW, 
        mz: float= inputs().mZ, 
        mt: float= inputs().mt, 
        mb: float= einital().initial()['mbos'], 
        muSUSY: float=inputs().mu, 
        beta: float=np.arctan(inputs().tanb), 
        At: float=inputs().At, 
        Ab: float=inputs().Ab, 
        dMLb: float=0.0)->Dict:

        tanb = self.tanb

        sthetaw = np.sqrt(1.0 - (mw/mz)**2)
        mb2 = mb*mb
        mt2 = mt*mt
        tanb = np.tan(beta) 
        M3SU2 = M3SU*M3SU
        M3SQ2 = M3SQ*M3SQ
        M3SD2 = M3SD*M3SD   
        if (M3SU < 0.0):
            M3SU2 = - M3SU2
        if (M3SQ < 0.0):
            M3SQ2 = - M3SQ2
        if (M3SD < 0.0):
            M3SD2 = - M3SD2 

        MstopEig = M3SQ2 + mt2 + mz**2*np.cos(2.0*beta)*(1/2.0 - 2/3.0*sthetaw**2)  
        # calculate mstop1 and mstop2 with mstop1 < mstop2
        mstop1=np.sqrt(M3SQ2/2.0 + M3SU2/2.0 + mt2 + mz**2/4.0*np.cos(2.0*beta) -np.sqrt((M3SQ2 - M3SU2 + mz**2*np.cos(2.0*beta)*(1/2.0-4/3.0*sthetaw**2))**2 + 4.0*mt2*(At-muSUSY/tanb)**2 )/2.0)
        mstop2=np.sqrt(M3SQ2/2.0 + M3SU2/2.0 + mt2 +mz**2/4.0*np.cos(2.0*beta) +np.sqrt((M3SQ2 - M3SU2 + mz**2 *np.cos(2.0*beta)*(1/2.0-4/3.0*sthetaw**2))**2+ 4.0*mt2*(At-muSUSY/tanb)**2 )/2.0)
        mst12 = mstop1**2
        mst22 = mstop2**2   
        thetat = np.arcsin(2.0*mt/(mst12-mst22)*(At-muSUSY/tanb))/2.0
        if (2.0*mt/(mst12-mst22)*(At-muSUSY/tanb)>1.0 and 2.0*mt/(mst12-mst22)*(At-muSUSY/tanb)<1.0001):
            thetat = np.arcsin(1.0)/2.0
        if (2.0*mt/(mst12-mst22)*(At-muSUSY/tanb) < -1.0 and 2.0*mt/(mst12-mst22)*(At-muSUSY/tanb)>-1.0001):
            thetat = np.arcsin(-1.0)/2.0    
        sthetat = np.sin(thetat)
        cthetat = np.cos(thetat)    
        if(np.abs(cthetat**2 *mst12 +sthetat**2 *mst22 - MstopEig)>10**-3):
            thetat = Pi/2.0 - thetat
            sthetat = np.sin(thetat)
            cthetat = np.cos(thetat)    
        c2thetat = cthetat**2 - sthetat**2
        s2thetat = 2*sthetat*cthetat    
        t2thetat = np.tan(2*thetat)

        
        MsbotEig = M3SQ2 + dMLb + mb2 + mz**2 * np.cos(2.0*beta) * (-1/2.0 + 1/3.0 * sthetaw**2)    
        #calculate msbot1 and msbot2 with msbot1 < msbot2
        msbot1=np.sqrt(M3SQ2/2.0 + M3SD2/2.0 + dMLb/2.0 + mb2 - mz**2/4.0*np.cos(2*beta)-np.sqrt((M3SQ2 - M3SD2 + dMLb + mz**2*np.cos(2*beta)*(-1/2.0+2/3.0*sthetaw**2))**2 + 4.0*mb2*(Ab-muSUSY*tanb)**2 )/2.0)
        msbot2=np.sqrt(M3SQ2/2.0 + M3SD2/2.0 + dMLb/2.0 + mb2 - mz**2/4.0*np.cos(2*beta)+np.sqrt((M3SQ2 - M3SD2 + dMLb + mz**2*np.cos(2*beta)*(-1/2.0+2/3.0*sthetaw**2))**2 + 4.0*mb2*(Ab-muSUSY*tanb)**2 )/2.0)	
        msb12 = msbot1**2
        msb22 = msbot2**2   
        
        thetab = np.arcsin(2.0*mb/(msb12-msb22)*(Ab-muSUSY*tanb))/2.0
        
        if(2.0*mb/(msb12-msb22)*(Ab-muSUSY*tanb)>1.0 and 2.0*mb/(msb12-msb22)*(Ab-muSUSY*tanb)<1.00010):
            thetab = np.arcsin(1.0)/2.0 
        
        if(2.0*mb/(msb12-msb22)*(Ab-muSUSY*tanb)<-1.0 and 2.0*mb/(msb12-msb22)*(Ab-muSUSY*tanb)>-1.00010):
            thetab = np.arcsin(-1.0)/2.0    
        
        sthetab = np.sin(thetab)
        cthetab = np.cos(thetab)    
        
        if(np.abs(cthetab**2 * msb12 + sthetab**2 * msb22 - MsbotEig) > 10**-3):
            thetab = Pi/2.0 - thetab
            sthetab = np.sin(thetab)
            cthetab = np.cos(thetab)    
        
        t2thetab = np.tan(2*thetab)
        c2thetab = cthetab**2 - sthetab**2
        s2thetab = 2*sthetab*cthetab    
    
        return {'mst12':mst12, 'mst22':mst22, 
            'sthetat':sthetat, 'cthetat':cthetat,
            's2thetat':s2thetat, 'c2thetat':c2thetat, 
            't2thetat':t2thetat, 'msb12':msb12, 'msb22':msb22, 
            'sthetab':sthetab, 'cthetab':cthetab, 's2thetab':s2thetab, 
            'c2thetab':c2thetab, 't2thetab':t2thetab}

    # mstop1, mstop2, msb1, msb2 should be on some benchmark scenario.
    def deltaML(self,  
        mstop1:float=1.0, 
        mstop2:float=1.5, 
        cthetat:float=0.0, 
        sthetat:float=0.0,
        msb1:float=1.0, 
        msb2:float=1.5, 
        cthetab:float=0.0, 
        sthetab:float=0.0,
        mu: float=einital().initial()['muRggh'],
        mtop: float=inputs().mt,
        mb: float= einital().initial()['mbos'],
        mgluino: float=inputs().M3, 
        )->float:
        tanb = self.tanb

        c2thetat = cthetat**2 - sthetat**2
        s2thetat = 2*sthetat*cthetat
        c2thetab = cthetab**2 - sthetab**2
        s2thetab = 2*sthetab*cthetab

        lb = 2.0*np.log(mu/mb)
        lb1 = 2.0*np.log(mu/msb1)
        lb2 = 2.0*np.log(mu/msb2)
        lt = 2.0*np.log(mu/mtop)
        lt1 = 2.0*np.log(mu/mstop1)
        lt2 = 2.0*np.log(mu/mstop2)
        lg = 2.0*np.log(mu/mgluino)

        deltaML =(-8*mb**2 - 4*lb*mb**2 + 4*msb1**2 + 3*c2thetab*msb1**2 +
            2*lb1*msb1**2 + c2thetab*lb1*msb1**2 + 4*msb2**2 -
            3*c2thetab*msb2**2 + 2*lb2*msb2**2 - c2thetab*lb2*msb2**2 -
            4*mstop1**2 - 3*c2thetat*mstop1**2 - 2*lt1*mstop1**2 -
            c2thetat*lt1*mstop1**2 - 4*mstop2**2 + 3*c2thetat*mstop2**2 -
            2*lt2*mstop2**2 + c2thetat*lt2*mstop2**2 + 8*mtop**2 +
            4*lt*mtop**2 + (mb**2 + mgluino**2 - msb1**2 -
            2*mb*mgluino*s2thetab)*b0fin(mb**2,mgluino,msb1,mu) +
            (mb**2 + mgluino**2 - msb2**2 + 2*mb*mgluino*s2thetab)*
            b0fin(mb**2,mgluino,msb2,mu) +
            mb**2*b0fin(msb1**2,mb,mgluino,mu) +
            c2thetab*mb**2*b0fin(msb1**2,mb,mgluino,mu) +
            mgluino**2*b0fin(msb1**2,mb,mgluino,mu) +
            c2thetab*mgluino**2*b0fin(msb1**2,mb,mgluino,mu) -
            msb1**2*b0fin(msb1**2,mb,mgluino,mu) -
            c2thetab*msb1**2*b0fin(msb1**2,mb,mgluino,mu) -
            2*mb*mgluino*s2thetab*b0fin(msb1**2,mb,mgluino,mu) +
            mb**2*b0fin(msb2**2,mb,mgluino,mu) -
            c2thetab*mb**2*b0fin(msb2**2,mb,mgluino,mu) +
            mgluino**2*b0fin(msb2**2,mb,mgluino,mu) -
            c2thetab*mgluino**2*b0fin(msb2**2,mb,mgluino,mu) -
            msb2**2*b0fin(msb2**2,mb,mgluino,mu) +
            c2thetab*msb2**2*b0fin(msb2**2,mb,mgluino,mu) +
            2*mb*mgluino*s2thetab*b0fin(msb2**2,mb,mgluino,mu) -
            mgluino**2*b0fin(mstop1**2,mtop,mgluino,mu) -
            c2thetat*mgluino**2*b0fin(mstop1**2,mtop,mgluino,mu) +
            mstop1**2*b0fin(mstop1**2,mtop,mgluino,mu) +
            c2thetat*mstop1**2*b0fin(mstop1**2,mtop,mgluino,mu) -
            mtop**2*b0fin(mstop1**2,mtop,mgluino,mu) -
            c2thetat*mtop**2*b0fin(mstop1**2,mtop,mgluino,mu) +
            2*mgluino*mtop*s2thetat*b0fin(mstop1**2,mtop,mgluino,mu) -
            mgluino**2*b0fin(mstop2**2,mtop,mgluino,mu) +
            c2thetat*mgluino**2*b0fin(mstop2**2,mtop,mgluino,mu) +
            mstop2**2*b0fin(mstop2**2,mtop,mgluino,mu) -
            c2thetat*mstop2**2*b0fin(mstop2**2,mtop,mgluino,mu) -
            mtop**2*b0fin(mstop2**2,mtop,mgluino,mu) +
            c2thetat*mtop**2*b0fin(mstop2**2,mtop,mgluino,mu) -
            2*mgluino*mtop*s2thetat*b0fin(mstop2**2,mtop,mgluino,mu) -
            mgluino**2*b0fin(mtop**2,mgluino,mstop1,mu) +
            mstop1**2*b0fin(mtop**2,mgluino,mstop1,mu) -
            mtop**2*b0fin(mtop**2,mgluino,mstop1,mu) +
            2*mgluino*mtop*s2thetat*b0fin(mtop**2,mgluino,mstop1,mu) -
            mgluino**2*b0fin(mtop**2,mgluino,mstop2,mu) +
            mstop2**2*b0fin(mtop**2,mgluino,mstop2,mu) -
            mtop**2*b0fin(mtop**2,mgluino,mstop2,mu) -
            2*mgluino*mtop*s2thetat*b0fin(mtop**2,mgluino,mstop2,mu)
            )/3.0
        return deltaML

    def get_json(self)->Dict:
        return json.dumps(self.calcsquarkmass(), cls=NumpyEncoder, indent=4)
