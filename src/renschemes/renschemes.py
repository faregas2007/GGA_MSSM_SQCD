# renschemes.py
# renormalize schemes from sushi
# Micheal ???

import sys
import logging
from typing import Dict, List

from interface import *
from init import *
from consts import *
from config import *
from utils import *
from evalsusy import *
from init import *


from renschemes.Abrenorm import *
from renschemes.mbrenorm import *
from renschemes.sbottomrenorm import *
from renschemes.squarkmatrix import *

# variation of facade pattern desgin

class renschemes:
    def __init__(self,
        inputs, 
        einital: inputs, 
        renormalizeSM: renorm, 
        quarkhiggscoup: coupling, 
        squarkhiggscoupCMSSM: coupling, 
        renormalizetest: renorm,
        #renormalize: renorm
        ):
        self._inputs = inputs()
        self._einital = einital().initial()
        self._renormalizeSM = renormalizeSM()
        self._quarkhiggscoup = quarkhiggscoup()
        self._squarkhiggscoupCMSSM = squarkhiggscoupCMSSM()
        self._renormalize = renormalize()
        self._renormalizetest = renormalizetest()

    def get_params(self)->Dict:
        logger.info('renscheme starts ...')
        temp0 = self._einital
        temp = self._renormalizeSM.get_params()
        
        logger.info(f'eintial: {temp0}')
        logger.info(f'renormalizeSM: {temp} ...')
        temp.update(temp0)
        
        # hard copy into logs, printing each stages of computations.
        #temp0 = self._renormalize.get_params()
        temp0 = self._renormalizetest.get_params()
        temp.update(temp0)
        logger.info(f'renormalize params: {temp0}')

        temp0 = self._quarkhiggscoup.get_coup()
        logger.info(f'quarkhiggscoup: {temp0}')
        temp.update(temp0)

        thetab = np.arcsin(temp['sthetab'])
        thetat = np.arcsin(temp['sthetat'])
        temp0 = self._squarkhiggscoupCMSSM.get_coup(mb=temp['mbsb'], thetat=thetat, thetab=thetab)
        logger.info(f'squarkhiggscoup: {temp0}')
        temp.update(temp0)
        
        if(self._inputs.yukfac[6] == 0.0):
            temp['gb1'] = complex(0.0,0.0)
            temp['gb2'] = complex(0.0,0.0)
        if(self._inputs.yukfac[5] == 0.0):
            temp['gt1'] = complex(0.0,0.0)
            temp['gt2'] = complex(0.0,0.0)
        
        temp['gb'] = temp['gbh'] * temp['gb']
        temp['gc'] = temp['gth'] * self._inputs.yukfac[0]
        temp['gt'] = temp['gth'] * self._inputs.yukfac[1]
        temp['gb1'] = temp['gbh11'] * self._inputs.yukfac[6]
        temp['gb2'] = temp['gbh22'] * self._inputs.yukfac[6]
        temp['gt1'] = temp['gth11'] * self._inputs.yukfac[5]
        temp['gt2'] = temp['gth22'] * self._inputs.yukfac[5]
        
        logger.info(f'final: {temp}')
        logger.info('renscheme ends...')
        return temp

    def get_json(self, obj:str)->Dict:
        return json.dumps(self.get_params(), cls=NumpyEncoder, indent=4)

    def to_file(self):
        return save_dict(self.get_params(), Path(model_dir, 'models.json'), cls=NumpyEncoder)

class renormalizeSM(renorm, einital):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cinit = super().initial() 
        self.mbos = self.cinit['mbos']      
        self.alphasggh = self.cinit['alphasggh']
        
        self.mbMSbarmuR = self.cinit['mbMSbarmuR']       
        self.mbrunloop = self.cinit['mbrunloop']
    
    def get_params(self)->Dict:
        gb1 = np.array([0.0,0.0])
        gb2 = np.array([0.0,0.0])
        gt1 = np.array([0.0,0.0])
        gt2 = np.array([0.0,0.0])
        if(self.mbrunloop == 0 and self.mbrunyuk == 0):
            mb = self.mbos
            mbyuk = self.mbmb
            dmb = 0.0
            gb = 1.0
            dgb = 0.0
        elif(self.mbrunloop == 0 and self.mbrunyuk == 1):
            mb = self.mbos
            mbyuk = self.mbmb
            dmb = 0.0
            gb = self.mbmb/self.mbos * yukfac[2]
            dgb = cf * self.alphasggh/Pi * Pi / self.alphasggh # 1-loop
        elif(self.mbrunloop == 0 and self.mbrunyuk == 2):
            mb = self.mbos
            mbyuk = self.mbMSbarmuR
            dmb = 0.0
            gb = self.mbMSbarmuR/self.mbos * yukfac[2]
            dgb = (-cf + 2.0*np.log(self.mbmb/self.muRggh))*self.alphasggh/Pi * Pi/self.alphasggh # 1-loop
        else:
            print("Renormalization scheme for bottom")
            print("Sector not known in the SM/H2DM")
        
        return {'mb':float(mb), 'mbyuk':float(mbyuk), 'gb':float(gb), 'dmb':float(dmb), 'dgb':float(dgb), 'gb1':gb1, 'gb2':gb2, 'gt1':gt1, 'gt2':gt2}

    def get_json(self):
        return json.dumps(self.get_params(), cls=NumpyEncoder, indent=4)

"""
test = renormalizeSM()
print(test.get_json())
"""

class renormalize(renorm, einital):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cinit = super().initial() 
        self.mbos = self.cinit['mbos']
    
        self.alphasggh = self.cinit['alphasggh']
        self.mbMSbarmuR = self.cinit['mbMSbarmuR']
        self.muRggh = self.cinit['muRggh']
        self.mbrunloop = self.cinit['mbrunloop']

    def get_params(self)->Dict:
        mbMSbarmuR = self.mbMSbarmuR
        alphasggh = self.alphasggh
        muRggh = self.muRggh

        mgl = self.M3
        alphasmz = self.alphasmz
        mZ = self.mZ
        mbos = self.mbos
        mt = self.mt

        mbsbch = self.mbsbch
        Abch = self.Abch
        thetach = self.thetach
        tanbresum = self.tanbresum
        beta = np.arctan(self.tanb)

        Ab = self.Ab
        At = self.At

        tanb = self.tanb
        muSUSY = self.mu

        sqmass = squarkmatrix().calcsquarkmass()
        msb12 = sqmass['msb12'] 
        msb22 = sqmass['msb22']
        mst12 = sqmass['mst12']
        mst22 = sqmass['mst22']

        cthetat = sqmass['cthetat']
        sthetat = sqmass['sthetat']
        cthetab = sqmass['cthetab']
        sthetab = sqmass['sthetab']

        c2thetat = sqmass['c2thetat']
        s2thetat = sqmass['s2thetat']
        c2thetab = sqmass['c2thetab']
        s2thetab = sqmass['s2thetab']

        # reconstruct mbsb
        mbsb = s2thetab/((Ab -muSUSY*tanb))*(msb12-msb22)/2.0
        if(mbsb < 0.0):
            print('sbottom sector input is inconsistent!')
        
        muD = (mgl+np.sqrt(msb12)+np.sqrt(msb22))/3.0
        mbMSbarmuD = mbrenorm().mbMSbarmuX(mu_in=muD)
        mbDRSMmuD = mbrenorm().mbMSDRtrans(mb=mbMSbarmuD, mu_in=muD)
        
        alphasren = SUSHI_alphas(alphasmz, rmurl=muD, mZ=mZ)
        temp = mbrenorm().deltab(
            msb12=msb12, 
            msb22=msb22, 
            mu_in=muD, 
            alphasren=alphasren)
        delmb = temp['delmb']
        delAb = temp['delAb']

        # calculate OS-sbottom masses, 
        dMLb = (alphasren/Pi)*squarkmatrix().deltaML(
            mstop1 = np.sqrt(mst12),
            mstop2 = np.sqrt(mst22),
            cthetat = cthetat,
            sthetat = sthetat,
            mb = mbos,
            msb1 = np.sqrt(msb12),
            msb2 = np.sqrt(msb22),
            cthetab = cthetab,
            sthetab = sthetab,
            mu=muD)

        # mb is chosen to be mbos--> should put as a class input
        mb = mbos
        mb2 = mb*mb
        mt2 = mt*mt
        mgl2 = mgl*mgl
        muD2 = muD*muD

        temp = sbottomrenorm().sbottomrenormHM(
            mbsb = mbsb,
            msb12 = msb12,
            msb22 = msb22,
            sthetab = sthetab,
            cthetab = cthetab,
            mu_in = muD,
            mb=mb
        )

        mbsb = sbottomrenorm().mb_sbottom(
            mb = mb,
            mbsb_var1=mbsb,
            mbsb_var2=mbDRSMmuD,
            msb12=msb12,
            msb22=msb22,
            dmsb11=temp['dmsb11'],
            dmsb21=temp['dmsb21'],
            s2thetab=s2thetab,
            c2thetab=c2thetab,
            dthetabos1=temp['dthetab1'],
            muD=muD,
            prefac=alphasren/Pi
        )

        # recalculating squark masses and mixing angles with shifted dMLb
        sqmass = squarkmatrix().calcsquarkmass(mb=mbsb, dMLb=dMLb)
        msb12 = sqmass['msb12'] 
        msb22 = sqmass['msb22']
        mst12 = sqmass['mst12']
        mst22 = sqmass['mst22']

        cthetat = sqmass['cthetat']
        sthetat = sqmass['sthetat']
        cthetab = sqmass['cthetab']
        sthetab = sqmass['sthetab']

        
        c2thetat = sqmass['c2thetat']
        s2thetat = sqmass['s2thetat']
        c2thetab = sqmass['c2thetab']
        s2thetab = sqmass['s2thetab']

        t2thetab = sqmass['t2thetab']

        thetat = np.arcsin(sthetat)
        thetab = np.arcsin(sthetab)

        mb = mbsb
        mb2 = mb*mb

        # caculate counter-terms dAb, dthetab, dmbsb
        temp = sbottomrenorm().sbottomrenormHM(
            mbsb=mbsb,
            msb12=msb12,
            msb22=msb22,
            sthetab=sthetab,
            cthetab=cthetab,
            mu_in=muD,
            mb=mb
        )

        dthetabos = np.array([temp['dthetab1'], temp['dthetab2']])
        dmsb1 = np.array([temp['dmsb11'], temp['dmsb21']])
        dmsb2 = np.array([temp['dmsb21'], temp['dmsb22']])
        dmbsbos = np.array([temp['dhbhb1'], temp['dhbhb2'], temp['dhbhb3']])

        # 2-0-0 is a default renshemes
        # various parameters choice
        dep = 3
        dAb = np.array([0.0,0.0,0.0])
        dthetab = np.array([0.0,0.0,0.0])
        dmbsb = np.array([0.0,0.0,0.0])
        if(mbsbch == 0):
            dmbsb[0] = 0.0
            dmbsb[1] = dmbsbos[0]
            dmbsb[2] = dmbsbos[1]
        elif(mbsbch == 1):
            dmbsb[0] = 0.0
            dmbsb[1] = 0.0
            dmbsb[2] = 0.0
        else:
            dep = 0

        if(thetach == 0):
            # thetab OS
            dthetab[0] = 0.0
            dthetab[1] = dthetabos[0]
            dthetab[2] = dthetabos[1]
        elif(thetach == 1): 
            dthetab[0] = 0.0
            dthetab[1] = 0.0
            dthetab[2] = 0.0
        else:
            if(dep == 0):
                print("error:two dependent quantities.")
                exit
            dep = 1

        if(Abch == 0):
            if(dep == 0):
                dAb[2] = Abrenorm().dAbosfintbHM(
                    msb12=msb12, 
                    msb22=msb22, 
                    sthetab=sthetab, 
                    cthetab=cthetab, 
                    dmsb1_1=dmsb1[0], 
                    dmsb2_1=dmsb2[0], 
                    deltab=dthetab[1], 
                    mu_in=muRggh, 
                    beta=beta)
            elif(dep == 1): 
                
                dAb[2] = Abrenorm().dAbosfinmbHm(
                    msb12=msb12, 
                    msb22=msb22, 
                    deltamb=dmbsb[1], 
                    mu_in=muRggh,
                    beta = beta, 
                    mb = mb)
            else:
                print("error:no depdendent quantity.")
                exit
        elif(Abch == 1):
            dAb[2] = 0.0
        else:
            if(dep < 2):
                print('error: two dependent quantities.')
            dep = 2
        
        if(dep == 0):
            # mb dependent
            dmbsb[0] = 0.0
            dmbsb[1] = 2.0*mb/(msb12 - msb22)*(dmsb1[0]*msb12 - dmsb2[0]*msb22)
            if(Abch == 0):
                dmbsb[1] = dmbsb[1] + 2.0/t2thetab*mb*dthetab[1] - 2.0*mb2/(msb12-msb22)/s2thetab*dAb[2]
            elif(np.s2thetab > 10**-8):
                dmbsb[1] = dmbsb[1] - (c2thetab*(msb12-msb22)*np.cos(beta)*np.sin(beta))/muSUSY*dthetab[1]- 2.0*mb2/(msb12-msb22)/s2thetab*Abrenorm().dAbosfintbHM(msb12=msb12, 
                        msb22=msb22, 
                        sthetab=sthetab, 
                        cthetab=cthetab, 
                        dmsb1_1=dmsb1[0], 
                        dmsb2_1=dmsb2[0], 
                        deltab=0.0, 
                        mu_in=muRggh)
            else:
                dmbsb[1] = sbottomrenorm().mb_sbottom(
                    mb = mb,
                    mbsb_var1=dmbsb[1], 
                    mbsb_var2=dmbsb[1], 
                    msb12=msb12, 
                    msb22=msb22, 
                    dmsb11=dmsb1[0], 
                    dmsb22=dmsb2[0], 
                    s2thetab=s2thetab, 
                    c2thetab=c2thetab, 
                    dthetabos1=dthetab[1], 
                    muD=muD, 
                    prefac=1)
            dmbsb[2] = 0.0
        elif(dep == 1):
            print("error: thetab cannot be chosen as a dependent quantity.")
            exit
            #dthetab[0] = t2thetab/mb/2.0 * dmbsb[0]
            #dthetab[1] = -t2thetab/(msb12-msb22)*(dmsb1[0]*msb12 - dmsb[0]*msb22) + t2thetab/mb/2.0 * dmbsb[1]
            #dthetab[2] = dthetab[1] + mb/(msb12 - msb22)/c2thetab*dAb[2]
        elif(dep == 2):
            dAb[0] = -s2thetab*(msb12-msb22)/2.0/mb2*dmbsb[0]
            dAb[1] = -s2thetab/mb* (dmsb1[0]*msb12 - dmsb2[0]*msb22) - s2thetab*(msb12-msb22)/2.0/mb2 * dmbsb[1] + c2thetab*(msb12- msb22)/mb * dthetab[1]
            dAb[2] =  s2thetab/mb * (dmsb1[1]*msb12 - dmsb2[1]*msb22) - s2thetab*(msb12-msb22)/2.0/mb2* dmbsb[2] + c2thetab*(msb12-msb22)/mb * dthetab[2]
        else:
            print('error: no dependent quantity.')
            exit

        dthetabstd = dthetabos[0]
        dmbstd = 2.0*mb/(msb12-msb22)*(dmsb1[0]*msb12 - dmsb2[0]*msb22)
        dAbstd = Abrenorm().dAbosfintbHM(
                    msb12=msb12, 
                    msb22=msb22, 
                    sthetab=sthetab, 
                    cthetab=cthetab, 
                    dmsb1_1=dmsb1[0], 
                    dmsb2_1=dmsb2[0], 
                    deltab=dthetabstd, 
                    mu_in=muRggh)

        dmbstd =sbottomrenorm().mb_sbottom(
                    mb = mb,
                    mbsb_var1=dmbstd, 
                    mbsb_var2=dmbstd, 
                    msb12=msb12, 
                    msb22=msb22, 
                    dmsb11=dmsb1[0], 
                    dmsb21=dmsb2[0], 
                    s2thetab=s2thetab, 
                    c2thetab=c2thetab, 
                    dthetabos1=dthetabstd, 
                    muD=muD, 
                    prefac=1)

        Abstd = 0.0
        Abstd = Ab + alphasggh/Pi * (dAbstd - (dAb[0] + dAb[1] + dAb[2]))
        thetab = thetab + alphasggh/Pi * (dthetabstd - (dthetab[0] + dthetab[1]))
        mbsb = mbsb + alphasggh/Pi * (dmbstd - (dmbsb[0] + dmbsb[1]))
        dmbsb = dmbsb/mbsb

        # setting for the internal yukawa couplings
        # relevant for QCD calculation
        # mbrunloop will always be set to 0, i.e mb is chosen to be on-shell.
        mbyuk = mbos
        mb = mbos
        dgb = 0.0
        dmb = np.array([0.0,0.0,0.0])

        if(tanbresum == 1):
            gb = self.yukfac[2]/(1.0 + delmb)
        elif(tanbresum == 2):
            gb = self.yukfac[2]*(1.0 - (1.0/tanb**2)*delmb - Amx[Sind, 2]*vev/(Amx[Sind,1]*vevS*tanb)*delmb)*(1.0/(1.0+delmb))
        else:
            gb = self.yukfac[2]  
        
        cthetab = np.cos(thetab)
        sthetab = np.sin(thetab)
        #t2thetab = np.tan(2*thetab) 
        #c2thetab = cthetab**2 - sthetab**2
        #s2thetab = 2*sthetab*cthetab

        return {'mb': mb, 'mbyuk': mbyuk, 'mbsb': mbsb, 
            'sthetab': sthetab, 'cthetab': cthetab, 'Ab': Abstd, 
            'msb12':msb12, 'msb22':msb22, 
            'gb': gb, 'dgb': dgb, 'dmb': dmb, 
            'mt':mt, 'sthetat': sthetat, 'cthetat':cthetat, 'At':At, 
            'mst12':mst12, 'mst22':mst22}
    
    def get_json(self):
        return json.dumps(self.get_params(), cls=NumpyEncoder, indent=4)

"""
test = renormalize()
print(test.get_json())
"""

class renormalizetest(renorm, einital):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cinit = super().initial() 
        self.mbos = self.cinit['mbos']
    
        self.alphasggh = self.cinit['alphasggh'] 
        self.mbmsbarmur = self.cinit['mbMSbarmuR'] 
        self.mbrunloop = self.cinit['mbrunloop']

    def get_params(self, muD: float=1500.0)->dict:
        M3SU = inputs().MQ3
        M3SQ = inputs().MTR
        M3SD = inputs().MBR

        tanbresum = self.tanbresum

        Ab = self.Ab
        At = self.At
        # test
        Abtree = Ab
        Attree = At

        alphasmz = self.alphasmz         
        tanb = self.tanb
        muSUSY = self.mu

        mt = self.mt
        mb = self.mbos
        mbmb = self.mbmb
        mbos = self.mbos

        mgl = self.M3
        mZ = self.mZ
        mW = self.mW

        beta = np.arctan(self.tanb)

        # mb-dep extension
        # sushi negative cthetab, cthetat --> mbsb increases with tanb
        # Micheal positive cthetab, cthetat --> mbsb decreases with tanb. 
        sqmass = squarkmatrix().calcsquarkmass(mb=mbmb, dMLb=0.0)
        msb12 = sqmass['msb12'] 
        msb22 = sqmass['msb22']
        mst12 = sqmass['mst12']
        mst22 = sqmass['mst22']

        cthetat = sqmass['cthetat']
        sthetat = sqmass['sthetat']
        cthetab = sqmass['cthetab']
        sthetab = sqmass['sthetab']

        c2thetat = sqmass['c2thetat']
        s2thetat = sqmass['s2thetat']
        c2thetab = sqmass['c2thetab']
        s2thetab = sqmass['s2thetab']
        print(f's2thetab:{s2thetab}')        
        mbsb = s2thetab/((Ab - muSUSY*tanb))*(msb12 -msb22)/2.0
        mbMSbarmuD = mbrenorm().mbMSbarmuX(mu_in=muD)
        alphasren = SUSHI_alphas(alphasmz, rmurl=muD, mZ=mZ)
        
        # mb-> mbmb.
        #mb = mbmb
        mb = mbmb
        mb2 = mb*mb
        mt2 = mt*mt 
        mgl2 = mgl*mgl
        muD2 = muD*muD

        # MSbar (smflag=True, drflag=False)
        # DRbar (smflag=True, drflag=True)
        counterterms = sbottomrenorm().sbottomrenormHM(
            mbsb = mbsb,
            msb12 = msb12,
            msb22 = msb22,
            sthetab = sthetab,
            cthetab = cthetab,
            mu_in = muD,
            mgl = mgl,
            mb = mbos, # mbos
            smflag=True,
            drflag=False
        )
    
        dmsb1os = counterterms['dmsb11']
        dmsb2os = counterterms['dmsb21']
        dthetabstd = counterterms['dthetab1']
        dmb = counterterms['dhbhb2']        

        dAb = Abrenorm().dAbos(
            msb12=msb12,
            msb22=msb22,
            sthetab=sthetab,
            cthetab=cthetab,
            dmsb1_1=dmsb1os,
            dmsb2_1=dmsb2os,
            dthetab=dthetabstd,
            mu_in=muD,
            beta=beta,
            mb=mbos,
            mgl=mgl,
            Ab = Abtree,
            muSUSY=muSUSY
        ) 
        """
        dmbsb = sbottomrenorm().dmbdep(
            Ab = Ab, 
            mu = muSUSY,
            msb12 = msb12,
            msb22 = msb22,
            sthetab = sthetab,
            cthetab = cthetab,
            dthetab = dthetabstd, 
            dmsb1_1 = dmsb1os,
            dmsb2_1= dmsb2os
        ) 
        """
         
        dmbsb = sbottomrenorm().dmbdep2(
            Ab = Ab, 
            mu = muSUSY,
            tanb = tanb,
            mb = mbos, 
            msb12 = msb12,
            msb22 = msb22,
            sthetab = sthetab,
            cthetab = cthetab,
            dthetab = dthetabstd, 
            dmsb12 = dmsb1os,
            dmsb22= dmsb2os,
            dAb = dAb
        )        
       
        thetab = np.arcsin(sthetab)
        prefac = SUSHI_alphas(alphasmz, muD, mZ)/Pi
        mbsb = mbMSbarmuD + prefac*(-dmbsb + dmb)
        #Ab = Ab + prefac*dAb 
        print(f"alphas:{prefac*Pi}")
        print(f'mbsb:{mbsb}')
        #print(f'Ab, dAb: {Ab}, {dAb}')
        print(f'thetab:{thetab}')

        dMLb = (alphasren/Pi)*squarkmatrix().deltaML(
            mstop1 = np.sqrt(mst12) ,
            mstop2 = np.sqrt(mst22),
            cthetat = cthetat,
            sthetat = sthetat,
            msb1 = np.sqrt(msb12),
            msb2 = np.sqrt(msb22),
            cthetab = cthetab,
            sthetab = sthetab,
            mu = muD,
            mtop = mt,
            mb = mbos, # mbmb, mbos
            mgluino = mgl)
        
        sqmass = squarkmatrix().calcsquarkmass(
            M3SU = M3SU,
            M3SQ = M3SQ,
            M3SD = M3SD,
            mw = mW,
            mz = mZ,
            mt = mt,
            mb = mbsb, # mbos
            muSUSY = muSUSY,
            beta = beta,
            At = At, 
            Ab = Abtree, 
            dMLb = dMLb
        )

        msb12 = sqmass['msb12'] 
        msb22 = sqmass['msb22']
        mst12 = sqmass['mst12']
        mst22 = sqmass['mst22']

        cthetat = sqmass['cthetat']
        sthetat = sqmass['sthetat']
        cthetab = sqmass['cthetab']
        sthetab = sqmass['sthetab']

        c2thetat = sqmass['c2thetat']
        s2thetat = sqmass['s2thetat']
        c2thetab = sqmass['c2thetab']
        s2thetab = sqmass['s2thetab']
        print(f's2thetab_2:{s2thetab}')
        thetabos = np.arcsin(sthetab)
        thetatos = np.arcsin(sthetat)

        # recalculate counter-terms with OS-values. 
        counterterms = sbottomrenorm().sbottomrenormHM(
            mbsb = mbsb,
            msb12 = msb12,
            msb22 = msb22,
            sthetab = sthetab,
            cthetab = cthetab,
            mu_in = muD,
            mgl = mgl,
            mb = mbos,
            smflag = True,
            drflag = False
        )  

        dmsb11 = counterterms['dmsb11']
        dmsb21 = counterterms['dmsb21']
        dthetabos = counterterms['dthetab1']
        dmb = counterterms['dhbhb2'] 
        
        dAb = Abrenorm().dAbos(
            msb12=msb12,
            msb22=msb22,
            sthetab=sthetab,
            cthetab=cthetab,
            dmsb1_1=dmsb1os,
            dmsb2_1=dmsb2os,
            dthetab=dthetabos,
            mu_in=muD,
            beta=beta,
            mb=mbos,
            mgl=mgl,
            Ab=Abtree,
            muSUSY=muSUSY
        )

        print(f'Ab, dAb:{Ab}, {dAb}')
        
        """
        dmbsb = sbottomrenorm().dmbdep(
            Ab = Ab, 
            mu = muSUSY,
            msb12 = msb12,
            msb22 = msb22,
            sthetab = sthetab,
            cthetab = cthetab,
            dthetab = dthetabos, 
            dmsb1_1 = dmsb1os,
            dmsb2_1= dmsb2os
        )
        """
                
        dmbsb = sbottomrenorm().dmbdep2(
            Ab = Ab, 
            mu = muSUSY,
            tanb = tanb,
            mb = mbos, 
            msb12 = msb12,
            msb22 = msb22,
            sthetab = sthetab,
            cthetab = cthetab,
            dthetab = dthetabos, 
            dmsb12 = dmsb1os,
            dmsb22= dmsb2os,
            dAb=dAb
        )        
        
        # mbsb != 1.8
        # fixed Ab
        thetab = np.arcsin(sthetab)
        mbsb = mbMSbarmuD + prefac*(-dmbsb + dmb)
        Abstd = Abtree
        #Abstd = Ab + prefac*dAb
        #thetabos = thetabos + prefac*dthetabos

        print(f'mbsb:{mbsb}')
        print(f'thetab:{thetab}')
        print(f'Ab:{Abstd}')

        # yukawa bottom quark
        mbyuk = mbos
        mb = mbos
        dgb = 0.0
        dmb = np.array([0.0,0.0,0.0])

        if(tanbresum == 1):
            gb = self.yukfac[2]/(1.0 + delmb)
        elif(tanbresum == 2):
            gb = self.yukfac[2]*(1.0 - (1.0/tanb**2)*delmb - Amx[Sind, 2]*vev/(Amx[Sind,1]*vevS*tanb)*delmb)*(1.0/(1.0+delmb))
        else:
            gb = self.yukfac[2]  
 
        return {
            'mt':mt,
            'mb': mb,
            'mbsb': mbsb,
            'mbos':mbos,
            'mbMSbarmuD':mbMSbarmuD,

            'msb12':msb12,
            'msb22':msb22,
            'mst12':mst12,
            'mst22':mst22,
        
            'Abtree':Abtree,
            'Ab':Abstd,
            'Attree':Attree,
            'At':At,

            'thetabos':thetabos,
            'thetatos':thetatos,
            'sthetab':sthetab,
            'cthetab':cthetab,
            'sthetat':sthetat,
            'cthetat':cthetat,
           
            's2thetab':s2thetab,
            'c2thetab':c2thetab,
            's2thetat':s2thetat,
            'c2thetat':c2thetat,

            'tanbeta':tanb,
            'mbyuk':mbyuk,
            'dgb':dgb,
            'dmb':dmb
                }

    def get_json(self):
        return json.dumps(self.get_params(), cls=NumpyEncoder, indent=4)
 
