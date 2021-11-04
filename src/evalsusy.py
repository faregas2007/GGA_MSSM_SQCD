import numpy as np
import json

from include import *
from init import *

# factory pattern --> will move into include.py/template.py 
class quarkhiggscoup(coupling, einital):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cinit = super().initial()
        self.Sind = self.cinit['Sind']
        self.beta = np.arctan(self.tanb)
    
    def get_coup(self):
        tanb = self.tanb
        sbeta = np.sin(self.beta)
        cbeta = np.cos(self.beta)
        Amx = self.Amx
        Sind = self.Sind

        # SM
        if(self.model == 0):
            gth = 1.0
            gbh = 1.0
        
        # MSSM
        if(self.model == 1):
            gth = Amx[Sind, 1]*1.0/tanb
            gbh = Amx[Sind, 1]*tanb
        
        return {'gth':gth, 'gbh':gbh}

    """    
    def update(self, **kwargs)->Dict:
        return inputs().update(kwargs)
    """

    def get_json(self):
        return json.dumps(self.get_coup(), indent=4)


class squarkhiggscoupCMSSM(coupling, einital):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.cinit = super().initial()
        self.Sind = self.cinit['Sind']
        self.lpseudo = self.cinit['Sind']

    def get_coup(self, 
        mb: float=einital().initial()['mbos'], 
        thetat: float=0.0,
        thetab: float=0.0):
        Sqrt2 = np.sqrt(2)
        vev = 1.0/np.sqrt(Sqrt2*self.GF)
        beta = np.arctan(self.tanb)

        mt = self.mt
        mZ = self.mZ

        mu = self.mu
        At = self.At
        Ab = self.Ab

        SiB = np.sin(beta)
        CoB = np.cos(beta)
        TaB = np.tan(beta)

        SiB2 = np.sin(beta)*np.sin(beta)
        CoB2 = np.cos(beta)*np.cos(beta)
        TaB2 = np.tan(beta)*np.tan(beta)

        ctW = self.mW/self.mZ
        stW = np.sqrt(1.0 - ctW**2)
        SiW2 = stW*ctW
        CoW2 = ctW*ctW
        c2W = np.cos(2.0*np.arccos(ctW))
        s2W = np.sin(2.0*np.arccos(ctW))

        CoI = complex(0.0, 1.0)
        Cmu = np.conj(mu)
        CAt = np.conj(At)
        CAb = np.conj(Ab)

        cthetat = np.cos(thetat)
        sthetat = np.sin(thetat)
        cthetab = np.cos(thetab)
        sthetab = np.sin(thetab)

        RU11 = cthetat
        RU12 = sthetat
        RU21 = -sthetat
        RU22 = cthetat
        RD11 = cthetab
        RD12 = sthetab
        RD21 = -sthetab
        RD22 = cthetab

        CRU11 = np.conj(RU11)
        CRU12 = np.conj(RU12)
        CRU21 = np.conj(RU21)
        CRU22 = np.conj(RU22)
        CRD11 = np.conj(RD11)
        CRD12 = np.conj(RD12)
        CRD21 = np.conj(RD21)
        CRD22 = np.conj(RD22)        

        mb2 = mb*mb
        mt2 = mt*mt
        mZ2 = mZ*mZ

        Ht11Hmx = np.zeros(4, dtype=complex)
        Ht12Hmx = np.zeros(4, dtype=complex)
        Ht21Hmx = np.zeros(4, dtype=complex)
        Ht22Hmx = np.zeros(4, dtype=complex)

        Hb11Hmx = np.zeros(4, dtype=complex)
        Hb12Hmx = np.zeros(4, dtype=complex)
        Hb21Hmx = np.zeros(4, dtype=complex)
        Hb22Hmx = np.zeros(4, dtype=complex)

        Ht11Hmx[0] = (1/3.0*(-1.0*CRU12*(3.0*mt*mu*RU11 - 4.*CoB*mZ2*RU12*SiB*SiW2) - 1.0*CRU11*(3.0*Cmu*mt*RU12 + CoB*mZ2*RU11*SiB*(-3.0*CoW2 + SiW2))))/SiB
        Ht11Hmx[1] = (1/3.0*(3.0*CAt*CRU12*mt*RU11 + 2.0*CRU12*RU12*(3.0*mt2 - 2.0*mZ2*SiB2*SiW2) + CRU11*(6.*mt2*RU11 + 3.0*At*mt*RU12 - 3.0*CoW2*mZ2*RU11*SiB2 + mZ2*RU11*SiB2*SiW2)))/SiB
        Ht11Hmx[2] = (CoI*mt*(-1.0*CRU12*mu*RU11 + Cmu*CRU11*RU12))/SiB
        Ht11Hmx[3] = (CoI*mt*(-1.0*CAt*CRU12*RU11 + At*CRU11*RU12))/SiB

        Ht12Hmx[0] = (1/3.0*(-1.0*CRU12*(3.0*mt*mu*RU21 - 4.*CoB*mZ2*RU22*SiB*SiW2) - 1.0*CRU11*(3.0*Cmu*mt*RU22 + CoB*mZ2*RU21*SiB*(-3.0*CoW2 + SiW2))))/SiB
        Ht12Hmx[1] = (1/3.0*(3.0*CAt*CRU12*mt*RU21 + 2.0*CRU12*RU22*(3.0*mt2 - 2.0*mZ2*SiB2*SiW2) + CRU11*(6.*mt2*RU21 + 3.0*At*mt*RU22 - 3.0*CoW2*mZ2*RU21*SiB2 + mZ2*RU21*SiB2*SiW2)))/SiB
        Ht12Hmx[2] = (CoI*mt*(-1.0*CRU12*mu*RU21 + Cmu*CRU11*RU22))/SiB
        Ht12Hmx[3] = (CoI*mt*(-1.0*CAt*CRU12*RU21 + At*CRU11*RU22))/SiB

        Ht21Hmx[0] = (1/3.0*(-1.0*CRU22*(3.0*mt*mu*RU11 - 4.*CoB*mZ2*RU12*SiB*SiW2) - 1.0*CRU21*(3.0*Cmu*mt*RU12 + CoB*mZ2*RU11*SiB*(-3.0*CoW2 + SiW2))))/SiB
        Ht21Hmx[1] = (1/3.0*(3.0*CAt*CRU22*mt*RU11 + 2.0*CRU22*RU12*(3.0*mt2 - 2.0*mZ2*SiB2*SiW2) + CRU21*(6.*mt2*RU11 + 3.0*At*mt*RU12 - 3.0*CoW2*mZ2*RU11*SiB2 + mZ2*RU11*SiB2*SiW2)))/SiB
        Ht21Hmx[2] = (CoI*mt*(-1.0*CRU22*mu*RU11 + Cmu*CRU21*RU12))/SiB
        Ht21Hmx[3] = (CoI*mt*(-1.0*CAt*CRU22*RU11 + At*CRU21*RU12))/SiB

        Ht22Hmx[0] = (1/3.0*(-1.0*CRU22*(3.0*mt*mu*RU21 - 4.*CoB*mZ2*RU22*SiB*SiW2) - 1.0*CRU21*(3.0*Cmu*mt*RU22 + CoB*mZ2*RU21*SiB*(-3.0*CoW2 + SiW2))))/SiB
        Ht22Hmx[1] = (1/3.0*(3.0*CAt*CRU22*mt*RU21 + 2.0*CRU22*RU22*(3.0*mt2 - 2.0*mZ2*SiB2*SiW2) + CRU21*(6.*mt2*RU21 + 3.0*At*mt*RU22 - 3.0*CoW2*mZ2*RU21*SiB2 + mZ2*RU21*SiB2*SiW2)))/SiB
        Ht22Hmx[2] = (CoI*mt*(-1.0*CRU22*mu*RU21 + Cmu*CRU21*RU22))/SiB
        Ht22Hmx[3] = (CoI*mt*(-1.0*CAt*CRU22*RU21 + At*CRU21*RU22))/SiB

        Hb11Hmx[0] = (1/3.0*(3.0*CAb*CRD12*mb*RD11 + 2.0*CRD12*RD12*(3.0*mb2 - 1.0*CoB2*mZ2*SiW2) + CRD11*(6.*mb2*RD11 - 3.0*CoB2*CoW2*mZ2*RD11 + 3.0*Ab*mb*RD12 - 1.0*CoB2*mZ2*RD11*SiW2)))/CoB
        Hb11Hmx[1] = (1/3.0*(CRD12*(-3.0*mb*mu*RD11 + 2.0*CoB*mZ2*RD12*SiB*SiW2) + CRD11*(-3.0*Cmu*mb*RD12 + CoB*mZ2*RD11*SiB*(3.0*CoW2 + SiW2))))/CoB
        Hb11Hmx[2] = (CoI*mb*(-1.0*CAb*CRD12*RD11 + Ab*CRD11*RD12))/CoB
        Hb11Hmx[3] = (CoI*mb*(-1.0*CRD12*mu*RD11 + Cmu*CRD11*RD12))/CoB

        Hb12Hmx[0] = (1/3.0*(3.0*CAb*CRD12*mb*RD21 + 2.0*CRD12*RD22*(3.0*mb2 - 1.0*CoB2*mZ2*SiW2) + CRD11*(6.*mb2*RD21 - 3.0*CoB2*CoW2*mZ2*RD21 + 3.0*Ab*mb*RD22 - 1.0*CoB2*mZ2*RD21*SiW2)))/CoB
        Hb12Hmx[1] = (1/3.0*(CRD12*(-3.0*mb*mu*RD21 + 2.0*CoB*mZ2*RD22*SiB*SiW2) + CRD11*(-3.0*Cmu*mb*RD22 + CoB*mZ2*RD21*SiB*(3.0*CoW2 + SiW2))))/CoB
        Hb12Hmx[2] = (CoI*mb*(-1.0*CAb*CRD12*RD21 + Ab*CRD11*RD22))/CoB
        Hb12Hmx[3] = (CoI*mb*(-1.0*CRD12*mu*RD21 + Cmu*CRD11*RD22))/CoB

        Hb21Hmx[0] = (1/3.0*(3.0*CAb*CRD22*mb*RD11 + 2.0*CRD22*RD12*(3.0*mb2 - 1.0*CoB2*mZ2*SiW2) + CRD21*(6.*mb2*RD11 - 3.0*CoB2*CoW2*mZ2*RD11 + 3.0*Ab*mb*RD12 - 1.0*CoB2*mZ2*RD11*SiW2)))/CoB
        Hb21Hmx[1] = (1/3.0*(CRD22*(-3.0*mb*mu*RD11 + 2.0*CoB*mZ2*RD12*SiB*SiW2) + CRD21*(-3.0*Cmu*mb*RD12 + CoB*mZ2*RD11*SiB*(3.0*CoW2 + SiW2))))/CoB
        Hb21Hmx[2] = (CoI*mb*(-1.0*CAb*CRD22*RD11 + Ab*CRD21*RD12))/CoB
        Hb21Hmx[3] = (CoI*mb*(-1.0*CRD22*mu*RD11 + Cmu*CRD21*RD12))/CoB

        Hb22Hmx[0] = (1/3.0*(3.0*CAb*CRD22*mb*RD21 + 2.0*CRD22*RD22*(3.0*mb2 - 1.0*CoB2*mZ2*SiW2)+ CRD21*(6.*mb2*RD21 - 3.0*CoB2*CoW2*mZ2*RD21 + 3.0*Ab*mb*RD22 - 1.0*CoB2*mZ2*RD21*SiW2)))/CoB
        Hb22Hmx[1] = (1/3.0*(CRD22*(-3.0*mb*mu*RD21 + 2.0*CoB*mZ2*RD22*SiB*SiW2) + CRD21*(-3.0*Cmu*mb*RD22 + CoB*mZ2*RD21*SiB*(3.0*CoW2 + SiW2))))/CoB
        Hb22Hmx[2] = (CoI*mb*(-1.0*CAb*CRD22*RD21 + Ab*CRD21*RD22))/CoB
        Hb22Hmx[3] = (CoI*mb*(-1.0*CRD22*mu*RD21 + Cmu*CRD21*RD22))/CoB

        Ht11Hmx = Ht11Hmx/mt2
        Ht12Hmx = Ht12Hmx/mt2
        Ht21Hmx = Ht21Hmx/mt2
        Ht22Hmx = Ht22Hmx/mt2

        Hb11Hmx = Hb11Hmx/mb2
        Hb12Hmx = Hb12Hmx/mb2
        Hb21Hmx = Hb21Hmx/mb2
        Hb22Hmx = Hb22Hmx/mb2

        if(self.lpseudo == 1):
            gth11 = -CoI*(SiB*Ht11Hmx[2] + CoB*Ht11Hmx[3])
            gth12 = -CoI*(SiB*Ht12Hmx[2] + CoB*Ht12Hmx[3])
            gth21 = -CoI*(SiB*Ht21Hmx[2] + CoB*Ht21Hmx[3])
            gth22 = -CoI*(SiB*Ht22Hmx[2] + CoB*Ht22Hmx[3])

            gbh11 = -CoI*(SiB*Hb11Hmx[2] + CoB*Hb11Hmx[3])
            gbh12 = -CoI*(SiB*Hb12Hmx[2] + CoB*Hb12Hmx[3])
            gbh21 = -CoI*(SiB*Hb21Hmx[2] + CoB*Hb21Hmx[3])
            gbh22 = -CoI*(SiB*Hb22Hmx[2] + CoB*Hb22Hmx[3])

        return {'gth11':gth11, 'gth12':gth12, 
            'gth21':gth21, 'gth22':gth22, 
            'gbh11':gbh11, 'gbh12':gbh12, 
            'gbh21':gbh21, 'gbh22':gbh22}

    def get_json(self):
        return json.dumps(self.get_coup(), cls=NumpyEncoder, indent=4)
