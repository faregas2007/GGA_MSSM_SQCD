import numpy as np
from argparse import Namespace

from xs.reals import *
from xs.integrations import *

from renschemes.renschemes import *

from consts import *
from init import *
from lhapdfs import *

# facade design pattern or use params_fp path to model directory ?
# facade desgin pattern with decompose subclasses and params are called via params_fp in namespace.
class sigma:
    def __init__(self, qq, qg, gg):
        params = Namespace(**load_dict(Path(model_dir, 'models.json')))
        self.alphasggh = params.alphasggh
        self.norder = params.norder
        self.mh2 = params.mh2
        self.sigmanull = params.sigmanull
        
        self._qq = qq()
        self._qg = qg()
        self._gg = gg()

    def xs(self): 
        sigma = np.array([0.0,0.0])
        error = np.array([0.0,0.0])
        chi2 = np.array([0.0,0.0])
        alphasggh = self.alphasggh
        sigmanull = self.sigmanull
        mh2 = self.mh2

        F0A = real_amps().AMPLO(mh2)
        if(self.norder == 0):
            delta, error[0], chi2[0] = self._gg.deltas()
            sigma[0] = delta*alphasggh*alphasggh
            error[0] = error[0]*alphasggh*alphasggh*sigmanull*F0A
        
            return {'norder': norder, 'sigmaLO':sgima[0], 'errorLO':error[0]}
    
        elif(self.norder == 1):
            gg ,delta, errorr, errorv, errord = self._gg.gginteg()   
            qg, errorqg = self._qg.qginteg()
            qq, errorqq = self._qq.qqinteg()

            prefac = (alphasggh**3/Pi)*sigmanull * F0A

            gg = (delta*Pi/alphasggh + gg)*prefac
            qg = qg * prefac
            qq = qq * prefac

            # error
            errorgg = np.sqrt(errorr**2 + errorv**2 + (Pi/alphasggh*errord)**2)*prefac
            errorqg = errorqg*prefac
            errorqq = errorqq*prefac

            sigma[0] = delta*alphasggh*alphasggh*sigmanull*F0A
            error[0] = errord*alphasggh*alphasggh*sigmanull*F0A

            sigma[1] = (gg + qg + qq)
            error[1] =  np.sqrt(errorgg*errorgg + errorqq*errorqq + errorqg*errorqg)

            return {'norder': norder, 'sigmaLO': sigma[0], 'errorLO': error[0], 'sigmaNLO': sigma[1], 'errorNLO':error[1]}
    
    def get_json(self)->Dict:
        return json.dumps(self.xs(), cls=NumpyEncoder, indent=4)

    def to_file(self):
        return save_dict(self.xs(), Path(model_dir, 'xs.json'), cls=NumpyEncoder)
    