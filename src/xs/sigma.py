import sys, logging
import numpy as np

from typing import Dict, List
from argparse import Namespace

from consts import *
from config import *
from init import *
from lhapdfs import *
from interface import *

from xs.reals import *
from xs.integrations import *

from renschemes.renschemes import *

class sigma_batch:
    def __init__(self, 
        qqinteg_batch, 
        qginteg_batch, 
        gginteg_batch):
        params = Namespace(**load_dict(Path(model_dir, 'models.json')))
        self.alphasggh = params.alphasggh
        self.norder = params.norder
        self.mh2 = params.mh2
        self.sigmanull = params.sigmanull
        
        self._qq = qqinteg_batch()
        self._qg = qginteg_batch()
        self._gg = gginteg_batch()

    def xs(self): 
        sigma = np.array([0.0,0.0])
        error = np.array([0.0,0.0])
        chi2 = np.array([0.0,0.0])
        alphasggh = self.alphasggh
        sigmanull = self.sigmanull
        norder = self.norder
        mh2 = self.mh2

        # logger will move into coresushi for matainable code. 
        # problem with transplus in batch mode --> large contributions qg and gg.
        logger.info('xs batch evaluation starts ...')
        F0A = real_amps().AMPLO_batch(mh2)
        logger.info(f'LO: {F0A}')
        
        if(self.norder == 0):
            logger.info(f'order of computation: LO')
            ggdelta = ggbatch(dim=1)
            delta, error[0], chi2[0] = integrators.integrate(integrand=ggdelta, ndim=ggdelta.dim)
            sigma[0] = delta*alphasggh*alphasggh
            error[0] = error[0]*alphasggh*alphasggh*sigmanull*F0A
            logger.info(f'alphasggh, delta, sigmaLO, errorLO: {alphasggh, delta, sigma[0], error[0]}')
            logger.info(f'xs evaluation ends ...')
            return {'norder': norder, 'sigmaLO':sgima[0], 'errorLO':error[0], 'sigmaNLO':sigma[1], 'errorNLO':error[1]}
    
        elif(self.norder == 1):
            logger.info(f'order of computation: NLO')
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
            logger.info(f'gg,qg,qq channels: {gg, errorgg}, {qg, errorqg}, {qq, errorqq}')

            sigma[0] = delta*alphasggh*alphasggh*sigmanull*F0A
            error[0] = errord*alphasggh*alphasggh*sigmanull*F0A

            sigma[1] = (gg + qg + qq)
            error[1] =  np.sqrt(errorgg*errorgg + errorqq*errorqq + errorqg*errorqg)
            logger.info(f'alphasggh, delta, sigmaLO, errorLO, sigmaNLO, errorNLO: {alphasggh, delta, sigma[0], error[0], sigma[1], error[1]}')
            logger.info(f'xs batch evaluation ends ...')
            return {'norder': norder, 'sigmaLO': sigma[0], 'errorLO': error[0], 'sigmaNLO': sigma[1], 'errorNLO':error[1]}
    
    def get_json(self)->Dict:
        return json.dumps(self.xs(), cls=NumpyEncoder, indent=4)

    def to_file(self):
        return save_dict(self.xs(), Path(model_dir, 'xs.json'), cls=NumpyEncoder)



class sigma:
    def __init__(self, 
        qq: integ, 
        qg: integ, 
        gg: integ):
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
        norder = self.norder
        mh2 = self.mh2

        # logger will move into coresushi
        # problem with transplus in batch mode --> numerical instabilities 
        logger.info('xs evaluation starts ...')
        F0A = real_amps().AMPLO_batch(mh2)
        logger.info(f'LO: {F0A}')
        
        if(self.norder == 0):
            logger.info(f'order of computation: LO')
            ggdelta = ggbatch(dim=1)
            delta, error[0], chi2[0] = integrators.integrate(integrand=ggdelta, ndim=ggdelta.dim)
            sigma[0] = delta*alphasggh*alphasggh
            error[0] = error[0]*alphasggh*alphasggh*sigmanull*F0A
            logger.info(f'alphasggh, delta, sigmaLO, errorLO: {alphasggh, delta, sigma[0], error[0]}')
            logger.info(f'xs evaluation ends ...')
            return {'norder': norder, 'sigmaLO':sgima[0], 'errorLO':error[0], 'sigmaNLO':sigma[1], 'errorNLO':error[1]}
    
        elif(self.norder == 1):
            logger.info(f'order of computation: NLO')
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
            logger.info(f'gg,qg,qq channels: {gg, errorgg}, {qg, errorqg}, {qq, errorqq}')

            sigma[0] = delta*alphasggh*alphasggh*sigmanull*F0A
            error[0] = errord*alphasggh*alphasggh*sigmanull*F0A

            sigma[1] = (gg + qg + qq)
            error[1] =  np.sqrt(errorgg*errorgg + errorqq*errorqq + errorqg*errorqg)
            logger.info(f'alphasggh, delta, sigmaLO, errorLO, sigmaNLO, errorNLO: {alphasggh, delta, sigma[0], error[0], sigma[1], error[1]}')
            logger.info(f'xs evaluation ends ...')
            return {'norder': norder, 'sigmaLO': sigma[0], 'errorLO': error[0], 'sigmaNLO': sigma[1], 'errorNLO':error[1]}
    
    def get_json(self)->Dict:
        return json.dumps(self.xs(), cls=NumpyEncoder, indent=4)

    def to_file(self):
        return save_dict(self.xs(), Path(model_dir, 'xs.json'), cls=NumpyEncoder)


class sigma_sqcd:
    def __init__(self, 
        qq: integ, 
        qg: integ, 
        gg: integ):
        params = Namespace(**load_dict(Path(model_dir, 'models.json')))
        self.alphasggh = params.alphasggh
        self.norder = params.norder
        self.mh2 = params.mh2
        self.sigmanull = params.sigmanull
        
        self._qq = qq()
        self._qg = qg()
        self._gg = gg()

    def xs(self): 
        sigma0, sigma1 = 0.0, 0.0
        error0, error1 = 0.0, 0.0
        chi2 = np.array([0.0,0.0])
        alphasggh = self.alphasggh
        sigmanull = self.sigmanull
        norder = self.norder
        mh2 = self.mh2

        logger.info('xs sqcd evaluation starts ...')
        F0A = real_amps().AMPLO_batch(mh2)
        logger.info(f'LO: {F0A}')
        
        if(self.norder == 0):
            logger.info(f'order of computation: LO')
            ggdelta = ggbatch(dim=1)
            #delta, error[0], chi2[0] = integrators.integrate(integrand=ggdelta, ndim=ggdelta.dim)
            delta, error0, chi20 = integrators.integrate(integrand=ggdelta, ndim=ggdelta.dim)
            sigma0 = delta*alphasggh*alphasggh
            #error[0] = error[0]*alphasggh*alphasggh*sigmanull*F0A
            error0 = error0*alphasggh*alphasggh*sigmanull*F0A
            logger.info(f'alphasggh, delta, sigmaLO, errorLO: {alphasggh, delta, sigma0, error[0]}')
            logger.info(f'xs sqcd evaluation ends ...')
            #return {'norder': norder, 'sigmaLO':sigma0, 'errorLO':error[0], 'sigmaNLO':sigma1, 'errorNLO':error[1]}
            return {'norder':norder, 'sigmaLO': sigmaLO, 'errorLO':error0, 'sigmaNLO':sigma1, 'errorNLO':error1}

        elif(self.norder == 1):
            logger.info(f'order of computation: NLO')
            gg ,delta, errorr, errorv, errord = self._gg.gginteg_sqcd()  
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
            logger.info(f'gg,qg,qq channels: {gg, errorgg}, {qg, errorqg}, {qq, errorqq}')

            sigma0 = delta*alphasggh*alphasggh*sigmanull*F0A
            #error[0] = errord*alphasggh*alphasggh*sigmanull*F0A
            error0 = errord*alphasggh*alphasggh*sigmanull*F0A

            sigma1 = gg + qg + qq
            #error[1] =  np.sqrt(errorgg*errorgg + errorqq*errorqq + errorqg*errorqg)
            error1 = np.sqrt(errorgg*errorgg + errorqq*errorqq + errorqg*errorqg)
            logger.info(f'alphasggh, delta, sigmaLO, errorLO, sigmaNLO, errorNLO: {alphasggh, delta, sigma0, error0, sigma1, error1}')
            logger.info(f'xs sqcd evaluation ends ...')
            return {'norder': norder, 'sigmaLO': sigma0, 'errorLO': error0, 'sigmaNLO': sigma1, 'errorNLO':error1}