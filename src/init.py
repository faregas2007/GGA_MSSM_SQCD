# init.py
#from __future__ import annotations

from argparse import Namespace
from numpyencoder import NumpyEncoder
from pathlib import Path

from abc import ABC, abstractmethod
from typing import Dict, List
import numpy as np
import pyslha
import lhapdf
import json

from config import *
from consts import *
from utils import *
from runalphas import *

class inputs:
    def __init__(self,
        params_fp: str=data_dir):
    
        self.params_fp = params_fp
        islha_obj = pyslha.read(str(Path(self.params_fp, "MSSM_SLHA.in")))
        
        self.model = islha_obj.blocks['SUSHI'][1]
        self.pseudo = islha_obj.blocks['SUSHI'][2]
        self.ppbar = islha_obj.blocks['SUSHI'][3]
        self.sqrtS = islha_obj.blocks['SUSHI'][4]
        self.norder = islha_obj.blocks['SUSHI'][5]

        self.alpha_em = islha_obj.blocks['SMINPUTS'][1]
        self.GF = islha_obj.blocks['SMINPUTS'][2]
        self.alphasmz = islha_obj.blocks['SMINPUTS'][3]
        self.mZ = islha_obj.blocks['SMINPUTS'][4]
        self.mW = islha_obj.blocks['SMINPUTS'][5]
        self.mbmb = islha_obj.blocks['SMINPUTS'][6]
        self.mbpole = islha_obj.blocks['SMINPUTS'][7]
        self.mt = islha_obj.blocks['SMINPUTS'][8]
        self.mcmc = islha_obj.blocks['SMINPUTS'][9]

        #self.size = islha_obj.blocks['MINPAR'][1]
        self.tanb = islha_obj.blocks['MINPAR'][3]

        self.M3 = islha_obj.blocks['EXTPAR'][3]
        self.At = islha_obj.blocks['EXTPAR'][11]
        self.Ab = islha_obj.blocks['EXTPAR'][12]
        self.mu = islha_obj.blocks['EXTPAR'][23]
        self.MQ3 = islha_obj.blocks['EXTPAR'][43]
        self.MTR = islha_obj.blocks['EXTPAR'][46]
        self.MBR = islha_obj.blocks['EXTPAR'][49]

        self.alpha = islha_obj.blocks['ALPHA'][1]

        self.mh = islha_obj.blocks['MASS'][1]
        self.mH = islha_obj.blocks['MASS'][2]
        self.mA = islha_obj.blocks['MASS'][3]

        self.pdflo = islha_obj.blocks['PDFSPEC'][1]
        self.pdfnlo = islha_obj.blocks['PDFSPEC'][2]
        self.pdfnnlo = islha_obj.blocks['PDFSPEC'][3]

        # distribution
        self.distrb = islha_obj.blocks['DISTRIB'][1]
        self.ptcut = islha_obj.blocks['DISTRIB'][2]
        self.rapcut = islha_obj.blocks['DISTRIB'][3]
        self.minrap = islha_obj.blocks['DISTRIB'][31]
        self.maxrap = islha_obj.blocks['DISTRIB'][32]
        self.minptc = islha_obj.blocks['DISTRIB'][21]
        self.maxptc = islha_obj.blocks['DISTRIB'][22]
        self.rap_c = islha_obj.blocks['DISTRIB'][4]

        # scales
        self.muRfacggh = islha_obj.blocks['SCALES'][1]
        self.muFfacggh = islha_obj.blocks['SCALES'][2]

        # renromalization of bottom
        self.mbrunyuk = islha_obj.blocks['RENORMBOT'][1]
        self.tanbresum = islha_obj.blocks['RENORMBOT'][2]

        # renormalization of sbottom
        self.mbsbch = islha_obj.blocks['RENORMSBOT'][1]
        self.Abch = islha_obj.blocks['RENORMSBOT'][2]
        self.thetach = islha_obj.blocks['RENORMSBOT'][3]

        # Block Vegas
        self.iters = islha_obj.blocks['VEGAS'][1]
        self.nstart = islha_obj.blocks['VEGAS'][2]
        self.nend = islha_obj.blocks['VEGAS'][3]
        self.alpha_vegas = islha_obj.blocks['VEGAS'][4]
        self.eps_tol = islha_obj.blocks['VEGAS'][5]

        # Block FACTORS
        self.yukfacc = islha_obj.blocks['FACTOR'][1]
        self.yukfact = islha_obj.blocks['FACTOR'][2]
        self.yukfacb = islha_obj.blocks['FACTOR'][3]
        self.yukfacst = islha_obj.blocks['FACTOR'][4]
        self.yukfacsb = islha_obj.blocks['FACTOR'][5]

        self.yukfac = np.zeros(9)
        self.yukfac[0] = self.yukfacc
        self.yukfac[1] = self.yukfact
        self.yukfac[2] = self.yukfacb
        self.yukfac[3] = 0.0
        self.yukfac[4] = 0.0
        self.yukfac[5] = self.yukfacst
        self.yukfac[6] = self.yukfacsb
        self.yukfac[7] = 0.0
        self.yukfac[8] = 0.0

        self.Hmx = np.zeros((2,)*2)
        self.Amx = np.zeros((2,)*2)
        self.CHmx = np.zeros((2,)*2, dtype=complex)
        
        # default value in MSSM
        # will move to initlize
        self.Amx[1,1] = 1.0

class einital(inputs):
    def initial(self):
        if(self.pseudo == 21):
            Sind = 1
            lpseudo = 1
            Mh = self.mA

        lam = 1
        muRggh = Mh*self.muRfacggh
        muFggh = Mh*self.muFfacggh
        
        # always set to zero in sushi
        mbrunloop = 0
        # must be in config file
        subtr= True
        ptcut= False
        rapcut= False    

        if(self.ppbar == 0):
            ppcoll = True
        else:
            ppcoll = False
        
        if(self.rap_c == 0):
            pseudorap = True
        else:
            pseudorap = False 

        mh2 = Mh*Mh
        z = mh2/self.sqrtS**2
        sigmanull = self.GF/288.0 / np.sqrt(2.0)/Pi * gev2pb

        if(self.distrb == 0):
            maxrap = -np.log(z)/2.0
            minrap = 0.0
        else:
            maxrap = self.maxrap
            minrap = self.minrap

        lfh = np.log(muFggh**2/mh2)
        lrh = np.log(muRggh**2/mh2)
        lfr = np.log(muFggh**2/muRggh**2)

        if(self.norder == 0):
            pdfnamein = self.pdflo
        elif(self.norder == 1):
            pdfnamein = self.pdfnlo
        else:
            pdfnamein = self.pdfnnlo
        
        apimz = alphasmzpdf(pdfnamein, self.mZ)
        apimur = runalphas(apimz, self.mZ, muRggh, nf, self.norder+1, size)
        alphasggh = apimur*Pi

        mcos, mbos, mbMSbarmuR = getmass(self.alphasmz, muRggh, self.mZ, self.mcmc, self.mbmb)
        
        return {
            'Sind':Sind, 'lpseudo':lpseudo, 'lam':lam, 
            'muRggh':muRggh, 'muFggh':muFggh,
            'subtr':subtr, 'ptcut':ptcut, 'rapcut':rapcut,
            'z':z, 'sigmanull':sigmanull, 
            'lfh':lfh, 'lrh':lrh, 'lfr':lfr,
            'pdfnamein':pdfnamein,
            'mcos':mcos, 'mbos':mbos, 'mbMSbarmuR':mbMSbarmuR,
            'alphasggh':alphasggh, 'mbrunloop':mbrunloop, 'mh2':mh2,
            'maxrap':maxrap, 'minrap':minrap, 'pseudorap':pseudorap, 
            'ppcoll':ppcoll, 'norder':self.norder, 'mgl':self.M3,
            'mu':self.mu, 'tanb':self.tanb, 'mA':self.mA
        }

    def get_json(self):
        return json.dumps(self.initial(), cls=NumpyEncoder, indent=4)

