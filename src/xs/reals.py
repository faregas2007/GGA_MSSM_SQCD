# reals.py
import numpy as np

from typing import Dict, List
from pathlib import Path
from argparse import Namespace

from utils import *
from config import *

from xs.integrals import *
from xs.utils import *

from renschemes.renschemes import renschemes

class real_amps(renschemes):
    def __init__(self, params_fp: str=Path(model_dir, 'models.json'), **kwargs):
        
        self.params_fp = params_fp
        params = Namespace(**load_dict(params_fp))
        """
        self.mc = self._renschemes.get_params['mc']
        self.mb = self._renschemes.get_params['mb']
        self.mt = self._renschemes.get_params['mt']
        self.mh2 = self._renschemes.get_params['mh2']

        self.gc = self.renschemes.get_params['gc']
        self.gb = self.renschemes.get_params['gb']
        self.gt = self.renschemes.get_params['gt']

        self.mc2 = self.mc*self.mc
        self.mb2 = self.mb*self.mb
        self.mt2 = self.mt*self.mt
        """
        self.mc = params.mcos
        self.mb = params.mb
        self.mt = params.mt
        self.mh2 = params.mh2

        self.gc = params.gc
        self.gb = params.gb
        self.gt = params.gt

        self.mc2 = self.mc*self.mc
        self.mb2 = self.mb*self.mb
        self.mt2 = self.mt*self.mt

    # qq channel
    def AMPqq(self, s, tpu):
        return self.AMPqqpure(s, tpu)/self.AMPLO(s+tpu)
    
    def AMPqq_batch(self, s, tpu):
        return self.AMPqqpure_batch(s, tpu)/self.AMPLO_batch(s+tpu)

    def AMPqqpure(self, s, tpu):
        
        mc2 = self.mc2
        mb2 = self.mb2
        mt2 = self.mt2
        
        gc = self.gc
        gb = self.gb
        gt = self.gt

        A = 0.0
        if(np.abs(gc) != 0.0):
            A = A + mc2*gc*ASMA(s, tpu, mc2)
        
        if(np.abs(gb) != 0.0):
            A = A + mb2*gb*ASMA(s, tpu, mb2)
        
        if(np.abs(gt) != 0.0):
            A = A + mt2*gt*ASMA(s, tpu, mt2)
        
        return np.abs(A)**2

    def AMPqqpure_batch(self, s, tpu):
        mc2 = self.mc2
        mb2 = self.mb2
        mt2 = self.mt2
        
        gc = self.gc
        gb = self.gb
        gt = self.gt
        A = np.zeros(len(tpu),)
        if(np.abs(gc) != 0.0):
            A = A + mc2*gc*np.array(ASMA_batch(s, tpu, mc2))
        
        if(np.abs(gb) != 0.0):
            A = A + mb2*gb*np.array(ASMA_batch(s, tpu, mb2))
        
        if(np.abs(gt) != 0.0):
            A = A + mt2*gt*np.array(ASMA_batch(s, tpu, mt2))
        
        return np.abs(A)**2

    
    # qg channel
    def AMPqg(self, s, t, u):
        return self.AMPqgpure(s, t, u)/self.AMPLO(s+t+u)
    
    def AMPqg_batch(self, s, t, u):
        return self.AMPqgpure_batch(s, t, u)/self.AMPLO_batch(s+t+u)
    
    def AMPqgpure(self, s, t, u):
        mc2 = self.mc2
        mb2 = self.mb2
        mt2 = self.mt2

        gc = self.gc
        gb = self.gb
        gt = self.gt

        A = 0.0
        if(np.abs(gc) != 0.0):
            A = A + mc2*gc*ASMA(t, s+u, mc2)
        if(np.abs(gb) != 0.0):
            A = A + mb2*gb*ASMA(t, s+u, mb2)
        if(np.abs(gt) != 0.0):
            A = A + mt2*gt*ASMA(t, s+u, mt2)
        
        return -(s**2 + u**2)/t*np.abs(A)**2/(s+t+u)

    def AMPqgpure_batch(self, s, t, u):
        mc2 = self.mc2
        mb2 = self.mb2
        mt2 = self.mt2

        gc = self.gc
        gb = self.gb
        gt = self.gt

        A = np.zeros((len(u),))
        if(np.abs(gc) != 0.0):
            A = A + mc2*gc*ASMA_batch(t, s+u, mc2)
        if(np.abs(gb) != 0.0):
            A = A + mb2*gb*ASMA_batch(t, s+u, mb2)
        if(np.abs(gt) != 0.0):
            A = A + mt2*gt*ASMA_batch(t, s+u, mt2)
        
        return -(s**2 + u**2)/t*np.abs(A)**2/(s+t+u)
    
    # gg channel =
    def AMPgg(self, s, t, u):
        return self.AMPggpure(s,t,u)/self.AMPLO(s+t+u)

    def AMPgg_batch(self, s, t, u):
        return self.AMPggpure_batch(s,t,u)/self.AMPLO_batch(s+t+u)

    def AMPggpure(self, s, t, u):
        charm1, charm2, charm3, charm4 = complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0)
        bottom1, bottom2, bottom3, bottom4 = complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0)
        top1, top2, top3, top4 = complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0)
        c1, c2, c3, c4 = 0.0,0.0,0.0,0.0
    
        mc2 = self.mc2
        mb2 = self.mb2
        mt2 = self.mt2

        gc = self.gc
        gb = self.gb
        gt = self.gt

        if(np.abs(gc) != 0.0):
            charm1 = C1SMA(s,t,u,mc2)
            charm2, charm3, charm4 = C2SMA(s,t,u, mc2)
            c1 = c1 + mc2*gc*charm1 
            c2 = c2 + mc2*gc*charm2
            c3 = c3 + mc2*gc*charm3
            c4 = c4 + mc2*gc*charm4

        if(np.abs(gb) != 0.0):
            bottom1 = C1SMA(s,t,u,mb2)
            bottom2, bottom3, bottom4 = C2SMA(s,t,u, mb2)
            c1 = c1 + mb2*gb*bottom1 
            c2 = c2 + mb2*gb*bottom2
            c3 = c3 + mb2*gb*bottom3
            c4 = c4 + mb2*gb*bottom4

        if(np.abs(gt) != 0.0):
            top1 = C1SMA(s,t,u,mt2)
            top2, top3, top4 = C2SMA(s,t,u, mt2)
            c1 = c1 + mt2*gt*top1 
            c2 = c2 + mt2*gt*top2
            c3 = c3 + mt2*gt*top3
            c4 = c4 + mt2*gt*top4

        return (np.abs(c1)**2 + np.abs(c2)**2 + np.abs(c3)**2 + np.abs(c4)**2)/s/t/u*9/32.0/(s+t+u)

    def AMPggpure_batch(self, s, t, u):
        #s = np.array(s)
        #t = np.array(t)
        #u = np.array(u)
        
        charm1 = np.zeros(len(s), dtype=complex) 
        charm2 = np.zeros(len(s), dtype=complex) 
        charm3 = np.zeros(len(s), dtype=complex) 
        charm4 = np.zeros(len(s), dtype=complex) 
        bottom1 = np.zeros(len(s), dtype=complex) 
        bottom2 = np.zeros(len(s), dtype=complex) 
        bottom3 = np.zeros(len(s), dtype=complex) 
        bottom4 = np.zeros(len(s), dtype=complex) 
        top1 = np.zeros(len(s), dtype=complex) 
        top2 = np.zeros(len(s), dtype=complex) 
        top3 = np.zeros(len(s), dtype=complex) 
        top4 = np.zeros(len(s), dtype=complex) 

            
        c1 = np.zeros(len(s),)
        c2 = np.zeros(len(s),)
        c3 = np.zeros(len(s),)
        c4 = np.zeros(len(s),)
        
        mc2 = self.mc2
        mb2 = self.mb2
        mt2 = self.mt2

        gc = self.gc
        gb = self.gb
        gt = self.gt

        if(np.abs(gc) != 0.0):
            charm1 = C1SMA_batch(s,t,u,mc2)
            charm2, charm3, charm4 = C2SMA_batch(s,t,u, mc2)
            c1 = c1 + mc2*gc*charm1 
            c2 = c2 + mc2*gc*charm2
            c3 = c3 + mc2*gc*charm3
            c4 = c4 + mc2*gc*charm4

        if(np.abs(gb) != 0.0):
            bottom1 = C1SMA_batch(s,t,u,mb2)
            bottom2, bottom3, bottom4 = C2SMA_batch(s,t,u, mb2)
            c1 = c1 + mb2*gb*bottom1 
            c2 = c2 + mb2*gb*bottom2
            c3 = c3 + mb2*gb*bottom3
            c4 = c4 + mb2*gb*bottom4

        if(np.abs(gt) != 0.0):
            top1 = C1SMA_batch(s,t,u,mt2)
            top2, top3, top4 = C2SMA_batch(s,t,u, mt2)
            c1 = c1 + mt2*gt*top1 
            c2 = c2 + mt2*gt*top2
            c3 = c3 + mt2*gt*top3
            c4 = c4 + mt2*gt*top4

        return (np.abs(c1)**2 + np.abs(c2)**2 + np.abs(c3)**2 + np.abs(c4)**2)/s/t/u*9/32.0/(s+t+u)

    def AMPLO(self, mh2):
        gc = self.gc
        gb = self.gb
        gt = self.gt

        mc2 = self.mc2
        mb2 = self.mb2
        mt2 = self.mt2

        A = 0.0
        if(np.abs(gc) != 0.0):
            A = A + mc2*gc*ALOSMA(mh2, mc2)
        if(np.abs(gb) != 0.0):
            A = A + mb2*gb*ALOSMA(mh2, mb2)
        if(np.abs(gt) != 0.0):
            A = A + mt2*gt*ALOSMA(mh2, mt2)
        return np.abs(A)**2

    def AMPLO_batch(self, mh2):
        gc = self.gc
        gb = self.gb
        gt = self.gt

        mc2 = self.mc2
        mb2 = self.mb2
        mt2 = self.mt2

        if(np.array(mh2).shape == ()):
            mh2 = np.array(mh2)
            A = np.zeros(np.array(mh2).size,)
        else:
            A = np.zeros(mh2.size,)
        if(np.abs(gc) != 0.0):
            A = A + mc2*gc*ALOSMA_batch(mh2, mc2)
        if(np.abs(gb) != 0.0):
            A = A + mb2*gb*ALOSMA_batch(mh2, mb2)
        if(np.abs(gt) != 0.0):
            A = A + mt2*gt*ALOSMA_batch(mh2, mt2)
        return np.abs(A)**2