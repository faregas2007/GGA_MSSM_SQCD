import sys
from pathlib import Path
from argparse import Namespace

import pyslha
import numpy as np
from utils import *
from config import *
from ast import literal_eval

param_fp = Path(config_dir, "params.json")
params = Namespace(**load_dict(param_fp))
islha_obj = pyslha.read(str(Path(data_dir , params.name)))
squark = islha_obj.blocks['INTEGER'][1]
if(Micheal==True):
    lam = islha_obj.blocks["IBPPARA"][1]
    mA = islha_obj.blocks["MASS"][1] + np.complex(0.,0.)
    mt = islha_obj.blocks["MASS"][2]*np.sqrt(1. - np.complex(0., lam))
    mtmsbar = islha_obj.blocks["MASS"][3] + np.complex(0.,0.)
    mG = islha_obj.blocks["MASS"][4] + np.complex(0.,0.)
    msq1 = islha_obj.blocks["MASS"][5]*np.sqrt(1. - np.complex(0., lam))
    msq2 = islha_obj.blocks["MASS"][6]*np.sqrt(1. - np.complex(0., lam))
    mu = islha_obj.blocks["REAL"][1]
    At = islha_obj.blocks["REAL"][2] + np.complex(0.,0.)
    beta = islha_obj.blocks["REAL"][3] + np.complex(0.,0.)
    t = islha_obj.blocks["REAL"][4]
    a = islha_obj.blocks["REAL"][5] + np.complex(0.,0.)
    b = islha_obj.blocks["REAL"][6] + np.complex(0.,0.)
else:
    # parameters from renschemes of sushi
    # there is top/stop and bottom/sbottom input files. 
    # Since amplitude only takes a generic Aq, msq1, msq2, mtmsbar
    # Becareful with the Ab<>At coversion. 
    param_fp = Path(model_dr, "models.json")
    params = Namespace(**load_dict(params_fp))
    mA = params.mA
    mG = params.M3
    # modify here for bidirectional commands from coresushi and from sqcd.
    if(squark == 1):
        mt = params.mt * np.sqrt(1. - np.complex(0., lam))
        mtmsbar = params.mt
        At = params.At
        msq1 = np.sqrt(params.mst12) * np.sqrt(1. - np.complex(0., lam))
        msq2 = np.sqrt(params.mst22) * np.sqrt(1. - np.complex(0., lam))
    elif(squark==2):
        mt = params.mbos * np.sqrt(1. - np.complex(0., lam))
        mtmsbar = params.mbMSmuR
        At = params.Ab
        msq1 = np.sqrt(params.msb12) * np.sqrt(1. - np.complex(0., lam))
        msq2 = np.sqrt(params.msb22) * np.sqrt(1. - np.complex(0., lam))
    
    mu = params.mu
    beta = params.beta
    a = params.sthetab
    b = params.cthetab
