import sys
import numpy as np
import pyslha
from ast import literal_eval
from utils import *
from config import*
from argparse import Namespace

param_fp = Path(config_dir, "params.json")
params = Namespace(**load_dict(param_fp))

def set_precision(variable):
	return literal_eval("{:.12f}".format(variable))

islha_obj = pyslha.read(str(Path(input_dir, params.name)))
#islha_obj = pyslha.read(sys.argv[1])
lam = islha_obj.blocks['IBPPARA'][1]
eps = islha_obj.blocks['IBPPARA'][2]

squark = islha_obj.blocks['INTEGER'][1]

mA = set_precision(islha_obj.blocks["MASS"][1])
mt_ns = set_precision(islha_obj.blocks["MASS"][2])#*np.sqrt(1. - np.complex(0., lam))
mtmsbar = set_precision(islha_obj.blocks["MASS"][3])
mG = set_precision(islha_obj.blocks["MASS"][4])
msq1_ns = set_precision(islha_obj.blocks["MASS"][5])#*np.sqrt(1. - np.complex(0., lam))
msq2_ns = set_precision(islha_obj.blocks["MASS"][6])#*np.sqrt(1. - np.complex(0., lam))

mu = set_precision(islha_obj.blocks["REAL"][1])
At = set_precision(islha_obj.blocks["REAL"][2])
#At = set_precision(islha_obj.blocks['REAL'][9])
beta = set_precision(islha_obj.blocks["REAL"][3])
t = set_precision(islha_obj.blocks["REAL"][4])
a = set_precision(islha_obj.blocks["REAL"][5])
b = set_precision(islha_obj.blocks["REAL"][6])
corrf = set_precision(islha_obj.blocks['REAL'][7])

iters = islha_obj.blocks['VEGAS'][1]
startdiv = islha_obj.blocks['VEGAS'][2]
enddiv = islha_obj.blocks['VEGAS'][3]
startfin =islha_obj.blocks['VEGAS'][4]
endfin = islha_obj.blocks['VEGAS'][5]
alphafix = islha_obj.blocks['VEGAS'][6]
integrator = str(islha_obj.blocks['VEGAS'][8])

alphaqs = set_precision(islha_obj.blocks['QCD'][1])

tanb = np.tan(beta)


