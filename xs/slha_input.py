import sys
import numpy as np
import pyslha
import warnings
from const import *

path = "/users/tp/dat/xs_full/"
warnings.filterwarnings('ignore', category=RuntimeWarning)
islha_obj = pyslha.read(path+'input.in')

# collider inputs
pesudo = islha_obj.blocks['SUSHI'][1]
ppbar = islha_obj.blocks['SUSHI'][2]
sqrtS = islha_obj.blocks['SUSHI'][3]
norder = islha_obj.blocks['SUSHI'][4]

# intial SM inputs
alpha_em = islha_obj.blocks['SMINPUTS'][1]
GF = islha_obj.blocks['SMINPUTS'][2]
alphasmz = islha_obj.blocks['SMINPUTS'][3]
mz = islha_obj.blocks['SMINPUTS'][4]
mbmb = islha_obj.blocks['SMINPUTS'][5]
mt = islha_obj.blocks['SMINPUTS'][6]
mcmc = islha_obj.blocks['SMINPUTS'][7]

tanb = islha_obj.blocks['MINPAR'][3]

# soft-susy parameter blocks
M3 = islha_obj.blocks['EXTPAR'][3]
At = islha_obj.blocks['EXTPAR'][11]
Ab = islha_obj.blocks['EXTPAR'][12]
muSUSY = islha_obj.blocks['EXTPAR'][23]
MQ3 = islha_obj.blocks['EXTPAR'][43]
MTR = islha_obj.blocks['EXTPAR'][46]
MBR =  islha_obj.blocks['EXTPAR'][49]

# Higgs parameter
alpha = islha_obj.blocks['ALPHA'][1]

mhc = islha_obj.blocks['MASS'][1]
mH = islha_obj.blocks['MASS'][2]
mA = islha_obj.blocks['MASS'][3]

# Distribution
ptcut = islha_obj.blocks['DISTRIB'][2]
rapcut = islha_obj.blocks['DISTRIB'][3]
minrap = islha_obj.blocks['DISTRIB'][31]
maxrap = islha_obj.blocks['DISTRIB'][32]

rap_pesudo = islha_obj.blocks['DISTRIB'][4] # choice between rapidity or pesudo-rapidity

# Scales
renormscale = islha_obj.blocks['SCALES'][1]
factorscale = islha_obj.blocks['SCALES'][2]

# renormalization of bottom
mbrunyuk = islha_obj.blocks['RENORMBOT'][1]
tanbresum = islha_obj.blocks['RENORMBOT'][2]

# renormalization sbottom
dep = islha_obj.blocks['RENORMSBOT'][1]

# PDFSEPC
pdflo = islha_obj.blocks['PDFSPEC'][1]
pdfnlo = islha_obj.blocks['PDFSPEC'][2]
pdfnnlo =  islha_obj.blocks['PDFSPEC'][3]

# Block Vegas
iters = islha_obj.blocks['VEGAS'][1]
nstart = islha_obj.blocks['VEGAS'][2]
nend =  islha_obj.blocks['VEGAS'][3]
alpha_vegas = islha_obj.blocks['VEGAS'][4]
eps_tol = islha_obj.blocks['VEGAS'][5]

# Block FACTORS
yukfacc = islha_obj.blocks['FACTORS'][1]
yukfact = islha_obj.blocks['FACTORS'][2]
yukfacb = islha_obj.blocks['FACTORS'][3]
yukfacst = islha_obj.blocks['FACTORS'][4]
yukfacsb = islha_obj.blocks['FACTORS'][5]

yukfac = np.zeros(9)
yukfac[0] = yukfacc # charm yukawa coupling
yukfac[1] = yukfact # top-yukawa coupling
yukfac[2] = yukfacb # bottom-yukawa coupling
yukfac[3] = 0.0 # top-yukawa coupling 4-gen
yukfac[4] = 0.0 # bottom-yukawa coupling 4-gen
yukfac[5] = yukfacst # stop-yukawa coupling
yukfac[6] = yukfacsb # sbottom-yukawa coupling
yukfac[7] = 0.0 # 4-gen-top mass
yukfac[8] = 0.0 # 4 gen-bottom mass


mb = mbmb
mb2 = mb**2
mt2 = mt**2
mc2 = mcmc**2

# running mb
mb = 4.904589444519764
mb2 = mb**2
# distribution
if(ptcut == 0):
    ptcut = False
else:
    ptcut = True

if(rapcut == 0):
    rapcut = False
else:
    rapcut = True

if (rap_pesudo == 0):
    # pick rapidity
    pesudorap = True
else:
    # pick pesudo_rapidity
    pesudorap = False

# order of calcuation
norderggh = norder

# after initcouplings
gt = 1.0/(tanb) 
gb = tanb
gc =0.0

mh = mA
mhiggs = mA
mh2 = mh**2
z = mh2/sqrtS**2

sigmanullh = GF/(288.0*np.sqrt(2.0)*np.pi)*gev2pb
sigmanullA = GF/(128.0*np.sqrt(2.0)*np.pi)*gev2pb

muRggh = mh*renormscale
muFggh = mh*factorscale

if (mbrunyuk == 1):
    muB = mbmb
elif (mbrunyuk == 2):
    muB = muRggh
else:
    muB = mh

lfh = np.log(muFggh**2/mh2)
lrh = np.log(muRggh**2/mh2)
lfr = np.log(muFggh**2/muRggh**2)
lrt = np.log(muRggh**2/mt**2)

subr = True
ptcut = False
rapcut = False
ppcoll = True