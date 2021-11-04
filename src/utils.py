import os, sys
import lhapdf
import numpy as np
import pandas as pd
import json, argparse

from typing import Dict, List

from ast import literal_eval
from pathlib import Path
from argparse import Namespace

from consts import *
from runalphas import *

def ifunc(a, b, c):
    return (a*b*np.log(a/b) + b*c*np.log(b/c) + a*c*np.log(c/a))/((a-b)*(b-c)*(a-c))

# should hide this into lhapdfs class
def SUSHI_alphas(alphasmz, rmurl, mZ):
	api0 = alphasmz/Pi
	apiout = runalphas(api0, mZ, rmurl, nf, 4, size)
	SUSHI_alphas = apiout*Pi
	#print(f'api0, apiout={api0},{apiout}')
	return SUSHI_alphas

def alphasmzpdf(pdfnamein, mZ):
    p = lhapdf.mkPDF(pdfnamein)
    return p.alphasQ(mZ)/Pi

def polemasshm(mqmq, apimq, nloop):
	lmum = np.float(0.0)
	mqpole = np.float(0.0)
	if(nloop == 0):
		mqpole = mqmq
	elif (nloop == 1):
		mqpole =  mqmq*( 1 + apimq*cf*(1 + (3*lmum)/4.0) )
	elif(nloop == 2):
		mqpole =  mqmq*( 1 + apimq*cf*(1 + (3*lmum)/4.0) + apimq**2*(cf**2*( -71/128.0 - (9*lmum)/32.0 + (9*lmum**2)/32.0 +(15*z2)/8.0  - 3*ln2*z2 + (3*z3)/4.0 )+ cf*ca*( 1111/384.0 + (185*lmum)/96.0 + (11*lmum**2)/32.0  - z2/2.0 + (3*ln2*z2)/2.0 - (3*z3)/8.0 )+ cf*tr*( -3/4.0 - (71*nf)/96.0 - (13*lmum*nf)/24.0 -(lmum**2*nf)/8.0 + (3*z2)/2.0 - (nf*z2)/2.0) ) )
	elif(nloop == 3):
		#mqpole =  mqmq*( 1 + apimq*cf*(1 + (3*lmum)/4.0) + apimq**2*(cf**2*( -71/128.0 - (9*lmum)/32.0 + (9*lmum**2)/32.0 +(15*z2)/8.0  - 3*ln2*z2 + (3*z3)/4.0 )+ cf*ca*( 1111/384.0 + (185*lmum)/96.0 + (11*lmum**2)/32.0  - z2/2.0+ (3*ln2*z2)/2.0 - (3*z3)/8.0 )+ cf*tr*( -3/4.0 - (71*nf)/96.0 - (13*lmum*nf)/24.0 -(lmum**2*nf)/8.0 + (3*z2)/2.0 - (nf*z2)/2.0) ) + apimq**3*(0.65270*(nf-1.0)**2 - 26.6550*(nf-1.0)+190.5950) ) 
		mqpole = mqmq*( 1 + apimq*cf*(1 + (3*lmum)/4.0)
			+ apimq**2*(
			+ cf**2*( -71/128.0 - (9*lmum)/32.0 + (9*lmum**2)/32.0
			+(15*z2)/8.0  - 3*ln2*z2 + (3*z3)/4.0 )+ cf*ca*( 1111
			/384.0 + (185*lmum)/96.0 + (11*lmum**2)/32.0  - z2/2.0
			+ (3*ln2*z2)/2.0 - (3*z3)/8.0 )+ cf*tr*( -3/4.0 - (71
			*nf)/96.0 - (13*lmum*nf)/24.0 -(lmum**2*nf)/8.0 + (3*z2
			)/2.0 - (nf*z2)/2.0) )
			+ apimq**3*(.65270*(nf-1.0)**2 - 26.6550*(nf-1.0)
			+190.5950) ) 
	else:
		print(f"<function polemass>: nloop={nloop} not implemented.")
	return mqpole

def getmass(alphasmz, muRggh, mZ, mcmc, mbmb):
	apimuR = SUSHI_alphas(alphasmz, muRggh, mZ)/Pi

	apim = SUSHI_alphas(alphasmz, mcmc, mZ)/Pi
	mcos = polemasshm(mcmc, apim, 3)

	apim = SUSHI_alphas(alphasmz, mbmb, mZ)/Pi
	mbos = polemasshm(mbmb, apim, 3)

	# mass mbmuR
	mbMSbarmuR = runmass(mbmb, apim, apimuR, nf, 4)

	return mcos, mbos, mbMSbarmuR

def Omega(a, b):
	Cbig = 0.5*(a+b) - 0.25*(a-b)*(a-b) - 0.25
	if (Cbig >= 0.0):
		sqCbig = np.sqrt(Cbig)
		Omega = sqCbig*(np.arctan((1 - a + b)/(2*sqCbig)) + np.arctan((1 - a + b)/(2*sqCbig)))
	elif (Cbig <= 0.0):
		sqCbig = np.sqrt(-Cbig)
		Omega = 0.5*sqCbig*np.log((a + b - 1.0 - 2.0*sqCbig)/(a + b - 1.0 + 2.0*sqCbig))
	else:
		Omega = 0.0
	return Omega

def b0fin(q, m1, m2, mu):
	"""
	used the b0fin in evalsusy.f
	change to new routine in ddsgluod.f
	only take the real part
	"""
	sm1 = m1**2
	sm2 = m2**2
	mu2 = mu**2
	if (q == 0.0):
		if (sm1 == 0.0 and sm2 != 0.0):
			b0fin = 1.0 - np.log(sm2/mu2)
		elif (m1 != 0.0 and m2 == 0.0):
			b0fin = 1.0 - np.log(sm1/mu2)
		elif (np.abs(sm1 - sm2) <= 10**-8):
			b0fin = -np.log(sm1/sm2)
		else:
			b0fin = 1.0 - np.log(sm2/mu2) + sm1/(sm1 - sm2)*np.log(sm2/sm1)
	else:
		if (sm1 == 0.0 and sm2 != 0.0):
			if (sm2 != q):
				b0fin = -(np.log(sm2/mu2) - 2.0 - (sm2/q - 1.0)*np.log(np.abs(1.0 - q/sm2)))
			else:
				b0fin = -(np.log(sm2/mu2) - 2.0)
		elif (sm2 == 0.0 and sm1 != 0.0):
			if (sm1 != q):
				b0fin = -(np.log(sm1/mu2) - 2.0 - (sm1/q - 1.0)*np.log(np.abs(1.0 - q/sm1)))
			else:
				b0fin = -(np.log(sm1/mu2) - 2.0)
		elif (sm2 == 0.0 and sm1 == 0.0):
			b0fin = - (np.log(q/mu2) - 2.0)
		else:
			b0fin = -(np.log(q/mu2) - 2.0 + 1.0/2.0*(1.0 + (sm1/q - sm2/q))*np.log(sm1/q) + 1.0/2.0*(1.0 - (sm1/q -sm2/q))*np.log(sm2/q) + 2.0*Omega(sm1/q, sm2/q))

	if (q <= (m1 + m2)**2):
		Imb0fin = 0.0
	else:
		Imb0fin = np.pi*np.sqrt(q - (m1 + m2)**2)*np.sqrt(q - (m1 - m2)**2)/q
	return b0fin

def set_precision(variable):
	return literal_eval("{:.12f}".format(variable))


def load_dict(filepath: str)->Dict:
	with open(filepath, 'r') as fp:
		d = json.load(fp)
	return d


def save_dict(d: Dict, filepath: str, cls=None, sortkeys: bool=False)->None:
	with open(filepath, 'w') as fp:
		json.dump(d, indent=2, fp=fp, cls=cls, sort_keys=sortkeys)

def parse_arguments():
    parser = argparse.ArgumentParser(description='gga_mssm_sqcd')
    parser.add_argument('-id', '--index_in', type=lambda x: int(float(x)), metavar='', help='diagram index',)
    parser.add_argument('-lam', '--lamb', type=lambda x: float(x), metavar='', help='lambda squark mass shifted',)

    parser.add_argument('-in', '--input', type=argparse.FileType('r'), metavar='INPUT_FILE', help='input directory',)
    parser.add_argument('-amp', '--amp_dir', type=dir_path, metavar='AMP_DIR', help='amplitude directory',)

    return parser.parse_args()

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentDefaultsHelpFormatter(f'readable_dir:{path} is not a valid path')


def combined_csv(allfiles, filename):
    combined = pd.concat([pd.read_csv(f) for f in allfiles])
    result = combined.sort_values(by=['lam'], ascending=True)
    result.to_csv(filename, index=False, encoding='utf-8-sig')
    return 0