import numpy as np
import pandas as pd
import lhapdf
import csv
from getdata import *
from const import *
from utils import writefile

from aggpdf import *
from virtuals import *

from scipy.integrate import quad
from runalphas import *
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)
lhapdf.setVerbosity(0)
# Get C data table, all input are in the array form

def delta(x, tau, muf, pdf):
	return pdf.xfxQ(21, x, muf)*pdf.xfxQ(21, tau/x, muf)/x

def delta_integ(tau, muf, pdf):
	return quad(delta, tau, 1, limit=1000, args=(tau,muf, pdf))

def sigmaBV(tau, mt, mb, mG, alphaggh, lpesudo, tanbresum, tantresum, muf, mur, pdf, norder):
	vec_delta_integ = np.vectorize(delta_integ)
	deltagg, error = vec_delta_integ(tau, muf, pdf)

	sigmalo, sigmanlo, errorsq = AMPQ(mt, mb, mG, lpesudo, tanbresum, tantresum, alphaggh, norder)
	sigmalo = sigmalo*alphaggh**2*sigmanullA*gev2picobarns*deltagg

	# relative error
	sigmanlo = np.real(sigmanlo)
	errorges = (errorsq*deltagg + error*sigmanlo)*alphaggh**2*sigmanullA
	sigmanlo = sigmanlo*alphaggh**2*sigmanullA
	sigmanlo = sigmanlo*deltagg*gev2picobarns
	errorges = errorges*gev2picobarns

	return sigmalo, sigmanlo, errorges

if __name__ == "__main__":

	params = getdata('summaryfile.csv')
	mt = 172.5000000000
	mb = 4.915765932164
	mG = 2500.000000000

	norder = 1
	lpesudo = True
	tanbresum = 0
	tantresum = 0

	# pdf
	pdflo = "MMHT2014lo68cl"
	pdfnlo = "PDF4LHC15_nlo_30"
	if (norder == 0):
		pdfnamein = pdflo
	elif (norder == 1):
		pdfnamein = pdfnlo

	pdf = lhapdf.mkPDF(pdfnamein)

	muf = params.mA/2.0
	mur = params.mA/2.0
	sqrtS = 13000.0
	S = sqrtS**2
	v = 246.0
	b0 = (33 - 2*nf)/12.0
	tau0 = params.mA**2/S
	beta = np.sqrt(1 - tau0)
	tol = 1e-10

	apimur = np.zeros(len(params.mA))
	apimz = alphasmzpdf(pdfnamein)
	nloop = 4
	for  idx, mur_idx in enumerate(mur):
		apimur[idx]= runalphas(apimz, mZ, mur_idx, nf, nloop, nsize)

	alphaggh = apimur*np.pi
	sigmalo, sigmanlo, error = sigmaBV(tau0, mt, mb, mG, alphaggh, lpesudo, tanbresum, tantresum, muf, mur, pdf, norder)

	fields = ['who', 'norder', 'tanbresum', 'tantresum', 'sigmaLO', 'sigmaNLO', 'errorNLO', 'mu', 'tanb', 'mA', 'mG', 'Ab', 'At', 'mst1', 'mst2', 'msb1', 'msb2', 'alphaqsb', 'alphaqst']
	who = 'KIT'
	writefile('xs_summary.csv', params, sigmalo, sigmanlo, error, norder, tanbresum, tantresum, mG)

