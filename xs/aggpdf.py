import lhapdf
import numpy as np
from const import *

def alphasmzpdf(pdfnamein):
	pdf = lhapdf.mkPDF(pdfnamein)
	apimz = pdf.alphasQ(mZ)/np.pi
	return apimz

def pdfs(pdfnamein, sign, xx, rmuf):
	ff = np.zeros(11)
	pdf = lhapdf.mkPDF(pdfnamein)
	for i in range(-5,6):
		ff[i] = pdf.xfxQ(sign*i, xx, rmuf)/xx

	return ff
