import lhapdf
import numpy as np
from slha_input import *

# this is part of Sushi code
# written in python by N.N.T.Dat

def intitalize_pdf():
    if(norderggh==0):
        pdfnamein = pdflo
    elif(norderggh==1):
        pdfnamein = pdfnlo
    
    p = lhapdf.mkPDF(pdfnamein)
    return p

def alphasmzpdf(pdfnamein):
    p = lhapdf.mkPDF(pdfnamein)
    apimz = p.alphasQ(mz)/np.pi
    return apimz

def pdfs(x, rmuf):
    p = intitalize_pdf()
    ff = [p.xfxQ(i, x, rmuf) for i in range(-5, 6)]
    return ff

def GetPDFs(x, rmuf):
    temp = pdfs(x, rmuf)
    bot = temp[10]
    chm = temp[9]
    stra = temp[8]
    up = temp[7]
    dn = temp[6]
    glu = temp[5]
    dsea = temp[4]
    usea = temp[3]
    sbar = temp[2]
    cbar = temp[1]
    bbar = temp[0]
    return up, dn, usea, dsea, stra, sbar, chm, cbar, bot, bbar, glu


def PDFgg(x1, x2, rmuf):
    upv1, dnv1, usea1, dsea1, str1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = GetPDFs(x1, rmuf)
    upv2, dnv2, usea2, dsea2, str2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = GetPDFs(x2, rmuf)

    PDFgg = glu1*glu2
    return PDFgg


def PDFqg(x1, x2, rmuf):
    upv1, dnv1, usea1, dsea1, str1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = GetPDFs(x1, rmuf)
    upv2, dnv2, usea2, dsea2, str2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = GetPDFs(x2, rmuf)

    q1 = upv1 + dnv1+usea1 + dsea1 + str1 + sbar1 + chm1 + cbar1 + bot1 + bbar1
    PDFqg = q1*glu2
    return PDFqg

def PDFqq(x1, x2, rmuf):
    upv1, dnv1, usea1, dsea1, str1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = GetPDFs(x1, rmuf)
    upv2, dnv2, usea2, dsea2, str2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = GetPDFs(x2, rmuf)

    if(ppcoll):
        PDFqq  =upv1*usea2+upv2*usea1+dnv1*dsea2+dnv2*dsea1+str1*sbar2+str2*sbar1+chm1*cbar2+chm2*cbar1+bot1*bbar2+bot2*bbar1
    else:
        PDFqq = upv1*upv2+usea1*usea2+dnv1*dnv2+dsea1*dsea2+ 2 * (str1*str2+chm1*chm2+bot1*bot2)
    return PDFqq