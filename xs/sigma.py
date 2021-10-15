# this is part of Sushi
# written in python by N.T.T.Dat

from integration import *
from integrals import *
from reals import *
from virtuals import *
from utils import *
from slha_input import *
from const import *
from runalphas import *
from lhapdfs import *

def xs(norderggh):
    
    mtop = mt
    mh2 = mh**2
    F02 = AMPLO(mh2)

    gg = 0.0
    qg = 0.0
    qq = 0.0
    error = [0.0,0.0]
    sigma = [0.0,0.0]
    # only focus on inclusive cross-section
    # leading order
    if (norderggh == 0):
        # will move into main.py
        pdfnamein = pdflo

        apimz = alphasmzpdf(pdfnamein)
        apimur = runalphas(apimz, mz, muRggh, nf,0, size)
        alphaggh = apimur*np.pi

        delta, error[0], chi2 = integrate(ndim, eps_tol,ggdelta, iters, nstart, nend, alpha_vegas)
        sigma[0] = delta*alphaggh**2*sigmanullA*F02
        error[0] = error[0]*alphaggh**2*sigmanullA*F02

        return sigma[0], error[0]

    elif (norderggh == 1):
        # next-to leading order
        # call alphas
        # will move into main.py
        pdfnamein = pdfnlo
        
        apimz = alphasmzpdf(pdfnamein)
        apimur = runalphas(apimz, mz, muRggh, nf, 3, size)
        alphaggh = apimur*np.pi

        # parallel here ?
        gg, delta, errorr, errorv, errord = intgg(ggdelta, ggsub, ggreal)
        qg, errorqg = intqg(qgreal, qgsub)
        qq, errorqq = intqq(qqreal)

        prefac = (alphaggh**3/np.pi) * sigmanullA * F02

        gg = (delta*np.pi/alphaggh + gg)*prefac
        qg = qg * prefac
        qq = qq * prefac

        # error
        errorgg = np.sqrt(errorr**2 + errorv**2 + (np.pi/alphaggh*errord)**2)*prefac 
        errorqg = errorqg*prefac
        errorqq = errorqq*prefac

        sigma[0] = delta*alphaggh**2*sigmanullA*F02
        error[0] = errord*alphaggh**2*sigmanullA*F02

        sigma[1] = (gg + qg + qq)
        error[1] = np.sqrt(errorgg**2 + errorqq**2 + errorqg**2)

        return sigma[0], error[0], sigma[1], error[1]

def intgg(ggdelta, ggsub, ggreal):

    error = [0.0,0.0,0.0]
    chi2 = [0.0,0.0,0.0]
    #virtual, virtual_susy, errorv_susy = c_virt()
    virtual = c_virt()

    delta, error[0], chi2[0] = integrate(1, eps_tol, ggdelta, iters, nstart, nend, alpha_vegas)
    sub, error[2], chi2[2] = integrate(2, eps_tol, ggsub, iters, nstart, nend, alpha_vegas)
    hard, error[1], chi2[1] = integrate(3, eps_tol, ggreal, iters, nstart, nend, alpha_vegas)
    
    errord = error[0]

    reell = ca*(hard+sub)
    errorr = ca*np.sqrt(error[1]**2 + error[2]**2)
    
    if (subtr):
        #virt = delta*(virtual+ virtual_susy +  np.pi**2-2*betanull*lfr + ca*dtermsz())
        virt = delta*(virtual + np.pi**2-2*betanull*lfr + ca*dtermsz())
        # relative error
        errorv = error[0]*(virtual+np.pi**2 - 2*betanull*lfr + ca*dtermsz())
        #errorv_v = errorv +  errorv_susy*delta
    else:
        virt = 0.0
        errorv = 0.0

    intgg = reell + virt
    return intgg, delta, errorr, errorv, errord

def intqg(qgreal, qgsub):
    error = [0.0,0.0]
    chi2 = [0.0,0.0]
    hard, error[0], chi2[0] = integrate(3, eps_tol, qgreal, iters, nstart, nend, alpha_vegas)
    if (subtr):
        sub, error[1], chi2[1] = integrate(2, eps_tol, qgsub, iters, nstart, nend, alpha_vegas)
    else:
        sub = 0.0
    intqg = (hard+sub)*cf/2.0
    error = np.sqrt(error[0]**2 + error[1]**2)*cf/2.0
    return intqg, error

def intqq(qqreal):
    error = 0.0
    chi2 = 0.0
    hard, error, chi2 = integrate(3, eps_tol, qqreal, iters, nstart, nend, alpha_vegas)

    intqq = hard*32/27.0
    error = error*32/27.0
    return intqq, error

if __name__ == "__main__":
    subtr = True
    ptcut = False
    rapcut = False
    ppcoll = True

    norderggh = 1
    minpt = 0.0
    maxpt = sqrtS*(1-z)/2.0
    minrap = 0.0
    maxrap = -np.log(z)/2.0

    sigmalo, errorlo, sigmanlo, errornlo = xs(norderggh=norderggh)
    print(f'sigmalo, errorlo, sigmanlo, errornlo: {sigmalo}, {errorlo},{sigmanlo},{errornlo}')
    

    """
    norderggh = 0
    ndim = 1
    sigmalo, errorlo = xs(norderggh=norderggh)
    print(f'sigmalo, errorlo: {sigmalo}, {errorlo}')
    """