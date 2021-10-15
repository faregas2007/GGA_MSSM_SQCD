# this is part of Sushi
# written by N.T.T.Dat

import lhapdf
import numpy as np
from timeit import default_timer as timer
from pylooptools.ltools import *
import vegas

from slha_input import *
from reals import *
from lhapdfs import *

lhapdf.setVerbosity(0)

def integrate(ndim, eps, integrand, iters, nstart, nend, alpha_vegas):
    np.random.seed((1,2,3))
    integ = vegas.Integrator(ndim*[[eps, 1-eps]], sync_ran=True)

    # adapt vegas into integrand
    integ(integrand, nitn=iters, neval=nstart, alpha=alpha_vegas)
    result = integ(integrand, nitn=iters, neval=nend, alpha=alpha_vegas)

    return result, result.sdev, result.chi2

def trans(integrand, xx):
    """
    performs the integrals
    int_0^1 dv int_0^1 dx1 int_0^1 dx2 Heaviside(z/x1/x2)Heaviside(1-z/x1/x2) integrand(z/x1/x2, s,t,u,x1,x2)
    with z = mh**2/cme**2, cme=hadronic center of mass energy
    |d(x1,x2)/d(tau1, tau2)|=(1-z)**2*tau1*x2/w2**2
    """
    tau1 = xx[0]
    tau2 = xx[1]
    tau3 = xx[2]

    w1 = tau1*tau2*(1-z) + z
    w2 = tau1*(1-z) + z

    x1 = w1/w2
    x2 = z/w1
    jdet = (1-z)**2*tau1*x2/w2**2

    # for space space integration
    s = mh2/w2
    t = mh2*(1-1/w2)*tau3
    u = mh2*(1-1/w2)*(1-tau3)

    if(pesudorap):
        # pesudo-rapidity
        y = np.log(x1*t/x2/u)/2.0
    else:
        # rapidity
        y = np.log(x1*(mh2-t)/x2/(mh2-u))/2.0
    pt = np.sqrt(t*u/s)
    #print(f'pesudorap, y: {pesudorap}, {y}')
    trans = integrand(w2,s,t,u,x1,x2, pt,y)*jdet
    return trans

def transplus(integrand1, integrand2, xx):
    """
    To handle the plus distribution integration
    integrand1 is a Catani-Seymour subtraction term
    integrand2 contains only plus distribution from integrand1
    int_0^1 dv int_0^1 dx1 int_0^1 dx2 Heaviside_theta(z/x1/x2) Heaviside_theta(1-z/x1/x2) integrand1(z/x1/x2,s,t,u,x1,x2)

    mapping into unit hypercube
    x1->exp(y)*sqrt(z)
    x2->exp(-y)*sqrt(z)/w2
    with
    w2 = tau1*(1-z) + z
    y = ymin + (ymax-ymin)*tau2
    with the Jacobian determinant
    |d(x1,x2)/d(tau1,tau2)| = z/w2**2*(1-z)(ymax-ymin)
    # enforce minrap and maxrap for now, inclusive cross-section
    """
    minrap = 0.0
    maxrap = -np.log(z)/2.0
    tau1 = xx[0]
    tau2 = xx[1]

    jdet1 = 1-z
    jdet2 = maxrap - minrap
    y = minrap + jdet2*tau2
    x1 = np.exp(y)*np.sqrt(z)
    x2 = z/x1
    w2 = tau1*(1-z)+z
    transplus = -integrand2(w2,x1,x2)*z*jdet1*jdet2
    transplus = transplus*2.0 # add subtraction terms for t->0 and u->0
    
    if(x1 <= w2):
        # subtraction term for t->0
        transplus = transplus + integrand1(w2, x1/w2, x2)*z/w2**2*jdet1*jdet2

    if(x1>= z/w2):
        # subtraction term for u->0
        transplus = transplus + integrand1(w2, x2/w2, x1)*z/w2**2*jdet1*jdet2
    
    transplus = transplus*2.0 # include negative y
    return transplus

def transdelta(integrand, xx):
    """ 
    to deal wit delta dirac terms
    int_0^1 dx1 int_0^1 dx2 Heaviside_theta(z/x1/x2) Heaviside_theta(1 - z/x1/x2) pdf(x1)/x1 pdf(x2)/x2 * Dirac_delta(1 - z/x1/x2)
    =int_0^1 dx1 int_0^1 dx2 Heaviside_theta(z/x1/x2) Heaviside_theta(1 - z/x1/x2) pdf(x1)/x1 pdf(x2)/x2 * x2 * Dirac_delta(x2 - z/x1)
    x1 -> exp(y)*sqrt(z)
    and
    y -> ymin + (ymax - ymin)*tau
    with the Jacobian determinant
    |dx1/dy|*|dy/tau| = x1*(ymax-ymin)
    enforce maxrap and minrap for now
    """
    maxrap = -np.log(z)/2.0
    minrap = 0.0
    tau = xx[0]
    jdet1 = maxrap - minrap
    y = minrap + jdet1*tau
    x1 = np.exp(y)*np.sqrt(z)
    jdet2 = x1
    x2 = z/x1
    transdelta = x2*jdet1*jdet2*integrand(x1,x2)
    transdelta = transdelta * 2.0 # to include negative y
    return transdelta

# delta terms in gg channel
def deltaggint(x1,x2):
    return PDFgg(x1, x2, muFggh)/(x1*x2)

def deltagg(x1):
    return transdelta(deltaggint , [x1])

# Catani-Seymour subtraction (gg channel)
def intsubgg(x1, x2):
    return  transplus(subggint1, subggint2, [x1,x2])

def subggint1(x,x1,x2):
    return (ggterms(x) + dterms(x))/2.0 * PDFgg(x1, x2, muFggh)/(x1*x2)

def subggint2(x,x1,x2):
    return  (dterms(x)/2.0) * PDFgg(x1, x2, muFggh)/(x1*x2)

def ggterms(x):
    return (4*np.log(1-x) - 2*lfh)*(-2*x + x**2 - x**3) -2.0*np.log(x)*(1 - x + x**2)**2/(1.0 - x)

def dterms(x):
    dd0 = 1/(1.0 -x)
    dd1 = np.log(1.0 -x)*dd0
    return -2*lfh * dd0 + 4*dd1

def dtermsz():
    dd0 = np.log(1-z)/(1-z)
    dd1 = (np.log(1-z)/2.0)*dd0
    return -2*lfh * dd0  + 4 *dd1

def subggt(s,t,u):
    return  (s**2 + s*(t+u) + (t+u)**2)**2/(s*t*(t+ u)*(s+t+u))

def subggu(s,t,u):
    return (s**2 + s*(t+ u) + (t+u)**2)**2/(s*u*(t+u)*(s+t+u))

# real correction gg channel (finite terms)
# set 
def realggint(x,s,t,u,x1,x2,pt,y):
    if (t >= 0.0):
        t = -s*10**-6
    if (u >= 0.0):
        u = -s*10**-16

    tmp = AMPgg(s,t,u)*jet(pt,y)
    if (subtr):
        if ((-s/t > 10**5) or (-s/u > 10**5)):
            return 0.0
        tmp = tmp - subggt(s,t,u)*jetsub(x1*x, x2) - subggu(s,t,u)*jetsub(x1, x2*x)

    return x*(1-x)*tmp*PDFgg(x1, x2, muFggh)/(x1*x2)

def realgg(x, x1, x2):
    return trans(realggint, [x1, x2, x])


# qg channel
# Catani-Seymour subtraction term (qg and gq channel)
# 2 dimensional integration
def intsubqg(x1, x2):
    return  transplus(subqgint1, subqgint2, [x1,x2])

def subqgint1(x, x1, x2):
    # integrated Catani Seymour subtraction term
    return  qgterms(x)*PDFqg(x1, x2, muFggh)/(x1*x2)

def subqgint2(x, x1, x2):
    return 0.0

def qgterms(x):
    return ((1 + (1 - x)**2)*(2*np.log(1 - x)-np.log(x) - lfh) + x**2)

# qg main terms
def realqgint(x,s,t,u,x1,x2, pt, y):
    tmpqg = AMPqg(s,t,u)*jet(pt,y)
    tmpgq = AMPqg(s,u,t)*jet(pt,y)

    if(subtr):
        tmpqg = tmpqg - subqg(s,t,u)*jetsub(x1*x, x2)
        tmpgq = tmpgq - subqg(s,u,t)*jetsub(x1, x2*x)
    return x*(1-x)*(tmpqg*PDFqg(x1, x2, muFggh)/(x1*x2) + tmpgq*PDFqg(x2, x1, muFggh)/(x2*x1))

def subqg(s, t, u):
    return -(s**2 + (t+u)**2)/(t*(s+t+u))

# total xs
def realqg(x, x1,x2):
    return trans(realqgint, [x1, x2, x])

# qq channel
def AMPqqjet(s,t,u):
    return (t**2 + u**2 )/(t+u)**2*3.0/2.0*AMPqq(s,t+u)

def realqqint(x,s,t,u,x1,x2, pt,y):
    return (1-x)**3*AMPqqjet(s,t,u)*jet(pt,y)*PDFqq(x1, x2, muFggh)/(x1*x2)

def realqq(x,x1,x2):
    return trans(realqqint,[x1,x2,x])

# jet functions
def jetsub(x1,x2):
    if(pesudorap):
        jetsub = jet(0.0, 1.0)
    else:
        jetsub = jet(0.0, np.log(x1/x2)/2.0)

    return jetsub

def jet(pt, y):
    jet = 1.0
    if (rapcut and (np.abs(y) > maxrap or np.abs(y)<minrap)): 
        jet = 0.0
    if (ptcut and (np.abs(pt) > maxptc or np.abs(pt) <minptc)):
        jet = 0.0
    return jet



########## integrand for integration for gg channel
def ggreal(x):
    return realgg(x[0], x[1], x[2])

def ggsub(x):
    return intsubgg(x[0], x[1])

def ggdelta(x):
    return deltagg(x[0])

########## integrand for integration qg channel
def qgreal(x):
    return realqg(x[0], x[1], x[2])

def qgsub(x):
    return intsubqg(x[0], x[1])

########## integrand for integration for qq channel
def qqreal(x):
    return realqq(x[0], x[1], x[2])


