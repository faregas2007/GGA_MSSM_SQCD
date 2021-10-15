# This is part of Sushi
# written in python by N.T.T.Dat
import numpy as np
from integrals import *
from slha_input import *

def C1SMA(s, t, u, m2):
    stu4 = Integral4(s,t,u,m2)
    sut4 = Integral4(s,u,t,m2)
    tsu4 = Integral4(t,s,u,m2)
    s3 = Integral3(s,s+t+u,m2)
    t3 = Integral3(t,s+t+u,m2)
    u3 = Integral3(u,s+t+u,m2)

    return -(s+t+u)*(2*((t+u)*s3 + (s+u)*t3 + (s+t)*u3) + s*u*stu4 + s*t*sut4 + t*u*tsu4)

def C2SMAcalc(s,t,u,m2):
    stu4 = Integral4(s,t,u,m2)
    sut4 = Integral4(s,u,t,m2)
    tsu4 = Integral4(t,s,u,m2)
    s3 = Integral3(s,s+t+u,m2)
    t3 = Integral3(t,s+t+u,m2)
    u3 = Integral3(u,s+t+u,m2)

    return -s * (2*((t+u)*s3 + (s-u)*t3 + (s-t)*u3) + s*u*stu4 + s*t*sut4 - t*u*tsu4)

def C2SMA(s, t, u, m2):
    stu4 = Integral4(s,t,u,m2)
    sut4 = Integral4(s,u,t,m2)
    tsu4 = Integral4(t,s,u,m2)
    s3 = Integral3(s,s+t+u,m2)
    t3 = Integral3(t,s+t+u,m2)
    u3 = Integral3(u,s+t+u,m2)

    lstu4 = stu4
    lsut4 = sut4
    ltsu4 = tsu4
    ls3 = s3
    lt3 = t3
    lu3 = u3

    c2 = C2SMAcalc(s,t,u,m2)

    stu4 = ltsu4
    tsu4 = lstu4
    s3 = lt3
    t3 = ls3
    c3 = C2SMAcalc(t,s, u, m2)

    sut4 = lstu4
    tsu4 = lsut4
    s3 = lu3
    u3 = lt3
    c4 = C2SMAcalc(u,s,t,m2)

    return c2, c3, c4

# squared amplitude for gg->gh
def AMPgg(s,t,u):
    # only focus on full mass dependence
    return AMPggpure(s,t,u)/AMPLO(s+t+u)

def AMPggpure(s,t,u):
    # squared amplutude for gg->gh
    charm1, charm2, charm3, charm4 = complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0)
    bottom1, bottom2, bottom3, bottom4 = complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0)
    top1, top2, top3, top4 = complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0), complex(0.0,0.0)
    c1, c2, c3, c4 = 0.0,0.0,0.0,0.0
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

# squared amplitude for qg->qh
def AMPqg(s,t,u):
    # squared amplitude for qg->qh
    # normalized to the LO amplitude
    return AMPqgpure(s,t,u)/AMPLO(s+t+u)

def AMPqgpure(s,t,u):
    # squared amplitude for qg->qh
    # without normalization
    A = 0.0
    if(np.abs(gc) != 0.0):
        A = A + mc2*gc*ASMA(t, s+u, mc2)
    
    if(np.abs(gb) != 0.0):
        A = A + mb2*gb*ASMA(t,s+u, mb2)
    
    if(np.abs(gt) != 0.0):
        A = A + mt2*gt*ASMA(t, s+u, mt2)
    
    return  -(s**2 + u**2)/t*np.abs(A)**2/(s+t+u)


# squared amplitude for qq->qh
def AMPqq(s, tpu):
    # squared amplitude for qq->gh
    # normalized to the LO amplitude

    return AMPqqpure(s,tpu)/AMPLO(s+tpu)

def AMPqqpure(s, tpu):
    # squared amplitude for qq->qh
    # without normalization

    A = 0.0
    if(np.abs(gc) != 0.0):
        A = A + mc2*gc*ASMA(s, tpu, mc2)
    
    if(np.abs(gb) != 0.0): 
        A = A + mb2*gb*ASMA(s, tpu, mb2)
    
    if(np.abs(gt) != 0.0):
        A = A + mt2*gt*ASMA(s, tpu, mt2)

    return np.abs(A)**2

def AMPLO(mh2):
    AMP = 0.0
    if(np.abs(gc) != 0.0):
        AMP = AMP + mc2*gc*ALOSMA(mh2, mc2)
    if(np.abs(gb) != 0.0):
        AMP = AMP + mb2*gb*ALOSMA(mh2, mb2)
    if(np.abs(gt) != 0.0):
        AMP = AMP + mt2*gt*ALOSMA(mh2, mt2)
    
    AMPLO = np.abs(AMP)**2
    #return np.abs(AMP)**2
    return AMPLO

def ALOSMA(mh2, m2):
    return -3*Integral3(0.0, mh2, m2)

def ASMA(s, tpu, m2):
    return  3.0 *Integral3(s, s+tpu, m2)


"""
# sushi result 0.35522100952687402        
s = 151911348.39435318       
t= -123602374.64065668       
u= -28268973.753696509 
print(AMPgg(s,t,u))
"""