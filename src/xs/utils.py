# utils.py

import numpy as np
from xs.integrals import *

def C1SMA(s, t, u, m2):
    stu4 = Integral4(s,t,u,m2)
    sut4 = Integral4(s,u,t,m2)
    tsu4 = Integral4(t,s,u,m2)
    s3 = Integral3(s,s+t+u,m2)
    t3 = Integral3(t,s+t+u,m2)
    u3 = Integral3(u,s+t+u,m2)

    return -(s+t+u)*(2*((t+u)*s3 + (s+u)*t3 + (s+t)*u3) + s*u*stu4 + s*t*sut4 + t*u*tsu4)

def C1SMA_batch(s,t,u,m2):
    stu4 = Integral4_batch(s,t,u,m2)
    sut4 = Integral4_batch(s,u,t,m2)
    tsu4 = Integral4_batch(t,s,u,m2)
    s3 = Integral3_batch(s,s+t+u,m2)
    t3 = Integral3_batch(t,s+t+u,m2)
    u3 = Integral3_batch(u,s+t+u,m2)

    return -(s+t+u)*(2*((t+u)*s3 + (s+u)*t3 + (s+t)*u3) + s*u*stu4 + s*t*sut4 + t*u*tsu4)

def C2SMAcalc(s,t,u,m2):
    stu4 = Integral4(s,t,u,m2)
    sut4 = Integral4(s,u,t,m2)
    tsu4 = Integral4(t,s,u,m2)
    s3 = Integral3(s,s+t+u,m2)
    t3 = Integral3(t,s+t+u,m2)
    u3 = Integral3(u,s+t+u,m2)

    return -s * (2*((t+u)*s3 + (s-u)*t3 + (s-t)*u3) + s*u*stu4 + s*t*sut4 - t*u*tsu4)

def C2SMAcalc_batch(s,t,u,m2):
    stu4 = Integral4_batch(s,t,u,m2)
    sut4 = Integral4_batch(s,u,t,m2)
    tsu4 = Integral4_batch(t,s,u,m2)
    s3 = Integral3_batch(s,s+t+u,m2)
    t3 = Integral3_batch(t,s+t+u,m2)
    u3 = Integral3_batch(u,s+t+u,m2)

    return -s * (2*((t+u)*s3 + (s-u)*t3 + (s-t)*u3) + s*u*stu4 + s*t*sut4 - t*u*tsu4)


def C2SMA(s, t, u, m2):
    c2 = C2SMAcalc(s,t,u,m2)
    c3 = C2SMAcalc(t,s, u, m2)
    c4 = C2SMAcalc(u,s,t,m2)

    return c2, c3, c4

def C2SMA_batch(s,t,u,m2):
    c2 = C2SMAcalc_batch(s,t,u,m2)
    c3 = C2SMAcalc_batch(t,s, u, m2)
    c4 = C2SMAcalc_batch(u,s,t,m2)

    return c2, c3, c4

def ALOSMA(mh2, m2):
    return -3.0*Integral3(0.0, mh2, m2)

def ALOSMA_batch(mh2, m2):
    if(np.array(mh2).size == 1): 
        zeros = np.zeros((np.array(mh2).size), )
    else:
        zeros = np.zeros(mh2.size, )
    return -3.0*Integral3_batch(zeros, mh2, m2)

def ASMA(s, tpu, m2):
    return 3.0*Integral3(s, s+tpu, m2)

def ASMA_batch(s, tpu, m2):
    return 3.0*Integral3_batch(s, s+tpu, m2)