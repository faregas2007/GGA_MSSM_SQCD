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

def ALOSMA(mh2, m2):
    return -3*Integral3(0.0, mh2, m2)

def ASMA(s, tpu, m2):
    return 3.0*Integral3(s, s+tpu, m2)