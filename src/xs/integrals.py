# inptegrals.py

import numpy as np
import pylooptools
from pylooptools.ltools import *
from utils import *
from mpmath import *
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)
mp.pretty = True
mp.dps = 12

# with QCD loop --> clear
"""
def Integral2(m1, mb2):
    B0(m1,mb2, mb2)
    return B0(m1, mb2, mb2)
    
def Integral3(m1,m2,mb2):
    C0(m1,0.0, m2, mb2, mb2, mb2)
    return C0(m1, 0.0, m2, mb2, mb2, mb2)
"""

# without QCD loops --> not clear

# sushi used these integrals, without QCD loop and mq2 != 0
def Integral2(ma, mq2):
    m1 = ma/mq2
    if(m1 == 0.0):
        re = -2.0
        im = 0.0
    elif(m1 < 0.0):
        re = -2.0*np.sqrt(1.0 - 4.0/m1)*np.log((np.sqrt(-m1) + np.sqrt(-m1 + 4.0))/2.0)
        im = 0.0
    elif(m1 <= 4.0):
        re = -2.0*np.sqrt(4.0/m1 - 1.0)*np.arcsin(np.sqrt(m1)/2.0)
        im = 0.0
    elif(m1 > 4.0):
        re = -2.0*np.sqrt(1-4.0/m1)*(np.log((np.sqrt(m1) + np.sqrt(m1-4.0))/(2.0)))
        im = 2.0*np.sqrt(1.0 - 4.0/m1)*np.pi/2.0
    else:
        print('Sushi error in Integral2')
        
    Integral2 = complex(re, im)
    return Integral2

def Integral3(ma, mb, mq2):
    m1 = ma/mq2
    m2 = mb/mq2
    if (m1 == 0.0):
        re = 0.0
        im = 0.0
    elif (m1 < 0.0):
        re = (np.log(2.0/(2.0-m1 + np.sqrt(m1*(-4+m1))))**2)/4.0
        im = 0.0
    elif (m1 <= 4.0):
        re = - np.arcsin(np.sqrt(m1/4.0))**2
        im = 0.0
    elif (m1 > 4.0):
        re = (np.log((-2.0 + m1+np.sqrt(m1*(-4.0+m1)))/(2.0))**2 - np.pi**2)/4.0
        im = -np.pi*np.log((-2.0 + m1 + np.sqrt(m1*(-4.0 + m1)))/2.0)/(m1-m2)
    
    if (m2 == 0.0):
        re = 0.0
        im = 0.0
    elif (m2 < 0.0):
        re = re - (np.log(2.0/(2.0 - m2 + np.sqrt(m2*(-4+m2))))**2)/4.0
    elif(m2 <= 4.0):
        re = re + np.arcsin(np.sqrt(m2/4.0))**2
    elif(m2 > 4.0):
        re= re - (np.log((-2.0 + m2 + np.sqrt(m2*(-4.0+m2)))/2.0)**2 - np.pi**2)/4.0
        im = im + np.pi*np.log((-2.0 + m2 + np.sqrt(m2*(-4.0+m2)))/(2.0))/(m1-m2)
    else:
        print("Sushi error in Integral3")
    Integral3=(complex(re*2.0/(m1-m2)/mq2, im/mq2))
    #print(f'Integral3, ma, mb, mq2: {Integral3}, {ma}, {mb}, {mq2}')
    return Integral3

def Integral3_batch(ma, mb, mq2):
    # intialize
    ma = np.array(ma)
    mb = np.array(mb)
    mq2 = np.array(mq2)
    
    if(ma.shape == ()):
        ma = np.array([ma])
    if(mb.shape == ()):
        mb = np.array([mb])
    if(mq2.shape == ()):
        mq2 == np.array([mq2])

    Integral3 = []
    re = np.zeros(len(mb),)
    im = np.zeros(len(mb),)
    M1 = np.array(ma/mq2)
    M2 = np.array(mb/mq2)
    
    # slow implementation batch mode.
    for idx, (m1, m2) in enumerate(zip(M1, M2)):

        if (m1 == 0.0):
            re[idx] = 0.0
            im[idx] = 0.0
        elif (m1 < 0.0):
            re[idx] = (np.log(2.0/(2.0-m1 + np.sqrt(m1*(-4+m1))))**2)/4.0
            im[idx] = 0.0
        elif (m1 <= 4.0):
            re[idx] = - np.arcsin(np.sqrt(m1/4.0))**2
            im[idx] = 0.0
        elif (m1 > 4.0):
            re[idx] = (np.log((-2.0 + m1+np.sqrt(m1*(-4.0+m1)))/(2.0))**2 - np.pi**2)/4.0
            im[idx] = -np.pi*np.log((-2.0 + m1 + np.sqrt(m1*(-4.0 + m1)))/2.0)/(m1-m2)

        if (m2 == 0.0):
            re[idx] = 0.0
            im[idx] = 0.0
        elif (m2 < 0.0):
            re[idx] -= (np.log(2.0/(2.0 - m2 + np.sqrt(m2*(-4+m2))))**2)/4.0
        elif(m2 <= 4.0):
            re[idx] += np.arcsin(np.sqrt(m2/4.0))**2
        elif(m2 > 4.0):
            re[idx] -= (np.log((-2.0 + m2 + np.sqrt(m2*(-4.0+m2)))/2.0)**2 - np.pi**2)/4.0
            im[idx] += np.pi*np.log((-2.0 + m2 + np.sqrt(m2*(-4.0+m2)))/(2.0))/(m1-m2)
        else:
            print("Sushi error in Integral3")
    for idx in range(len(re)):
        Integral3.append(complex(re[idx]*2.0/(m1-m2)/mq2, im[idx]/mq2))
    return np.array(Integral3)

# from pylooptools, not numerical stable.
"""
def Integral4(m1, m2,m3,mb2):
    return D0(0.0,0.0,0.0,m1+m2+m3, m1,m3,mb2,mb2,mb2,mb2)
    
def Integral4_batch(m1,m2,m3,mb2):
    # slow implementation for batch mode. 
    Integral4 = []
    #mb2 = np.full(len(m2), mb2)

    if(np.array(m1).size ==1):
        Integral4.append(D0(0.0,0.0,0.0,m1+m2+m3, m1, m3, mb2, mb2, mb2, mb2))
        return np.array(Itegral4)
    else:
        
        for m11, m22, m33 in zip(m1,m2,m3):
            Integral4.append(D0(0.0,0.0,0.0, m11+m22+m33, m11, m33, mb2, mb2, mb2, mb2))
    return np.array(Integral4)
"""
def Integral4(m1,m2,m3, mb2):
    #m1 = ma/mq2
    #m2 = mb/mq2
    #m3 = mc/mq2

    return D0(0.0,0.0,0.0,m1+m2+m3,m1,m3,mb2,mb2,mb2,mb2)

def Integral4_batch(m1,m2,m3,mb2):
    #m1 = ma/mq2
    #m2 = mb/mq2
    #m3 = mc/mq2

    Integral4 = []
    if(np.array(m1).size == 1):
        Integral4.append(D0(0.0,0.0,0.0,m1+m2+m3,m1,m3,mb2,mb2,mb2,mb2))
        return np.array(Integral4)
    else:
        for m11, m22, m33 in zip(m1,m2,m3):
            Integral4.append(D0(0.0,0.0,0.0,m11+m22+m33,m11,m33,mb2,mb2,mb2,mb2))
        return np.array(Integral4)
# fixed Integral4 numerical unstability problem
"""
def Integral4_sushi(ma, mb, mc, mq2):
    m1 = ma/mq2
    m2 = mb/mq2
    m3 = mc/mq2
    pn = [m1, m1+m2+m3, m3]
    p = [0.0,0.0,0.0]
    tmp = 1 + m1*m3/m2/2.0
    if (tmp**2 == 1.0):
        Integral4 = 1/6.0/mq2**2
        return Integral4
    
    x = [0.0,0.0]
    x[1] = -tmp - np.sqrt(tmp**2-1)
    x[0] = 1/x[1]

    Integral4 = complex(2.0*(2.0*polylog(2, 1+x[1])) + np.log(-x[1])**2/2.0, 0.0)
    for i in range(1,3):
        wrz = np.sqrt(complex((1-pn[i]/2.0)**2-1, 0.0))
        if (pn[i] == 0.0):
            p[i] = 1.0
        elif(pn[i] > 0.0):
            p[i] = 1/(1-pn[i]/2.0 - wrz)
        else:
            p[i] = 1 - pn[i]/2.0 + wrz
        
        for k in range(1,2):
            if(pn[i] < 0.0):
                arg = np.real(p[i]*x[k])
                Integral4 = Integral4 + (-1)**(i+k)*(-2*polylog(2, 1+1/arg) 
                        - np.log(-arg)**2/2)
            elif(pn[i] >= 4.0):
                arg = np.real(p[i]*x[k])
                Integral4 = Integral4 + (-1)**(i+k)*(polylog(2, -1/arg) - zeta2 
                        - np.log(1+1/arg)*np.log(arg) 
                        -(np.log(arg) - complex(0.0, np.pi))**2/4.0
                        -complex(0.0, np.pi)*np.log(1+arg))*2.0
            else:
                Integral4 = Integral4 + (-1)**(i+k)*(polylog(2,1+p[i]*x[k]) + polylog(2,1+x[k]/p[i]))

    Integral4 = Integral4/np.sqrt(tmp**2-1)/m2/mq2**2/2.0
    #print(f'Integral4, ma, mb, mc, mq2 : {Integral4}, {ma}, {mb}, {mc}, {mq2}')
    return Integral4
"""