import numpy as np
from pylooptools.ltools import *
from utils import *
from mpmath import *
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)
mp.pretty = True

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
    #print('Integral2: {Integral2}')
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
    Integral3 = complex(re*2.0/(m1-m2)/mq2, im/mq2)
    #print(f'Integral3, ma, mb, mq2: {Integral3}, {ma}, {mb}, {mq2}')
    return Integral3

# from pylooptools, not numerical stable.
def Integral4(m1,m2,m3,mb2):
    D0(0.0,0.0,0.0, m1+m2+m3, m1,m3 ,mb2,mb2,mb2,mb2)
    return D0(0.0,0.0,0.0, m1+m2+m3, m1,m3 ,mb2,mb2,mb2,mb2)