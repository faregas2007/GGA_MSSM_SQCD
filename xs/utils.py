import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from mpmath import *
import os

mp.pretty = True
mp.dps = 12

base_dir = Path(__file__).parent.parent.absolute()
path = base_dir
os.chdir(path)
pd.set_option('display.precision', 12)
params = argparse.Namespace

def polylog2(x):
    return polylog(2,x)

def polylog3(x):
    return polylog(3,x)

def getdata(filename):
    df = pd.read_csv(filename)  
    df2 = df[df['squark']==1].reset_index()
    df3 = df[df['squark']==2].reset_index() 
    params.stop = df2['squark']
    params.sbot = df3['squark'] 
    params.tanb = df2['tanb']
    params.mu = df2['mu']
    params.mA = df2['mA']

    params.timet = df2['TIME']
    params.lamt = df2['lam']
    params.epst = df2['eps']
    params.reCt = df2['ReFinAvg']
    params.reCterr = df2['ReFinErr']
    params.imCt = df2['ImFinAvg']
    params.imCterr = df2['ImFinErr']
    params.mst1 = df2['msq1']
    params.mst2 = df2['msq2']
    params.At = df2['Aq']
    params.alphaqst = df2['alphaqs']

    params.timeb = df3['TIME']
    params.lamb = df3['lam']
    params.epsb = df3['eps']
    params.reCb = df3['ReFinAvg']
    params.reCberr = df3['ReFinErr']
    params.imCb = df3['ImFinAvg']
    params.imCberr = df3['ImFinErr']
    params.msb1 = df3['msq1']
    params.msb2 = df3['msq2']
    params.Ab = df3['Aq']
    params.alphaqsb = df3['alphaqs']

    return params

def Ifunc(a, b, c):
    return (a*b*np.log(a/b) + b*c*np.log(b/c) + a*c*np.log(c/a))/((a-b)*(b-c)*(a-c))

