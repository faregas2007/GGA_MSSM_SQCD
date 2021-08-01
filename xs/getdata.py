import pandas as pd
import numpy as np
import argparse
import os

path = '/xs'
os.chdir(path)
pd.set_option('display.precision', 12)
params = argparse.Namespace

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
