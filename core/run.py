#!/usr/bin/env python
from __future__ import division
import vegas
import numpy as np
import gvar as gv
import os, sys
import pyslha
import pycuba
from auxiliary import *
from importlib import import_module
from timeit import default_timer as timer
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

islha_obj = pyslha.read(sys.argv[1])
iters = islha_obj.blocks['VEGAS'][1]
startdiv = islha_obj.blocks['VEGAS'][2]
enddiv = islha_obj.blocks['VEGAS'][3]
startfin =islha_obj.blocks['VEGAS'][4]
endfin = islha_obj.blocks['VEGAS'][5]
alphafix = islha_obj.blocks['VEGAS'][6]
integrator = islha_obj.blocks['VEGAS'][7]
eps = islha_obj.blocks['IBPPARA'][2]
corrf = islha_obj.blocks["REAL"][7]

mA = islha_obj.blocks['MASS'][1]
beta = islha_obj.blocks['REAL'][3]

squark = islha_obj.blocks['INTEGER'][1]
lam = islha_obj.blocks['IBPPARA'][1]
eps = islha_obj.blocks['IBPPARA'][2]

tanb = np.round(np.tan(beta), 1)
if(np.floor(tanb)-tanb == 0.0):
	tanb = int(tanb)

# get module of function file
function_dir = sys.argv[2]
index_in = int(sys.argv[3])
#dia_index = sys.argv[3]

dia_index = [2, 3, 5, 7, 9, 11, 12, 13, 14, 15]
sys.path.insert(0, function_dir)

vegasbatch = False
	
def num_integrate(index_in):
	if (index_in == 9):
		module = 'ampibp_total_'+'CT'
	elif (index_in == 3 and vegasbatch == True):
		module = 'ampibp_total_dia7'
	else:
		module = 'ampibp_total_dia'+str(dia_index[index_in])
	fcn_module = import_module(module)
	
	ibpptotal = getattr(fcn_module, 'ibpptotal')
	ibpftotal = getattr(fcn_module, 'ibpftotal')

	np.random.seed((1,2,3))
	fname=['real', 'imag']
	integ = vegas.Integrator(5*[[eps, 1. - eps]], nhcube_batch=2000, sync_ran=True)
	
	pole = []
	poleerr = []
	fin = []
	finerr = []

	start = timer()
	if (index_in == 9):
		for i in range(2):
			pole.append(ibpptotal(i))
			poleerr.append(ibpptotal(i))
			fin.append(ibpftotal(i))
			finerr.append(ibpftotal(i))
	# speed up for diagram 7
	# Cuhre integrator is the best one.
	elif (index_in == 3 and vegasbatch == False):
		NDIM = 2
		NCOMP = 1
		KEY = 0
		VERBOSE = 0
		module = 'ampibp_total_dia7_cuba'
		fcn_module = import_module(module)
		
		freal = getattr(fcn_module, 'freal')
		fimag = getattr(fcn_module, 'fimag')
		preal = getattr(fcn_module, 'preal')
		pimag = getattr(fcn_module, 'pimag')
		
		for i in [preal, pimag]:
			temp = Integrator(integrator, i, NDIM, KEY, startfin, endfin, VERBOSE)
			result = float(temp['result'])
			error = float(temp['error'])
			pole.append(result)
			poleerr.append(error)
			
		for i in [freal, fimag]:
			temp = Integrator(integrator, i, NDIM, KEY, startdiv, enddiv, VERBOSE)
			result = float(temp['result'])
			error = float(temp['error'])
			fin.append(result)
			finerr.append(error)
	else:
		for i in range(2):
			integ(ibpptotal(i), nitn=iters, neval=startdiv, alpha=alphafix)
			result = integ(ibpptotal(i), nitn=iters, neval=enddiv, alpha=alphafix)
			pole.append(result.mean)
			poleerr.append(result.sdev)

		for i in range(2):
			integ(ibpftotal(i), nitn=iters, neval=startfin, alpha=alphafix)
			result = integ(ibpftotal(i), nitn=iters, neval=endfin, alpha=alphafix)
			
			if(index_in == 0 or index_in == 6 or index_in == 8):
				fin.append(corrf*result.mean)
				finerr.append(corrf*result.sdev)

			fin.append(result.mean)
			finerr.append(result.sdev)
	end = timer()
	return {'ReFinAvg':fin[0], 'ReFinErr':finerr[0] ,'ImDivAvg':fin[1], 'ImDivErr':finerr[1], 'ReDivAvg':pole[0], 'ReDivErr':poleerr[0],'ImDivAvg':pole[1], 'ImDivErr':poleerr[1], 'Time':end-start}


def csvwrite(name, index_in):
	import csv
	fields = ['diagrams', 'squark', 'who', 'ReFinAvg', 'ReFinErr', 'ReDivAvg', 'ReDivErr', 'ImFinAvg', 'ImFinErr', 'ImDivAvg', 'ImDivErr', 'lam', 'eps', 'TIME']
	temp =num_integrate(index_in)
	diagram = dia_index[index_in]
	DDS = float(0.0)
	who = 'KIT'
	ReFinAvg = temp['ReFinAvg']
	ReFinErr = temp['ReFinErr']
	ImFinAvg = temp['ImDivAvg']
	ImFinErr = temp['ImFinErr']
	ReDivAvg = temp['ReDivAvg']
	ReDivErr = temp['ReDivErr']
	ImDivAvg = temp['ImDivAvg']
	ImDivErr = temp['ImDivErr']
	TIME = temp['TIME']
	
	if (squark == 1):
		squarks = 'stop'
	elif (squark == 2):
		squarks = 'sbot'
		
	if (index_in == 9):
		filename = 'CT'+sys.argv[1]+'.csv'
		rows = [[str(diagrams), str(squark), who, str(ReFinAvg), str(ReFinErr), str(ReDivAvg), str(ReDivErr), str(ImFinAvg), str(ImFinErr), str(ImDivAvg), str(ImDivErr), str(lam), str(eps), str(TIME)]]
	else:
		filename = str(diagrams)+sys.argv[1]+'.csv'
		rows = [[str(diagrams), str(squark), who, str(ReFinAvg), str(ReFinErr), str(ReDivAvg), str(ReDivErr), str(ImFinAvg), str(ImFinErr), str(ImDivAvg), str(ImDivErr), str(lam), str(eps), str(TIME)]]
	with open(filename, 'w') as csvfile:
		csvwriter = csv.writer(csvfile)
		csvwriter.writerow(fields)
		csvwriter.writerows(rows)
	csvfile.close

if __name__ == '__main__':
	if True:
		name = sys.argv[1]
		csvwrite(name, index_in)
	else:
		import hotshot, hotshot.stats
		prof = hotshot.Profile('vegas.prof')
		prof.runcall(main)
		prof.close()
		stats = hotshot.stats.load('vegas.prof')
		stats.strip_dirs()
		stats.sort_status('time', 'calls')
		stats.print_stats(40)
