#!/usr/bin/env python
from __future__ import division
import vegas
import numpy as np
import gvar as gv
import os, sys
import pyslha
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
eps = islha_obj.blocks['IBPPARA'][2]
corrf = islha_obj.blocks["REAL"][7]

mA = islha_obj.blocks['MASS'][1]
beta = islha_obj.blocks['REAL'][3]

tanb = np.round(np.tan(beta), 1)
if(np.floor(tanb)-tanb == 0.0):
	tanb = int(tanb)

# get module of function file
function_dir = sys.argv[2]
index_in = int(sys.argv[3])
#dia_index = sys.argv[3]

dia_index = [2, 3, 5, 7, 9, 11, 12, 13, 14, 15]
sys.path.insert(0, function_dir)

	
def num_integrate(index_in):
	if (index_in == 9):
		fcn = 'ampibp_total_'+'CT'+'.py'
	else:
		fcn = 'ampibp_total_dia'+str(dia_index[index_in])+'.py'
	module = fcn.strip('.py')
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

	if (index_in == 9):
		for i in range(2):
			pole.append(ibpptotal(i))
			poleerr.append(ibpptotal(i))
			fin.append(ibpftotal(i))
			finerr.append(ibpftotal(i))
	else:
		for i in range(2):
			start = timer()
			integ(ibpptotal(i), nitn=iters, neval=startdiv, alpha=alphafix)
			result = integ(ibpptotal(i), nitn=iters, neval=enddiv, alpha=alphafix)
			end = timer()
			pole.append(result.mean)
			poleerr.append(result.sdev)

		for i in range(2):
			start = timer()
			integ(ibpftotal(i), nitn=iters, neval=startfin, alpha=alphafix)
			result = integ(ibpftotal(i), nitn=iters, neval=endfin, alpha=alphafix)
			end = timer()
			
			if(index_in == 0 or index_in == 6 or index_in == 8):
				fin.append(corrf*result.mean)
				finerr.append(corrf*result.sdev)

			fin.append(result.mean)
			finerr.append(result.sdev)
	
	return {'ReDivAvg':fin[0], 'ReDivErr':finerr[0] ,'ImDivAvg':fin[1], 'ImDivErr':finerr[1], 'ReDivAvg':pole[0], 'ReDivErr':poleerr[0],'ImDivAvg':pole[1], 'ImDivErr':poleerr[1]}


def csvwrite(name, index_in):
	import csv
	temp =main(index_in)
	fields = ['squark', 'tanb', 'mA', 'diagrams', 'DDS', 'who', 'ReFinAvg', 'ReFinErr', 'ReDivAvg', 'ReDivErr', 'ImFinAvg', 'ImFinErr', 'Lambda', 'Eps']
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

	filename = 'deltatb'+name+'.csv'
	if(squark == 1):
		squarks = 'stop'
	else(squark == 2):
		squarks = 'sbot'
		
	if (index_in == 9):
		filename = 'CT'+sys.argv[1]+'.csv'
		rows = [[squarks, str(tanb), str(mA), who, str(DDS), who, str(ReFinAvg), str(ReFinErr), str(ReDivAvg), str(ReDivErr), str(ImFinAvg), str(ImFinErr), str(ImDivAvg), str(ImDivErr), str(Lambda), str(eps)]]
	else:
		filename = str(diagrams)+sys.argv[1]+'.csv'
		rows = [[squarks, str(tanb), str(mA), who, str(DDS), who, str(ReFinAvg), str(ReFinErr), str(ReDivAvg), str(ReDivErr), str(ImFinAvg), str(ImFinErr), str(ImDivAvg), str(ImDivErr), str(Lambda), str(eps)]]
	with open(filename, 'w') as csvfile:
		csvwriter = csv.writer(csvfile)
		csvwriter.writerow(fields)
		csvwriter.writerows(rows)
	csvfile.close

if __name__ == '__main__':
	if True:
		name = sys.argv[1]
		num_integrate(index_in)
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
