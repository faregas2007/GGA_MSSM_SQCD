from __future__ import division
import vegas
import numpy as np
import gvar as gv
import os, sys
import csv
import pyslha
from auxiliary import *
from importlib import import_module
from timeit import default_timer as timer
from argparse import Namespace
from utils import *
from config import *
from slha_input import *
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

param_fp = Path(config_dir, "params.json")
params = Namespace(**load_dict(param_fp))

vegasbatch = params.vegasbatch
function_dir = amp_dir
index_in = params.index_in

dia_index = [2, 3, 5, 7, 9, 11, 12, 13, 14, 15]
sys.path.insert(0, str(function_dir))
sys.path.insert(0, str(input_dir))
def num_integrate(index_in=index_in, lam=lam):
	np.random.seed((1,2,3))
	if (index_in == 9):
		module = 'ampibp_total_diaCT'
	elif (index_in == 3 and vegasbatch == True):
		module = 'ampibp_total_dia7'
	else:
		module = 'ampibp_total_dia' + str(dia_index[index_in])

	fcn_module = import_module(module)
	ibpptotal=getattr(fcn_module, 'ibpptotal')
	ibpftotal=getattr(fcn_module, 'ibpftotal')
	integ = vegas.Integrator(5*[[eps, 1. - eps]], nhcube_batch=2000, sync_ran=True)
	pole = []
	poleerr = []
	fin = []
	finerr = []

	start = timer()
	if (index_in == 9):
		for i in range(2):
			pole.append(ibpptotal(i)[0])
			poleerr.append(0.0)
			fin.append(ibpftotal(i)[0])
			finerr.append(0.0)

	elif (index_in == 3 and vegasbatch == False):
		lam1 = lam
		NDIM = 2
		NCOMP =1
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
			integ(ibpptotal(i, lam), nitn=iters, neval=startdiv, alpha=alphafix)
			result = integ(ibpptotal(i,lam), nitn=iters, neval=enddiv, alpha=alphafix)

			if (index_in == 0 or index_in == 6 or index_in == 8):
				pole.append(corrf*result.mean)
				poleerr.append(corrf*result.sdev)
			else:
				pole.append(result.mean)
				poleerr.append(result.sdev)

		for i in range(2):
			integ(ibpftotal(i, lam), nitn=iters, neval=startfin, alpha=alphafix)
			result = integ(ibpftotal(i, lam), nitn=iters, neval=endfin, alpha=alphafix)
			if (index_in == 0 or index_in == 6 or index_in == 8):
				fin.append(corrf*result.mean)
				finerr.append(corrf*result.sdev)
			else:
				fin.append(result.mean)
				finerr.append(result.sdev)
	end = timer()
	print(fin[0], finerr[0], fin[1], finerr[1], lam)
	return {'ReFinAvg': fin[0], 'ReFinErr': finerr[0], 'ImFinAvg':fin[1], 'ImFinErr':finerr[1], 'ReDivAvg':pole[0], 'ReDivErr':poleerr[0], 'ImDivAvg':pole[1], 'ImDivErr':poleerr[1], 'TIME':end - start}

def csvwrite(working, name, index_in, lam1, richardson):
	temp = num_integrate(index_in=index_in, lam=lam1)
	ReFinAvg = temp['ReFinAvg']
	ReFinErr = temp['ReFinErr']
	ImFinAvg = temp['ImFinAvg']
	ImFinErr = temp['ImFinErr']
	ReDivAvg = temp['ReDivAvg']
	ReDivErr = temp['ReDivErr']
	ImDivAvg = temp['ImDivAvg']
	ImDivErr = temp['ImDivErr']
	time  = temp['TIME']
	fields = ['diagrams', 'squark', 'who', 'ReFinAvg', 'ReFinErr', 'ReDivAvg', 'ReDivErr', 'ImFinAvg', 'ImFinErr', 'ImDivAvg', 'ImDivErr','msq1', 'msq2', 'mA', 'tanb', 'Aq', 'lam', 'eps','startdiv', 'enddiv', 'startfin', 'endfin', 'TIME']
	diagrams = 'D'+str(dia_index[index_in])
	who = 'KIT'

	if richardson == True:
		filename = diagrams+'_'+str(lam1)+'.csv'
	else:
		filename = diagrams+'_'+str(name)+'.csv'
	
	if (squark == 1):
		squarks = 'stop'
	elif (squark == 2):
		squarks = 'sbot'

	rows = [[str(dia_index[index_in]), str(squark), who, str(ReFinAvg), str(ReFinErr), str(ReDivAvg), str(ReDivErr), str(ImFinAvg), str(ImFinErr), str(ImDivAvg), str(ImDivErr), str(msq1_ns), str(msq2_ns), str(mA.real), str(tanb), str(At.real), str(lam1), str(eps), str(startdiv), str(enddiv), str(startfin), str(endfin), str(time)]]
	with open(working+'/'+filename, 'w') as csvfile:
		csvwriter = csv.writer(csvfile)
		csvwriter.writerow(fields)
		csvwriter.writerows(rows)
	csvfile.close

if __name__ == '__main__':
	if True:
		args = parse_arguments()
		function_dir = args.amp_dir
		name = args.input.name
		index_in = args.index_in
		lam = args.lamb
		Richardson = params.richardson
		csvwrite(working=str(working_dir), name=str(name), index_in=index_in, lam1=lam, richardson=Richardson)
	else:
		import hotshot, hotshot.stats
		prof = hotshot.Profile('vegas.prof')
		prof.runcall(main)
		prof.close()
		stats = hotshot.stats.load('vegas.prof')
		stats.strip_dirs()
		stats.sort_status('time', 'calls')
		stats.print_stats(40)
