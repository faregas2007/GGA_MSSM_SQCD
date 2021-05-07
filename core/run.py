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

# get module of function file
function_dir = sys.argv[2]
index_in = int(sys.argv[3])
#dia_index = sys.argv[3]

dia_index = [2, 3, 5, 7, 9, 11, 12, 13, 14]
sys.path.insert(0, function_dir)
fcn = 'ampibp_total_dia' + str(dia_index[index_in]) + '.py'
module = fcn.strip('.py')
fcn_module = import_module(module)
ibpptotal=getattr(fcn_module, 'ibpptotal')
ibpftotal=getattr(fcn_module, 'ibpftotal')

def main():
	np.random.seed((1,2,3))
	fname=['real', 'imag']
	integ = vegas.Integrator(5*[[eps, 1. - eps]], nhcube_batch=2000, sync_ran=True)

	for i in range(2):
		start = timer()
		integ(ibpptotal(i), nitn=iters, neval=startdiv, alpha=alphafix)
		result = integ(ibpptotal(i), nitn=iters, neval=enddiv, alpha=alphafix)
		end = timer()
		print('diagram%s'%dia_index[index_in]+' '+fname[i]+' '+'pole=%s'%result+' '+'time=%s'%(end-start)+'\n')
		print(result.summary())

	for i in range(2):
		start = timer()
		integ(ibpftotal(i), nitn=iters, neval=startfin, alpha=alphafix)
		result = integ(ibpftotal(i), nitn=iters, neval=endfin, alpha=alphafix)
		end = timer()
		print('diagram%s'%dia_index[index_in]+' '+fname[i]+' '+'finite=%s'%result+' '+'time=%s'%(end-start)+'\n')
		print(result.summary())

if __name__ == '__main__':
	if True:
		main()
	else:
		import hotshot, hotshot.stats
		prof = hotshot.Profile('vegas.prof')
		prof.runcall(main)
		prof.close()
		stats = hotshot.stats.load('vegas.prof')
		stats.strip_dirs()
		stats.sort_status('time', 'calls')
		stats.print_stats(40)
