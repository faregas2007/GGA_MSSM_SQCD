from __future__ import absolute_import, unicode_literals, print_function
import sys
import pycuba

import ctypes
from ctypes import *

import numpy as np
from mpmath import *
from timeit import default_timer as timer

def polylog2(x):
	return polylog(2, x)

def polylog3(x):
	return polylog(3, x)


def mul(a, b):
	a = c_double(a)
	if isinstance(b, c_double):
		return c_double(a.value*b.value)
	else:
		return c_double(b*a.value)

def sum(a, b):
	a = c_double(a)
	if isinstance(b, c_double):
		return c_double(a.value + b.value)
	else:
		return c_double(b + a.value)

def transform(a, b, x):
	delta = b - a
	return a + delta*x

def jaco(a, b, d):
	return (b -a)**d


def print_header(name):
	return print('---------------------------%s---------------------------'% name)

def vegas_Integrator(module, ndim, userdata, epsrel, epsabs, verbose, ncomp, seed, mineval, maxeval, nstart,
			nincrease, nbatch, gridno, statefile, nvec):
	from os import environ as env
	if 'CUBAVERBOSE' in env:
		verbose = int(env['CUBAVERBOSE'])

	print_header('Vegas')
	start = timer()
	results = pycuba.Vegas(module, ndim, userdata, epsrel, epsabs, verbose, ncomp, seed,
			mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile, nvec)
	end = timer()
	return {'result':results['results'][0][u'integral'], 'error':results['results'][0]['error'], 'pvalue':results['results'][0][u'prob'], 'time':end-start}

def cuhre_Integrator(module, ndim, key, mineval, maxeval, ncomp, userdata, seed, epsrel, epsabs, verbose,
			statefile, nvec):
	from os import environ as env
	if 'CUBAVERBOSE' in env:
		verbose = int(env['CUBAVERBOSE'])
	print_header('Cuhre')
	start = timer()
	results = pycuba.Cuhre(module, ndim, key, mineval, maxeval, ncomp, userdata, seed,
			epsrel, epsabs, verbose, statefile, nvec)
	end = timer()
	return {'result':results['results'][0][u'integral'], 'error':results['results'][0]['error'], 'pvalue':results['results'][0][u'prob'], 'time':end-start}

def suave_Integrator(module, ndim, nnew, nmin, flatness, userdata, epsrel, epsabs, verbose, ncomp, seed, mineval,
			maxeval, statefile, nvec):
	from os import environ as env
	if 'CUBAVERBOSE' in env:
		verbose = int(env['CUBAVERBOSE'])

	print_header('Suave')
	start = timer()
	results = pycuba.Suave(module, ndim, nnew, nmin, flatness, userdata,
			epsrel, epsabs, verbose, ncomp, seed, mineval, maxeval, statefile, nvec)
	end = timer()
	return {'result':results['results'][0][u'integral'], 'error':results['results'][0]['error'], 'pvalue':results['results'][0][u'prob'], 'time':end-start}

def divonne_Integrator(module, ndim, key1, key2, key3, maxpass, border, maxchisq, mindeviation,
			mineval, maxeval, ncomp, ldxgiven, xgiven, nextra, peakfinder, userdata, seed, epsrel,
			epsabs, verbose, statefile, nvec):
	from os import environ as env
	if 'CUBAVERBOSE' in env:
		verbose = int(env['CUBAVERBOSE'])

	print_header('Divonne')
	start = timer()
	results = pycuba.Divonne(module, ndim, key1, key2, key3, maxpass, border,
		maxchisq, mindeviation, mineval, maxeval, ncomp, ldxgiven, xgiven, nextra, peakfinder,
		userdata, seed, epsrel, epsabs, verbose, statefile, nvec)
	end = timer()
	return {'result':results['results'][0][u'integral'], 'error':results['results'][0]['error'], 'pvalue':results['results'][0][u'prob'], 'time':end-start}

def Integrator(name, module, ndim, key, mineval, maxeval, verbose):

	#common flags
	NULL = ctypes.POINTER(c_int)()
	NCOMP = 1
	NVEC = 1
	STATEFILE = NULL
	USERDATA = NULL
	SEED = None
	EPSREL = 0.001
	EPSABS = 1e-12

	# vegas
	NCOMP = 1
	SEED = None
	NSTART = 1000
	NINCREASE = 500
	NBATCH = 1000
	GRIDNO = 0

	#Cuhre
	KEY=0
	USERDATA = NULL


	#SUAVE
	NNEW = 1000
	NMIN = 2
	FLATNESS = 50.0

	#Divonne
	KEY1 = 47
	KEY2 = 2
	KEY3 = 1
	MAXPASS = 5
	BORDER = 0.
	MAXCHISQ = 10.
	MINDEVIATION = 0.25
	XGIVEN = None
	LDXGIVEN = None
	NEXTRA = 0
	PEAKFINDER = None


	if name == 'Vegas':
		return vegas_Integrator(module, ndim, USERDATA, EPSREL, EPSABS, verbose, NCOMP, SEED,
				mineval, maxeval, NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE, NVEC)
	elif name == 'Cuhre':
		return cuhre_Integrator(module, ndim, key, mineval, maxeval, NCOMP, USERDATA, SEED,
				EPSREL, EPSABS, verbose, STATEFILE, NVEC)
	elif name == 'Suave':
		return suave_Integrator(module, ndim, NNEW, NMIN, FLATNESS, USERDATA, EPSREL, EPSABS, verbose, NCOMP,
				SEED, mineval, maxeval, STATEFILE, NVEC)
	elif name == 'Divonne':
		return divonne_Integrator(module, ndim, KEY1, KEY2, KEY3, MAXPASS, BORDER,
					MAXCHISQ, MINDEVIATION, mineval, maxeval, NCOMP, LDXGIVEN, XGIVEN, NEXTRA,
					PEAKFINDER, USERDATA, SEED, EPSREL, EPSABS, verbose, STATEFILE, NVEC)
