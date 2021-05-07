import numpy as np
from mpmath import *
mp.dps = 12
mp.pretty = True

def polylog2(x):
	polylog_array = np.frompyfunc(polylog, 2, 1)
	return polylog_array(2, x).astype('cdouble')

def polylog3(x):
	polylog_array = np.frompyfunc(polylog, 2, 1)
	return poylog_array(3, x).astype('cdouble')


