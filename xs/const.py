import numpy as np
global size, n_size, nf, gev2picobarns, CA, CF, GF, betanull, sigmanull, muf, mur, mZ

mZ = 9.11876000e+01

size = 1000.0
nsize = 1000
nf = 5
gev2picobarns = 0.389379660*10**9
CA = 3.0
CF = 4.0/3.0
GF = 1.16638*10**-5
betanull = (11/12.0)*CA - nf/6.0
sigmanullh = GF/(288.0*np.sqrt(2.0)*np.pi)
sigmanullH = GF/(288.0*np.sqrt(2.0)*np.pi)
sigmanullA = GF/(128.0*np.sqrt(2.0)*np.pi)
