# GGA_MSSM_SQCD

In order to use richardson flag, user has to set the richardson flag in .json file to true and pick a value for order k of the extrapolation. Then,

python3 main.py.

By default, the runs will submit to each order to cluster using sbatch.

All diagrams are stable below and above quark threshold up to lam=10^-5, where lam is the imaginary shift of quark and squark masses.

# Update: 

Real and QCD corrections part for xs from sushi.

Delete current version of sqcd corrections in xs.

Added summaryfile.csv file for CSQCD factor normalized to LO factor. 

# Install:

For pylooptools, 

python3 setup.py install

In ltools.py, DLL_NAME is replaced with a path to looptools.so. To test,

python3 -m import "from pylooptools.ltools import *"

# External package:

https://github.com/djukanovic/pylooptools

https://github.com/JohannesBuchner/PyMultiNest

https://sushi.hepforge.org/