from __future__ import division
import vegas
import numpy as np
import cmath as cm
import sys
import gvar as gv
import pyslha
from slha_input import *

"""
islha_obj = pyslha.read(sys.argv[1])
lam = islha_obj.blocks["IBPPARA"][1]
mA = islha_obj.blocks["MASS"][1] + complex(0.,0.)
mt = islha_obj.blocks["MASS"][2]*np.sqrt(1. - np.complex(0., lam))
mtmsbar = islha_obj.blocks["MASS"][3] + np.complex(0.,0.)
mG = islha_obj.blocks["MASS"][4] + np.complex(0.,0.)
msq1 = islha_obj.blocks["MASS"][5]*np.sqrt(1. - np.complex(0., lam))
msq2 = islha_obj.blocks["MASS"][6]*np.sqrt(1. - np.complex(0., lam))
mu = islha_obj.blocks["REAL"][1]
At = islha_obj.blocks["REAL"][2] + complex(0.,0.)
beta = islha_obj.blocks["REAL"][3]
t = islha_obj.blocks["REAL"][4]
a = islha_obj.blocks["REAL"][5]
b = islha_obj.blocks["REAL"][6]
"""

class ibpftotal(vegas.BatchIntegrand):
	def __init__(self, reim, lam):
		self.reim = reim
		self.lam = lam
	def __call__(self, x):
		mt = mt_ns*np.sqrt(1. - complex(0., self.lam))
        msq1 = msq1_ns*np.sqrt(1. - complex(0., self.lam))
        msq2 = msq2_ns*np.sqrt(1. - complex(0., self.lam))

		if (self.reim == 0):
			return np.ascontiguousarray(((complex(0.,1.)*mA**2*((((-3*complex(0.,1.))/2)*mt*complex(1.,0.)*x[:,0]**2*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(a**2*mt*(-1 + complex(1., 0.)*x[:,3]) + b**2*mt*(-1 + complex(1., 0.)*x[:,3]) + 2*a*b*mG*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*np.arctan((mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,0] - 2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1] - complex(1.,0.)*x[:,4] + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))/np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq1**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))))/((-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq1**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))) - (((3*complex(0.,1.))/2)*mt*complex(1.,0.)*x[:,0]**2*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(a**2*mt*(-1 + complex(1., 0.)*x[:,3]) + b**2*mt*(-1 + complex(1., 0.)*x[:,3]) + 2*a*b*mG*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*np.arctan((mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4]))))/np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq1**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))))/((-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq1**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))) - (((3*complex(0.,1.))/2)*mt*complex(1.,0.)*x[:,0]**2*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(a**2*mt*(-1 + complex(1., 0.)*x[:,3]) + b**2*mt*(-1 + complex(1., 0.)*x[:,3]) - 2*a*b*mG*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*np.arctan((mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,0] - 2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1] - complex(1.,0.)*x[:,4] + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))/np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq2**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))))/((-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq2**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))) - (((3*complex(0.,1.))/2)*mt*complex(1.,0.)*x[:,0]**2*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(a**2*mt*(-1 + complex(1., 0.)*x[:,3]) + b**2*mt*(-1 + complex(1., 0.)*x[:,3]) - 2*a*b*mG*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*np.arctan((mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4]))))/np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq2**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))))/((-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq2**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2)))))/(mt**2*((-complex(0.,1.))*np.pi + np.log((1 + np.sqrt(1 - (4*mt**2)/mA**2))/(1 - np.sqrt(1 - (4*mt**2)/mA**2))))**2)).real)
		elif (self.reim == 1):
			return np.ascontiguousarray(((complex(0.,1.)*mA**2*((((-3*complex(0.,1.))/2)*mt*complex(1.,0.)*x[:,0]**2*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(a**2*mt*(-1 + complex(1., 0.)*x[:,3]) + b**2*mt*(-1 + complex(1., 0.)*x[:,3]) + 2*a*b*mG*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*np.arctan((mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,0] - 2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1] - complex(1.,0.)*x[:,4] + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))/np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq1**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))))/((-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq1**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))) - (((3*complex(0.,1.))/2)*mt*complex(1.,0.)*x[:,0]**2*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(a**2*mt*(-1 + complex(1., 0.)*x[:,3]) + b**2*mt*(-1 + complex(1., 0.)*x[:,3]) + 2*a*b*mG*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*np.arctan((mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4]))))/np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq1**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))))/((-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq1**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))) - (((3*complex(0.,1.))/2)*mt*complex(1.,0.)*x[:,0]**2*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(a**2*mt*(-1 + complex(1., 0.)*x[:,3]) + b**2*mt*(-1 + complex(1., 0.)*x[:,3]) - 2*a*b*mG*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*np.arctan((mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,0] - 2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1] - complex(1.,0.)*x[:,4] + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))/np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq2**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))))/((-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq2**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))) - (((3*complex(0.,1.))/2)*mt*complex(1.,0.)*x[:,0]**2*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(a**2*mt*(-1 + complex(1., 0.)*x[:,3]) + b**2*mt*(-1 + complex(1., 0.)*x[:,3]) - 2*a*b*mG*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*np.arctan((mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4]))))/np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq2**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2))))/((-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*np.sqrt(mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(4*(1 + (-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*complex(1., 0.)*x[:,3])*(mt**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3]) - msq2**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1])*(-1 + complex(1., 0.)*x[:,3]) + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*complex(1., 0.)*x[:,3]*(mG**2*(-1 + complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]) + mA**2*(-1 + complex(1.,0.)*x[:,0])*(-1 + complex(1., 0.)*x[:,3])*complex(1.,0.)*x[:,4])) - mA**2*complex(1.,0.)*x[:,0]*complex(1.,0.)*x[:,1]*(-1 + complex(1., 0.)*x[:,3])*(-1 + complex(1., 0.)*x[:,3]*(1 + complex(1.,0.)*x[:,4] - complex(1.,0.)*x[:,0]*(1 + complex(1.,0.)*x[:,1]*complex(1.,0.)*x[:,4])))**2)))))/(mt**2*((-complex(0.,1.))*np.pi + np.log((1 + np.sqrt(1 - (4*mt**2)/mA**2))/(1 - np.sqrt(1 - (4*mt**2)/mA**2))))**2)).imag)

class ibpptotal(vegas.BatchIntegrand):
	def __init__(self, reim, lam):
		self.reim = reim
		self.lam = lam
	def __call__(self, x):
		mt = mt_ns*np.sqrt(1. - complex(0., self.lam))
        msq1 = msq1_ns*np.sqrt(1. - complex(0., self.lam))
        msq2 = msq2_ns*np.sqrt(1. - complex(0., self.lam))
		if (self.reim ==  0):
			return np.ascontiguousarray((0.*x[:,0]*complex(1.,0.)+ 0.*x[:,1]*complex(1.,0.)+0.*x[:,2]*complex(1.,0.)+0.*x[:,3]*complex(1.,0.)+0.*x[:,4]*complex(1.,0.)).real)
		elif (self.reim == 1):
			return np.ascontiguousarray((0.*x[:,0]*complex(1.,0.)+ 0.*x[:,1]*complex(1.,0.)+0.*x[:,2]*complex(1.,0.)+0.*x[:,3]*complex(1.,0.)+0.*x[:,4]*complex(1.,0.)).imag)
