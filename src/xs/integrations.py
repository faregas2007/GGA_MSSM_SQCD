# integrations.py
# slow convergence vegas intgrations --> optimzation with c++ code ? 
# avoid using . for reference variables. 

import vegas

from typing import Dict, List

from xs.integrals import *
from xs.reals import *

from renschemes.virtuals import *

#from renschemes.renschemes import virtuals
from include import *
from consts import *
from init import *
from lhapdfs import *

# facade desgin pattern
# using direct called params method instead of using class inheritance

#class integrators(einital):
class integrators(integ, inputs):
    _z = einital().initial()['z']
    _pseudorap = einital().initial()['pseudorap']
    _minrap = einital().initial()['minrap']
    _maxrap = einital().initial()['maxrap']
    _mh2 = einital().initial()['mh2']
    
    def integrate(self, ndim, integrand)->float:
        np.random.seed((1,2,3))
        integ = vegas.Integrator(ndim*[[eps, 1-eps]], sync_ran=True)

        # adapt vegas into integrand
        integ(integrand, nitn=self.iters, neval=self.nstart, alpha=self.alpha_vegas)
        result = integ(integrand, nitn=self.iters, neval=self.nend, alpha=self.alpha_vegas)

        return result, result.sdev, result.chi2


    def trans(self, integrand, xx)->float:
        z = self._z
        mh2 = self._mh2
        pseudorap = self._pseudorap

        tau1 = xx[0]
        tau2 = xx[1]
        tau3 = xx[2]

        w1 = tau1*tau2*(1-z) + z
        w2 = tau1*(1-z) + z

        x1 = w1/w2
        x2 = z/w1
        jdet = (1-z)**2*tau1*x2/w2**2

        # for space space integration
        s = mh2/w2
        t = mh2*(1-1/w2)*tau3
        u = mh2*(1-1/w2)*(1-tau3)

        if(pseudorap):
            # pesudo-rapidity
            y = np.log(x1*t/x2/u)/2.0
        else:
            # rapidity
            y = np.log(x1*(mh2-t)/x2/(mh2-u))/2.0
        pt = np.sqrt(t*u/s)
        #print(f'pseudorap, y: {pseudorap}, {y}')
        trans = integrand(w2,s,t,u,x1,x2, pt,y)*jdet
        return trans

    def transplus(self, integrand1, integrand2, xx)->float:
        z = self._z
        minrap = self._minrap
        maxrap = self._maxrap
        
        tau1 = xx[0]
        tau2 = xx[1]

        jdet1 = 1-z
        jdet2 = maxrap - minrap
        y = minrap + jdet2*tau2
        x1 = np.exp(y)*np.sqrt(z)
        x2 = z/x1
        w2 = tau1*(1-z)+z
        transplus = -integrand2(w2,x1,x2)*z*jdet1*jdet2
        transplus = transplus*2.0 # add subtraction terms for t->0 and u->0
        
        if(x1 <= w2):
            # subtraction term for t->0
            transplus = transplus + integrand1(w2, x1/w2, x2)*z/w2**2*jdet1*jdet2

        if(x1>= z/w2):
            # subtraction term for u->0
            transplus = transplus + integrand1(w2, x2/w2, x1)*z/w2**2*jdet1*jdet2
        
        transplus = transplus*2.0 # include negative y
        return transplus

    def transdelta(self, integrand, xx):
        z = self._z
        minrap = self._minrap
        maxrap = self._maxrap
        """
        maxrap = -np.log(z)/2.0
        minrap = 0.0
        """

        tau = xx[0]
        jdet1 = maxrap - minrap
        y = minrap + jdet1*tau
        x1 = np.exp(y)*np.sqrt(z)
        jdet2 = x1
        x2 = z/x1
        transdelta = x2*jdet1*jdet2*integrand(x1,x2)
        transdelta = transdelta * 2.0 # to include negative y
        return transdelta

    def jetsub(self, x1, x2):
        if(self._pseudorap):
            jetsub = self.jet(0.0,1.0)
        else:
            jetsub = self.jet(0.0,np.log(x1/x2)/2.0)
        
        return jetsub
    
    def jet(self, pt, y):
        jet = 1.0
        if(self.rapcut and (np.abs(y) > self._maxrap or np.abs(y)<self._minrap)):
            jet = 0.0
        if(self.ptcut and (np.abs(y) > self.maxptc or np.abs(pt)<self.minptc)):
            jet = 0.0
        return jet

    
class qq(integrators):

    def AMPqqjet(self, s, t, u):
        return (t**2 + u**2)/(t+u)**2*3.0/2.0 * real_amps().AMPqq(s, t + u)
    
    def realqqint(self, x,s,t,u,x1,x2,pt,y):
        return (1-x)**3*self.AMPqqjet(s,t,u)*self.jet(pt, y)*sushi_pdfs().PDFqq(x1, x2)/(x1*x2)
    
    def realqq(self, x, x1, x2):
        return self.trans(self.realqqint, [x1,x2,x])

    # should try with batch integrands
    # have to deal with C0 function issue that it does not recognize batchs input. 
    def qqreal(self, x):
        return self.realqq(x[0], x[1], x[2])

    def qqinteg(self):
        hard, error, chi2 = self.integrate(integrand=self.qqreal, ndim=3)
        intqq = hard*32/27.0
        error = error*32/27.0
        return intqq, error


class qg(integrators):
    _lfh = einital().initial()['lfh']
    _subtr = einital().initial()['subtr']
    def intsubqg(self, x1, x2):
        return self.transplus(self.subqgint1, self.subqgint2, [x1,x2])

    def subqgint1(self, x, x1, x2):
        # integrated Catani Seymour subtraction term
        return self.qgterms(x)*sushi_pdfs().PDFqg(x1, x2)/(x1*x2)

    def subqgint2(self, x, x1, x2):
        return 0.0

    def qgterms(self, x):
        return ((1 + (1 - x)**2)*(2*np.log(1 - x)-np.log(x) - self._lfh) + x**2)

    # qg main terms
    def realqgint(self, x,s,t,u,x1,x2, pt, y):
        tmpqg = real_amps().AMPqg(s,t,u)*self.jet(pt,y)
        tmpgq = real_amps().AMPqg(s,u,t)*self.jet(pt,y)

        if(self._subtr):
            tmpqg = tmpqg - self.subqg(s,t,u)*self.jetsub(x1*x, x2)
            tmpgq = tmpgq - self.subqg(s,u,t)*self.jetsub(x1, x2*x)
        return x*(1-x)*(tmpqg*sushi_pdfs().PDFqg(x1, x2)/(x1*x2) + tmpgq*sushi_pdfs().PDFqg(x2, x1)/(x2*x1))

    def subqg(self, s, t, u):
        return -(s**2 + (t+u)**2)/(t*(s+t+u))

    def realqg(self, x, x1,x2):
        return self.trans(self.realqgint, [x1, x2, x])

    def qgreal(self, x):
        return self.realqg(x[0], x[1], x[2])

    def qgsub(self, x):
        return self.intsubqg(x[0], x[1])

    def qginteg(self):
        error = [0.0,0.0]
        chi2 = [0.0,0.0]
        hard, error[0], chi2[0] = self.integrate(ndim=3, integrand=self.qgreal)
        if (self._subtr):
            sub, error[1], chi2[1] = self.integrate(ndim=2, integrand=self.qgsub)
        else:
            sub = 0.0
        intqg = (hard+sub)*cf/2.0
        error = np.sqrt(error[0]**2 + error[1]**2)*cf/2.0
        return intqg, error

class gg(integrators):
    _z = einital().initial()['z']
    _lfh = einital().initial()['lfh']
    _lfr = einital().initial()['lfr']
    _subtr = einital().initial()['subtr']
    # delta terms in gg channel
    def deltaggint(self, x1,x2):
        return sushi_pdfs().PDFgg(x1, x2)/(x1*x2)

    def deltagg(self, x1):
        return self.transdelta(self.deltaggint , [x1])

    # Catani-Seymour subtraction (gg channel)
    def intsubgg(self, x1, x2):
        return  self.transplus(self.subggint1, self.subggint2, [x1,x2])

    def subggint1(self, x,x1,x2):
        return (self.ggterms(x) + self.dterms(x))/2.0 * sushi_pdfs().PDFgg(x1, x2)/(x1*x2)

    def subggint2(self, x,x1,x2):
        return  (self.dterms(x)/2.0) * sushi_pdfs().PDFgg(x1, x2)/(x1*x2)

    def ggterms(self, x):
        return (4*np.log(1-x) - 2*self._lfh)*(-2*x + x**2 - x**3) -2.0*np.log(x)*(1 - x + x**2)**2/(1.0 - x)

    def dterms(self, x):
        dd0 = 1/(1.0 -x)
        dd1 = np.log(1.0 -x)*dd0
        return -2*self._lfh * dd0 + 4*dd1
    
    def dtermsz(self):
        z = self._z
        dd0 = np.log(1-z)/(1-z)
        dd1 = (np.log(1-z)/2.0)*dd0
        return -2*self._lfh * dd0  + 4 *dd1

    def subggt(self, s,t,u):
        return  (s**2 + s*(t+u) + (t+u)**2)**2/(s*t*(t+ u)*(s+t+u))

    def subggu(self, s,t,u):
        return (s**2 + s*(t+ u) + (t+u)**2)**2/(s*u*(t+u)*(s+t+u))

    # real correction gg channel (finite terms)
    def realggint(self, x,s,t,u,x1,x2,pt,y):
        if (t >= 0.0):
            t = -s*10**-6
        if (u >= 0.0):
            u = -s*10**-16

        tmp = real_amps().AMPgg(s,t,u)*self.jet(pt,y)
        if (self._subtr):
            if ((-s/t > 10**5) or (-s/u > 10**5)):
                return 0.0
            tmp = tmp - self.subggt(s,t,u)*self.jetsub(x1*x, x2) - self.subggu(s,t,u)*self.jetsub(x1, x2*x)

        return x*(1-x)*tmp*sushi_pdfs().PDFgg(x1, x2)/(x1*x2)

    def realgg(self, x, x1, x2):
        return self.trans(self.realggint, [x1, x2, x])

    def ggreal(self, x):
        return self.realgg(x[0], x[1], x[2])

    def ggsub(self, x):
        return self.intsubgg(x[0], x[1])

    def ggdelta(self, x):
        return self.deltagg(x[0])

    # for leading order cross-section
    def deltas(self):
        delta, error, chi2 = self.integrate(ndim=1, integrand=self.ggdelta)
        return delta, error, chi2

    def gginteg(self):
        lfr = self._lfr
        error = [0.0,0.0,0.0]
        chi2 = [0.0,0.0,0.0]
        #virtual, virtual_susy, errorv_susy = c_virt()
        virts = virtual().c_virt()

        delta, error[0], chi2[0] = self.integrate(ndim=1, integrand = self.ggdelta)
        sub, error[2], chi2[2] = self.integrate(ndim=2, integrand = self.ggsub)
        hard, error[1], chi2[1] = self.integrate(ndim=3, integrand= self.ggreal)
        errord = error[0]

        reell = ca*(hard+sub)
        errorr = ca*np.sqrt(error[1]**2 + error[2]**2)
        # subtr is defined in sushicore. Its default value is true. 
        if (self._subtr):
            #virt = delta*(virtual+ virtual_susy +  np.pi**2-2*betanull*lfr + ca*dtermsz())
            virt = delta*(virts + Pi*Pi - 2*betanull*lfr + ca*self.dtermsz())
            # relative error
            errorv = error[0]*(virts + Pi*Pi - 2*betanull*lfr + ca*self.dtermsz())
            #errorv_v = errorv +  errorv_susy*delta
        else:
            virt = 0.0
            errorv = 0.0

        intgg = reell + virt
        return intgg, delta, errorr, errorv, errord


