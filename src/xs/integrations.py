# integrations.py
# slow convergence vegas intgrations --> optimzation with c++ code ??? 
# avoid using . for referenced variables. 

import vegas
import gvar as gv

from typing import Dict, List

from xs.integrals import *
from xs.reals import *

from renschemes.virtuals import *

#from renschemes.renschemes import virtuals
from interface import *
from consts import *
from init import *
from lhapdfs import *

# facade desgin pattern

#class integrators(einital):
class integrators(integ, inputs):
    _z = einital().initial()['z']
    _pseudorap = einital().initial()['pseudorap']
    _minrap = einital().initial()['minrap']
    _maxrap = einital().initial()['maxrap']
    _mh2 = einital().initial()['mh2']
    _rapcut = einital().initial()['rapcut']
    _ptcut = einital().initial()['ptcut']
    
    def integrate(self, ndim, integrand):
        np.random.seed((1,2,3))
        integ = vegas.Integrator(ndim*[[eps, 1-eps]], sync_ran=True)

        # adapt vegas into integrand
        integ(integrand, nitn=self.iters, neval=self.nstart, alpha=self.alpha_vegas)
        result = integ(integrand, nitn=self.iters, neval=self.nend, alpha=self.alpha_vegas)

        return result.mean, result.sdev, result.chi2


    def trans(self, integrand, xx):
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

        trans = integrand(w2,s,t,u,x1,x2, pt,y)*jdet
        return trans

    def transplus(self, integrand1, integrand2, xx):
        z = self._z
        minrap = self._minrap
        maxrap = self._maxrap
        
        tau1 = np.array(xx[0])
        tau2 = np.array(xx[1])

        jdet1 = 1-z
        jdet2 = maxrap - minrap
        y = minrap + jdet2*tau2
        x1 = np.exp(y)*np.sqrt(z)
        x2 = z/x1
        w2 = tau1*(1-z)+z
        
        # batch mode modiying-> fast:
        if(tau1.shape == ()):
            transplus = -integrand2(w2, x1, x2)*z*jdet1*jdet2
            transplus = transplus*2.0 # add subraction terms for t->0 and u->0
            if(x1 <= w2):
                # subtraction term for t->0
                transplus += integrand1(w2, x1/w2, x2)*z/w2**2*jdet1*jdet2

            if(x1>= z/w2):
                # subtraction term for u->0
                transplus += integrand1(w2, x2/w2, x1)*z/w2**2*jdet1*jdet2
        else:
            # batch mode modifying
            transplus = np.zeros((tau1.size), )
            transplus = np.full(transplus.size ,-integrand2(w2, x1, x2)*z*jdet1*jdet2)
            transplus = transplus*2.0
            
            # subtraction term for t->0
            ids1 = np.where(x1-w2<=0.0)
            ids2 = np.where(x1-z/w2>=0.0)
            if(ids1[0].size>0): 
                X1 = x1[ids1]
                X2 = x2[ids1]
                W2 = w2[ids1] 
                transplus[ids1] += integrand1(W2, X1/W2, X2)*z/W2**2*jdet1*jdet2
            if(ids2[0].size > 0):
                W2 = w2[ids2]
                X1 = x1[ids2]
                X2 = x2[ids2]
                transplus[ids2] +=  integrand1(W2, X2/W2, X1)*z/W2**2*jdet1*jdet2

        transplus = transplus*2.0 # include negative y
        return transplus

    def transdelta(self, integrand, xx):
        z = self._z
        minrap = self._minrap
        maxrap = self._maxrap

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
        if(self._rapcut and (np.abs(y) > self._maxrap or np.abs(y)<self._minrap)):
            jet = 0.0
        if(self._ptcut and (np.abs(y) > self.maxptc or np.abs(pt)<self.minptc)):
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


class qqbatch(integrators, vegas.BatchIntegrand):
    def __init__(self, dim):
        self.dim = dim
        
    def AMPqqjet(self, s, t, u):
        return (t**2 + u**2)/(t+u)**2*3.0/2.0 * real_amps().AMPqq_batch(s, t + u)

    def realqqint(self, x,s,t,u,x1,x2,pt,y):
        return (1-x)**3*self.AMPqqjet(s,t,u)*self.jet(pt, y)*sushi_pdfs().PDFqq_batch(x1, x2)/(x1*x2)
    
    def realqq(self, x, x1, x2):
        return self.trans(self.realqqint, [x1,x2,x])
    
    def __call__(self, x):
        return self.realqq(x[:,0], x[:,1], x[:,2])

class qqinteg_batch(integrators):
    
    def qqinteg(self):
        qqreal = qqbatch(dim=3)
        hard, error, chi2 = self.integrate(integrand=qqreal, ndim=qqreal.dim)
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

class qgbatch(integrators, vegas.BatchIntegrand):
    _lfh = einital().initial()['lfh']
    _subtr = einital().initial()['subtr']
    
    def __init__(self, dim):
        self.dim = dim

    def intsubqg(self, x1, x2):
        return self.transplus(self.subqgint1, self.subqgint2, [x1,x2])

    def subqgint1(self, x, x1, x2):
        # integrated Catani Seymour subtraction term
        return self.qgterms(x)*sushi_pdfs().PDFqg_batch(x1, x2)/(x1*x2)

    def subqgint2(self, x, x1, x2):
        return 0.0

    def qgterms(self, x):
        return ((1 + (1 - x)**2)*(2*np.log(1 - x)-np.log(x) - self._lfh) + x**2)

    # qg main terms
    def realqgint(self, x,s,t,u,x1,x2, pt, y):
        tmpqg = real_amps().AMPqg_batch(s,t,u)*self.jet(pt,y)
        tmpgq = real_amps().AMPqg_batch(s,u,t)*self.jet(pt,y)

        if(self._subtr):
            tmpqg = tmpqg - self.subqg(s,t,u)*self.jetsub(x1*x, x2)
            tmpgq = tmpgq - self.subqg(s,u,t)*self.jetsub(x1, x2*x)
        return x*(1-x)*(tmpqg*sushi_pdfs().PDFqg_batch(x1, x2)/(x1*x2) + tmpgq*sushi_pdfs().PDFqg_batch(x2, x1)/(x2*x1))

    def subqg(self, s, t, u):
        return -(s**2 + (t+u)**2)/(t*(s+t+u))

    def realqg(self, x, x1,x2):
        return self.trans(self.realqgint, [x1, x2, x])

    def __call__(self, x):
        if(self.dim == 2):
            return self.intsubqg(x[:,0], x[:,1])
        elif(self.dim == 3):
            return self.realqg(x[:,0], x[:,1], x[:,2])
        

class qginteg_batch(integrators):
    _subtr = einital().initial()['subtr']
    def qginteg(self):
        qgreal = qgbatch(dim=3)
        qgsub =  qgbatch(dim=2)
        error = [0.0,0.0]
        chi2 = [0.0,0.0]
        hard, error[0], chi2[0] = self.integrate(integrand=qgreal, ndim=qgreal.dim)
        if (self._subtr):
            sub, error[1], chi2[1] = self.integrate(integrand=qgsub, ndim=qgsub.dim)
        else:
            sub = 0.0
        print(f'hard, sub: {hard}, {sub}')
        intqg = (hard+sub)*cf/2.0
        error = np.sqrt(error[0]**2 + error[1]**2)*cf/2.0
        return intqg, error

class gg(integrators):
    _z = einital().initial()['z']
    _lfh = einital().initial()['lfh']
    _lfr = einital().initial()['lfr']
    _subtr = einital().initial()['subtr']
    _alphasggh = einital().initial()['alphasggh']
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
            t = -s*10**-16
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
        virts = virtual().c_virts()
        #virts_sqcd, e_virts_sqcd = vritual().c_vrit_sqcd()

        delta, error[0], chi2[0] = self.integrate(ndim=1, integrand = self.ggdelta)
        sub, error[2], chi2[2] = self.integrate(ndim=2, integrand = self.ggsub)
        hard, error[1], chi2[1] = self.integrate(ndim=3, integrand= self.ggreal)
        errord = error[0]
        reell = ca*(hard+sub)
        errorr = ca*np.sqrt(error[1]**2 + error[2]**2)
        if (self._subtr):
            #virt = delta*(virtual+ virtual_susy +  np.pi**2-2*betanull*lfr + ca*dtermsz())
            virt = delta*(virts + Pi*Pi - 2*betanull*lfr + ca*self.dtermsz())
            # relative error
            errorv = error[0]*(virts + Pi*Pi - 2*betanull*lfr + ca*self.dtermsz())
            #errorv_v = errorv +  errorv_susy*delta
        else:
            virt = 0.0
            errorv = 0.0
        logger.info(f"delta, sub, hard, virts, dtermsz: {delta, sub, hard, virts, self.dtermsz()}")
        intgg = reell + virt
        return intgg, delta, errorr, errorv, errord

    # A single value + array values. 
    def gginteg_sqcd(self):
        lfr = self._lfr
        alphasggh = self._alphasggh

        error = np.array([0.0,0.0,0.0])
        chi2 = np.array([0.0,0.0,0.0])
        virts_qcd = virtual().c_virts()
        virts_sqcd, e_virts_sqcd = virtual().c_virts_sqcd(alphas=alphasggh)
        
        delta, error[0], chi2[0] = self.integrate(ndim=1, integrand=self.ggdelta)
        sub, error[2], chi2[2] = self.integrate(ndim=2, integrand=self.ggsub)
        hard, error[1], chi2[1] = self.integrate(ndim=3, integrand=self.ggreal)

        errord = error[0]
        reell = ca*(hard + sub)
        errorr = ca*np.sqrt(error[1]*error[1] + error[2]*error[2])

        if(self._subtr):
            print(f'virts_qcd, e_vrtis_sqcd: {virts_qcd, e_virts_sqcd}')
            virt = delta*(virts_qcd + virts_sqcd + Pi*Pi - 2*betanull*lfr + ca*self.dtermsz())
            print(f'virts_qcd + virts_sqcd, error[0]: {virts_qcd + virts_sqcd + Pi*Pi, error[0]}')
            errorv = error[0]*(virts_qcd + virts_sqcd + Pi*Pi - 2*betanull*lfr + ca*self.dtermsz())
            errorv = errorv +  e_virts_sqcd*delta
        else:
            vrit = 0.0
            errorv = 0.0
        
        logger.info(f'delta, sub, hard, virts_qcd, virts_sqcd, e_virts_sqcd: {delta, sub, hard, virts_qcd, virts_sqcd, e_virts_sqcd}')
        intgg = reell + virt
        return intgg, delta, errorr, errorv, errord

class ggbatch(integrators, vegas.BatchIntegrand):
    _z = einital().initial()['z']
    _lfh = einital().initial()['lfh']
    _lfr = einital().initial()['lfr']
    _subtr = einital().initial()['subtr']
    def __init__(self, dim):
        self.dim = dim

    # delta terms in gg channel
    def deltaggint(self, x1,x2):
        return sushi_pdfs().PDFgg_batch(x1, x2)/(x1*x2)

    def deltagg(self, x1):
        return self.transdelta(self.deltaggint , [x1])

    # Catani-Seymour subtraction (gg channel)
    def intsubgg(self, x1, x2):
        return  self.transplus(self.subggint1, self.subggint2, [x1,x2])

    def subggint1(self, x,x1,x2):
        return (self.ggterms(x) + self.dterms(x))/2.0 * sushi_pdfs().PDFgg_batch(x1, x2)/(x1*x2)

    def subggint2(self, x,x1,x2):
        return  (self.dterms(x)/2.0) * sushi_pdfs().PDFgg_batch(x1, x2)/(x1*x2)

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

    # real correction gg channel (finite terms) --> batch modifiying 
    # potential bugs --> got a differs result with non-batch version.
    def realggint(self, x,s,t,u,x1,x2,pt,y):
        ids = np.where(t>=0.0)
        t[ids] = -s[ids] * 1e-16
        ids = np.where(u>=0.0) 
        u[ids] = -s[ids] * 1e-16
        
        tmp = real_amps().AMPgg_batch(s,t,u)*self.jet(pt, y)
        # subtr always true --> move tmp inside if loop is okay.
        # problem when subtr --> false.
        if(self._subtr):
            tmp = np.zeros(tmp.size, )
            ids = np.where(tmp>=0.0)
            ids1 = np.where((-s/t >= 1e5) | (-s/u >= 1e5))
            if(ids1[0].size >= 1):
                tmp[ids1] = 0.0
                ids = np.setdiff1d(ids,ids1)

                s = s[ids] 
                t = t[ids]
                u = u[ids]

                tmp[ids] -= self.subggt(s,t,u)*self.jetsub(x1*x, x2) + self.subggu(s,t,u)*self.jetsub(x1, x2*x)
            else:
                tmp -= self.subggt(s,t,u)*self.jetsub(x1*x, x2) + self.subggu(s,t,u)*self.jetsub(x1, x2*x)
        
        return x*(1-x)*tmp*sushi_pdfs().PDFgg_batch(x1, x2)/(x1*x2)

    def realgg(self, x, x1, x2):
        return self.trans(self.realggint, [x1, x2, x])

    def __call__(self, x):

        if(self.dim == 1):
            return self.deltagg(x[:,0])
        elif(self.dim == 2 ):
            return self.intsubgg(x[:,0], x[:,1])
        elif(self.dim == 3 ):
            return self.realgg(x[:, 0], x[:, 1], x[:,2])



class gginteg_batch(integrators):
    _z = einital().initial()['z']
    _lfh = einital().initial()['lfh']
    _lfr = einital().initial()['lfr']
    _subtr = einital().initial()['subtr']
    def gginteg(self):
        lfr = self._lfr
    
        error = [0.0,0.0,0.0]
        chi2 = [0.0,0.0,0.0]
        #virtual, virtual_susy, errorv_susy = c_virt()
        virts = virtual().c_virts()
        #virts_sqcd, e_virts_sqcd = virtual().c_virt_sqcd()

        ggreal = ggbatch(dim=3)
        ggsub =  ggbatch(dim=2)
        ggdelta = ggbatch(dim=1)

        delta, error[0], chi2[0] = self.integrate(integrand = ggdelta, ndim=ggdelta.dim)
        sub, error[2], chi2[2] = self.integrate(integrand = ggsub, ndim=ggsub.dim)
        hard, error[1], chi2[1] = self.integrate(integrand= ggreal, ndim=ggreal.dim)
        errord = error[0]
        print(f"delta, sub, hard, virts:{delta, sub, hard, virts}")
        reell = ca*(hard+sub)
        errorr = ca*np.sqrt(error[1]**2 + error[2]**2)

        # subtr is defined in sushicore. Its default value is true. 
        if (self._subtr):
            dtermsz = ggreal.dtermsz()
            #virt = delta*(virtual+ virtual_susy +  np.pi**2-2*betanull*lfr + ca*dtermsz())
            virt = delta*(virts + Pi*Pi - 2*betanull*lfr + ca*dtermsz)
            # relative error
            errorv = error[0]*(virts + Pi*Pi - 2*betanull*lfr + ca*dtermsz)
            #errorv_v = errorv +  errorv_susy*delta
        else:
            virt = 0.0
            errorv = 0.0

        intgg = reell + virt
        return intgg, delta, errorr, errorv, errord

