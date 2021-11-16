import sys, os
from pathlib import Path

from init import *
from evalsusy import *
from lhapdfs import *

from timeit import default_timer as timer

start = timer()
#print(inputs().message())

print(einital().get_json())

# mb will be initialize as mbos depending on mbrunloop and mbrunyuk
# mb everywhere in program will be mbos, determine via getmass method.
# mb will always the crucial inputs parameters.
#test_config={'test':1}
print(quarkhiggscoup().get_json())
#print(quarkhiggscoup().update(**test_config))

thetat = 0.004
thetab = 0.001
mb = 4.908

print(squarkhiggscoupCMSSM().get_coup(mb=mb, thetat=thetat, thetab=thetab))

x1 = 0.9976
x2 = 3.3366e-3

print(sushi_pdfs().PDFgg(x1,x2))

import sys
#sys.path.insert(0, '/users/tp/dat/aggrenorm2/src/renschemes/')


from renschemes.squarkmatrix import *
print(squarkmatrix().get_json())
print(squarkmatrix().deltaML())

from renschemes.sbottomrenorm import *
print(f'sbottomrenormHM: {sbottomrenorm().get_json()}')
print(f"sbottomrenormDS: {sbottomrenorm().get_json(obj='sbottomrenormDS')}")
print(f'dmbfin: {sbottomrenorm().dmbfin()}')
print(f"mb_sbottom: {sbottomrenorm().mb_sbottom()}")


from renschemes.mbrenorm import *
print(f'dmbnontanb: {mbrenorm().dmbnontanb()}')
print(f'dmbtanb: {mbrenorm().dmbtanb()}')
print(f'mbMSDRtrans: {mbrenorm().mbMSDRtrans(mu_in=1400.0)}')
print(f'mbMSbarmuX: {mbrenorm().mbMSbarmuX(mu_in=1400.0)}')
print(f'deltab:{mbrenorm().deltab()}')
#print(f'mb_sbottom: {mbrenorm().mb_sbottom()}')

from renschemes.Abrenorm import *
print(f'dAbosfintbHM: {Abrenorm().dAbosfintbHM()}')
print(f'dAbosfinmbHM: {Abrenorm().dAbosfinmbHM()}')

from renschemes.renschemes import *
print(f'renormalize: {renormalize().get_json()}')
print(f'renormalizeSM: {renormalizeSM().get_json()}')
#print(f'renschemes: {renschemes(inputs, einital, renormalizeSM, quarkhiggscoup, squarkhiggscoupCMSSM, renormalize).get_params()}')
print(f'renschemes: {renschemes(inputs, einital, renormalizeSM, quarkhiggscoup, squarkhiggscoupCMSSM, renormalize).to_file()}')

"""
from xs.reals import *
s = 13000**2
tpu = 12
x = 0.1
x1 = 0.11
x2 = 1e-9
pt = 0.0
y = 0.0
s = 151911348.39435318
t = -123602374.64065668
u = -28268973.753696509
"""
"""
s = np.full(288, s)
t = np.full(288, t)
u = np.full(288, u)
print(f'AMPgg: {real_amps().AMPgg_batch(s,t,u)}')
"""
"""
print(f'AMPqq: {real_amps().AMPqq(s, tpu)}')
print(f'AMPqg: {real_amps().AMPqg(s, t, u)}')
print(f'AMPgg: {real_amps().AMPgg(s, t, u)}')
"""


from xs.integrations import *
"""
print(f'AMPqqjet: {qq().AMPqqjet(s=s, t=t, u=u)}')
"""
"""
print(f'qqinteg: {qq().qqinteg()}')
print(f'qginteg: {qg().qginteg()}')
print(f'gginteg: {gg().gginteg()}')
"""

"""
from xs.integrations import *

print(f'qqinteg: {qqinteg_batch().qqinteg()}')
"""


from xs.sigma import *
print(f'sigmasqcd: {sigma_sqcd(qq, qg, gg).get_json()}')

"""
print(f'sigma: {sigma(qq, qg, gg).get_json()}')
"""
"""
from xs.integrals2 import *
m1 = 10.0
m2 = 20.0
m3 = 95.0
mb2 = 20.0
print(f'Integral4: {Integral4(m1,m2,m3,mb2)}')

ma = m1*mb2
mb = m2*mb2
mc = m3*mb2
print(f'Integral4_sushi:{Integral4_sushi(ma,mb,mc,mb2)}')
"""
"""
print(f'qginteg: {qginteg_batch().qginteg()}')
print(f'gginteg:{gginteg_batch().gginteg()}')
"""

"""
end = timer()
"""

"""
from xs.sigma import *
print(f'sigma: {sigma(qqinteg_batch, qginteg_batch, gginteg_batch).get_json()}')
"""

end = timer()
print(f'{end - start}s')


"""
 # non-batch mdoe
        if(t.shape == ()):
            if (t >= 0.0):
                t = -s*10**-16
            if (u >= 0.0):
                u = -s*10**-16
        # batch mode
        else:
            ids = np.where(t>0.0)
            #t[ids] = -s[ids]*10**-16
            t[ids] = 0.0
            ids = np.where(u>0.0)
            #u[ids] = -s[ids]*10**-16
            u[ids] = 0.0
        
        tmp = real_amps().AMPgg_batch(s,t,u)*self.jet(pt,y)
        #print(f'AMPgg_batch before: {real_amps().AMPgg_batch(s,t,u)}')
        if (self._subtr):
            if(t.size == 1):
                if ((-s/t > 10**5) or (-s/u > 10**5)):
                    return 0.0
            else:
                total = np.vstack((tmp, -s/t, -s/u))
                ids = np.where(total[0,::]>0.0)
                ids1 = np.where(total[1,::] > 10**5)
                ids2 = np.where(total[2,::] > 10**5)
                print(f'ids1:{ids1}')
                print(f'ids2:{ids2}')
                if(ids1[0].size >= 1 or ids2[0].size >= 1):
                    ids3 = np.union1d(ids1[0],ids2[0])
                    ids4 = np.setdiff1d(ids[0],ids3[0])
                    print(f'ids3:{ids3}')
                    tmp[ids3] = 0.0
                    
                    s = s[ids4]
                    t = t[ids4]
                    u = u[ids4]
                    tmp[ids4] -= self.subggt(s,t,u)*self.jetsub(x1*x, x2) - self.subggu(s,t,u)*self.jetsub(x1, x2*x)
                else:
                    tmp -= self.subggt(s,t,u)*self.jetsub(x1*x, x2) - self.subggu(s,t,u)*self.jetsub(x1, x2*x)
                    #print(f'tmpzerosize:{tmp[tmp==0.0].size}')
            
        return x*(1-x)*tmp*sushi_pdfs().PDFgg_batch(x1, x2)/(x1*x2)
"""