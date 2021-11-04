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

from xs.reals import *
s = 13000**2
tpu = 12
t = 1.0
u = 10.0
"""
print(f'AMPqq: {real_amps().AMPqq(s, tpu)}')
print(f'AMPqg: {real_amps().AMPqg(s, t, u)}')
print(f'AMPgg: {real_amps().AMPgg(s, t, u)}')

from xs.integrations import *
print(f'AMPqqjet: {qq().AMPqqjet(s=s, t=t, u=u)}')
print(f'qqinteg: {qq().qqinteg()}')

print(f'qginteg: {qg().qginteg()}')
print(f'gginteg: {gg().gginteg()}')
"""
from xs.sigma import *
print(f'sigma: {sigma(qq, qg, gg).get_json()}')
end = timer()
print(f'{end - start}s')