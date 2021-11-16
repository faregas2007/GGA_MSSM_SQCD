import numpy as np
from pylooptools.ltools import *
#from utils import *
#from const import *
#from slha_input import *
import warnings
from mpmath import *


warnings.filterwarnings('ignore', category=RuntimeWarning)
mp.pretty = True
mp.dps = 12

zeta2 = 1.644934066848226436472415166650
zeta3= 1.20205690315959428539970
zeta4= 1.082323233711138191516003696540
zeta5= 1.03692775514336992633140
zeta6= 1.017343061984449139714517929790
# with QCD loop --> deal with massless internal less from gluon ?

def Integral2_lt(m1, mb2):
    B0(m1,mb2, mb2)
    return B0(m1, mb2, mb2)
    
def Integral3_lt(m1,m2,mb2):
    C0(m1,0.0, m2, mb2, mb2, mb2)
    return C0(m1, 0.0, m2, mb2, mb2, mb2)

def Integral4_lt(m1,m2,m3,mb2):
    #D0(0.0,0.0,0.0, m1+m2+m3, m1,m3 ,mb2,mb2,mb2,mb2)
    return D0(0.0,0.0,0.0, m1+m2+m3, m1,m3 ,mb2,mb2,mb2,mb2)


# sushi used these integrals, without QCD loop and mq2 != 0
def Integral2(ma, mq2):
    m1 = ma/mq2
    if(m1 == 0.0):
        re = -2.0
        im = 0.0
    elif(m1 < 0.0):
        re = -2.0*np.sqrt(1.0 - 4.0/m1)*np.log((np.sqrt(-m1) + np.sqrt(-m1 + 4.0))/2.0)
        im = 0.0
    elif(m1 <= 4.0):
        re = -2.0*np.sqrt(4.0/m1 - 1.0)*np.arcsin(np.sqrt(m1)/2.0)
        im = 0.0
    elif(m1 > 4.0):
        re = -2.0*np.sqrt(1-4.0/m1)*(np.log((np.sqrt(m1) + np.sqrt(m1-4.0))/(2.0)))
        im = 2.0*np.sqrt(1.0 - 4.0/m1)*np.pi/2.0
    else:
        print('Sushi error in Integral2')
        
    Integral2 = complex(re, im)
    #print('Integral2: {Integral2}')
    return Integral2

def Integral3(ma, mb, mq2):
    m1 = ma/mq2
    m2 = mb/mq2

    if (m1 == 0.0):
        re = 0.0
        im = 0.0
    elif (m1 < 0.0):
        re = (np.log(2.0/(2.0-m1 + np.sqrt(m1*(-4+m1))))**2)/4.0
        im = 0.0
    elif (m1 <= 4.0):
        re = - np.arcsin(np.sqrt(m1/4.0))**2
        im = 0.0
    elif (m1 > 4.0):
        re = (np.log((-2.0 + m1+np.sqrt(m1*(-4.0+m1)))/(2.0))**2 - np.pi**2)/4.0
        im = -np.pi*np.log((-2.0 + m1 + np.sqrt(m1*(-4.0 + m1)))/2.0)/(m1-m2)
    
    if (m2 == 0.0):
        re = 0.0
        im = 0.0
    elif (m2 < 0.0):
        re = re - (np.log(2.0/(2.0 - m2 + np.sqrt(m2*(-4+m2))))**2)/4.0
    elif(m2 <= 4.0):
        re = re + np.arcsin(np.sqrt(m2/4.0))**2
    elif(m2 > 4.0):
        re= re - (np.log((-2.0 + m2 + np.sqrt(m2*(-4.0+m2)))/2.0)**2 - np.pi**2)/4.0
        im = im + np.pi*np.log((-2.0 + m2 + np.sqrt(m2*(-4.0+m2)))/(2.0))/(m1-m2)
    else:
        print("Sushi error in Integral3")
    Integral3 = complex(re*2.0/(m1-m2)/mq2, im/mq2)
    #print(f'Integral3, ma, mb, mq2: {Integral3}, {ma}, {mb}, {mq2}')
    return Integral3


def Integral4_sushi(ma, mb, mc, mq2):
    m1 = ma/mq2
    m2 = mb/mq2
    m3 = mc/mq2
    pn = np.array([m1, m1+m2+m3, m3])
    p = [0.0,0.0,0.0]
    tmp = 1 + m1*m3/m2/2.0
    if (tmp**2 == 1.0):
        Integral4 = 1/6.0/mq2**2
        return Integral4
    
    x = [0.0,0.0]
    x[1] = -tmp - np.sqrt(tmp**2-1)
    x[0] = 1/x[1]

    Integral4 = complex(2.0*(2.0*polylog(2, 1+x[1])) + np.log(-x[1])**2/2.0, 0.0)
    for i in range(1,3):
        wrz = np.sqrt(complex((1-pn[i]/2.0)**2-1, 0.0))
        if (pn[i] == 0.0):
            p[i] = 1.0
        elif(pn[i] > 0.0):
            p[i] = 1/(1-pn[i]/2.0 - wrz)
        else:
            p[i] = 1 - pn[i]/2.0 + wrz
        
        for k in range(1,2):
            if(pn[i] < 0.0):
                arg = np.real(p[i]*x[k])
                Integral4 = Integral4 + (-1)**(i+k)*(-2*polylog(2, 1+1/arg) 
                        - np.log(-arg)**2/2)
            elif(pn[i] >= 4.0):
                arg = np.real(p[i]*x[k])
                Integral4 = Integral4 + (-1)**(i+k)*(polylog(2, -1/arg) - zeta2 
                        - np.log(1+1/arg)*np.log(arg) 
                        -(np.log(arg) - complex(0.0, np.pi))**2/4.0
                        -complex(0.0, np.pi)*np.log(1+arg))*2.0
            else:
                Integral4 = Integral4 + (-1)**(i+k)*(polylog(2,1+p[i]*x[k]) + polylog(2,1+x[k]/p[i]))

    Integral4 = Integral4/np.sqrt(tmp**2-1)/m2/mq2**2/2.0
    #print(f'Integral4, ma, mb, mc, mq2 : {Integral4}, {ma}, {mb}, {mc}, {mq2}')
    return Integral4

"""
"""
m1 = 40.0
m2 = 60.0
m3 = 95.0
mb2 = 120.0
print(f'Integral2_lt: {Integral2_lt(m1,mb2)}')
print(f'Integral3_lt: {Integral3_lt(m1,m2,mb2)}')
print(f'Integral4_lt: {Integral4_lt(m1,m2,m3,mb2)}')
"""

"""
m1 = 40.0
m2 = 60.0
m3 = 95.0
mb2 = 120.0
ma = m1*mb2
mb = m2*mb2
mc = m3*mb2
print(f'Integral2_sushi: {Integral2(m1, mb2)}')
print(f'Integral3_sushi: {Integral3(m1, m2, mb2)}')
print(f'Integral4_sushi: {Integral4_sushi(m1,m2,m3,mb2)}')


"""
def cli2(x):
    cli2 = 0.0
    zero = 10**-40
    xr = np.real(x)
    xi = np.imag(x)
    r2 = xr*xr+xi*xi
    if (r2 <= zero):
        cli2 = x
    
    rr = xr/r2
    if (r2 == 1.0 and xi == 0.0):
        if (xr == 1.0):
            cli2 = np.complex(zeta2, 0.0)
        else:
            cli2 = -np.complex(zeta2/2, 0.0)
    
    elif(r2 > 1.0 and rr > 0.5):
        y = (x - 1.0)/x
        cli2 = sushi_cli2(y) + zeta2 - np.log(x)*cdlog(1.0 - x) + 0.5*np.log(x)**2
    elif(r2 > 1.0 and rr < 0.5):
        y = 1.0/x
        cli2 = -sushi_cli2(y) - zeta2 - 0.5*np.log(-x)**2
    elif(r2 < 1.0 and xr > 0.5):
        y = 1.0-x
        cli2 = -sushi_cli2(y) + zeta2 - np.log(x)*np.log(1.0-x)**2
    elif(r2 < 1.0 and xr < 0.5):
        y = x
        cli2 = sushi_cli2(y)
    return cli2

def sushi_cli2(x):
    # talyor expansion for complex dilogarithm (spence function)
    n=nber-1
    z = -np.log(1.0-x)
    sushi_cli2 = B2(nber)
    for i in range(n,1,-1):
        sushi_cli2 = z*sushi_cli2+B2[i]

    sushi_cli2 = z**2*2*sushi_cli2+z
    return sushi_cli2
"""

# real dilogarithm (from HqT)
"""
def dli2(x):
    ZERO = set_precision(0.0)
    ONE = set_precision(1.0)
    HALF = set_precision(0.5)
    MALF = set_precision(-0.5)
    MONE = set_precision(-1.0)
    MTWO = set_precision(-2.0)

    PI6 = 1.6449340668482260
    PI3 = 3.2898681336964530

    c = [0.42996693560813700, 0.40975987533077110,-0.01858843665014600, 0.00145751084062270,
        -0.00014304184442340, 0.00001588415541880,-0.00000190784959390, 0.00000024195180850,
        -0.00000003193341270, 0.00000000434545060,-0.00000000060578480, 0.00000000008612100,
        -0.00000000001244330, 0.00000000000182260,-0.00000000000027010, 0.00000000000004040,
        -0.00000000000000610, 0.00000000000000090,-0.00000000000000010]

    if(x > 1.000000000010):                                    
        print('problems in LI2')
        print('x=',x)                                                                                          
        print('set x == 1')
        x = 1.0                              
    elif(x > ONE):                                          
      x = 1.0                                                      

    if(x == ONE):                                                
        LI2OLD=PI6
        dli2=LI2OLD                                                            
    elif(x == MONE):                                          
        LI2OLD=MALF*PI6
        dli2=LI2OLD                                                                                                               
    T=-x                                                               
    if(T <= MTWO):                                               
        Y=MONE/(ONE+T)                                                    
        S=ONE                                                             
        A=-PI3+HALF*(np.log(-T)**2-set_precision(np.log(ONE+ONE/T)**2))                        
    elif(T < MONE):                                          
        Y=MONE-T                                                          
        S=MONE                                                            
        A=LOG(-T)                                                         
        A=-PI6+A*(A+set_precision(np.log(ONE+ONE/T)))                                       
    elif(T <= MALF):                                          
        Y=(MONE-T)/T                                                      
        S=ONE                                                             
        A=np.log(-T)                                                         
        A=-PI6+A*(MALF*A+set_precision(np.log(ONE+T)))                                      
    elif(T <= ZERO):                                          
        Y=-T/(ONE+T)                                                      
        S=MONE                                                            
        A=HALF*set_precision(np.log(ONE+T)**2)                                              
    elif(T <= ONE):                                           
        Y=T                                                               
        S=ONE                                                             
        A=ZERO                                                            
    else:                                                               
        Y=ONE/T                                                           
        S=MONE                                                            
        A=PI6+HALF*set_precision(np.log(T)**2)                                                                                                           
                                                                       
    H=Y+Y-ONE                                                          
    ALFA=H+H                                                           
    B1=ZERO                                                            
    B2=ZERO                                                            
    for i in range(18,0,-1):                                                    
      B0=c[i]+ALFA*B1-B2                                               
      B2=B1                                                            
      B1=B0                                                                                                                          
    LI2OLD=-(S*(B0-H*B2)+A) 
    dli2=LI2OLD
    return dli2
"""