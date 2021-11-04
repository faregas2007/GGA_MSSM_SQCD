# include.py
from abc import ABC, abstractmethod

# factory pattern
class coupling(ABC):

    @abstractmethod
    def get_coup(self):
        pass

    @abstractmethod
    def get_json(self):
        pass
"""
class quarkhiggscoup: SM quark higgs couplings to top and bottom
class squarkhiggscoupCMSSM: real MSSM squark higgs couplings to stop/sbottom in mass-egienstate basis.
"""


# composite pattern/for flexible interface and added new features

class renorm(ABC):
    @abstractmethod
    def get_params(self):
        pass

    @abstractmethod
    def get_json(self):
        pass

# each classes follow basic factory desgin pattern
"""
class sbottomrenorm: sbottom renormalization from HM and DS papers
class Abrenorm: Ab renormalization from HM and DS papers
class mbrenorm: mb renormalization from HM and DS papers
class sqmassmatrix: calculate squark masses, mixing angles and the shifted dMLb
class virtuals: include qcd contributions from interpolate results and sqcd part from full mass depdendence.
""" 

class integ(ABC):

    @abstractmethod
    # add new method of integration here.
    def integrate(self):
        pass
    
    @abstractmethod
    def trans(self):
        pass

    @abstractmethod
    def transplus(self):
        pass

    @abstractmethod
    def transdelta(self):
        pass
    
    @abstractmethod
    def jetsub(self):
        pass

    @abstractmethod
    def jet(self):
        pass
    

# each class will follow basic factory desgin pattern
"""
class integrator: contain all phase space of the transformation into unit hyper cube
class gg: contain real contributions from gg channel
class qg: contain real contributions from qg channel
class qq: contain real contributions from qq channel
"""

