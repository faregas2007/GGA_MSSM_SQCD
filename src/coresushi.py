# coresushi.py
# include sqcd c factor amplitude

# first run renschemes to generate models.json
# then run sigma.py taken input from models.json to generate hadronic cross-section --> if qcd flags is on.
# sqcd parameters then are feeded into csqcd amplitude calcuations
# the results came back as the csqcd.csv datatable.
# then sigma.py taken input from models.json, cqcd.csv and csqcd.csv to evaluate hadronic cross-section.

# each computational stages could be tracked using logging handlers storing at logs directory. 
# for multiple runs with several parameters points --> class coresushi will have attributes
# for scanning relevant parameters (mA, tanb).
# for various plots: generate a full csv file with a format from full json files
# norder, sigma, error, tanb, mA, muRggh, muFggh, muD, alphasgg, deltamb, mst1, mst2, msb1, msb2, thetab, thetat, Ab, At, muSUSY, mb, mbyuk, mbos, mbsb, mt, mc, mcos, gb, gt, gc, dgb, tanbresum, mbsch, thetasch, Abch

from renschemes import *
from sigma import *

class coresushi:
    def __init__(renschemes, sigma):
        self._renschemes = renschemes
        self._sigma = sigma

    def operations(self):
        pass
    
    def scannings(self):
        pass
    
    def get_json(self):
        pass

    def to_file(self):
        pass

if (__name__ == "main"):
    coresushi()


