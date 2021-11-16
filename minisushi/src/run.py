import os, csv, sys, pyslha

from importlib import import_module
from timeit import default_timer as timer
from argparse import Namespace

from init import *
from config import *
from utils import *

from renschemes.renschemes import *
from xs.sigma import *


params = Namespace(**load_dict(Path(config_dir, 'params.json')))

def main():
    #return sigma(qqinteg_batch, qginteg_batch, gginteg_batch).to_file()
    #return sigma_sqcd(qq, qg, gg).to_file()
    if(params.sqcd == True):
        return sigma_sqcd(qq, qg, gg).xs()
    else:
        return sigma(qq, qg, gg).xs()

def csvwrite(
    working: str=working_dir, 
    name: str=Path(data_dir, 'MSSM_SLHA.in')):
    params = Namespace(**load_dict(Path(model_dir, 'models.json')))
    
    Ab = params.Ab
    At = params.At
    
    msb1 = np.sqrt(params.msb12)
    msb2 = np.sqrt(params.msb22)
    mst1 = np.sqrt(params.mst12)
    mst2 = np.sqrt(params.mst22)

    mA = params.mA
    tanb = params.tanb
    mgl = params.mgl

    sthetat = params.sthetat
    cthetat = params.cthetat
    sthetab = params.sthetab
    cthetab = params.cthetab

    gc = params.gc
    gt = params.gt

    # higgs-sbottom couplings mixing
    gbh11 = params.gbh11
    gbh12 = params.gbh12
    gbh21 = params.gbh21
    gbh22 = params.gbh22

    # higgs-stop couplings mixing
    gth11 = params.gth11
    gth12 = params.gth12
    gth21 = params.gth21
    gth22 = params.gth22

    mb = params.mb
    mt = params.mt

    tmp = main()
    norder = tmp['norder']
    sigmaLO = tmp['sigmaLO']
    errorLO = tmp['errorLO']
    sigmaNLO = tmp['sigmaNLO']
    errorNLO = tmp['errorNLO']

    filename = 'xs_'+str(name).strip('.in')+'.csv'

    fields = ['norder', 'sigmaLO', 'errorLO', 'sigmaNLO', 'errorNLO', 'mb', 'mt', 'mA', 'tanb', 'mgl,' 'sthetat', 'cthetat', 'sthetab', 'cthetab', 'msb1', 'msb2', 'mst1', 'mst2', 'Ab', 'At']
    rows = [[str(norder), str(sigmaLO), str(errorLO), str(sigmaNLO), str(errorNLO), str(mb), str(mt), str(mA), str(tanb), str(mgl), str(sthetat), str(cthetat), str(sthetab), str(cthetab), str(msb1), str(msb2), str(mst1), str(mst2), str(Ab), str(At)]]
    with open(Path(working_dir, str(filename)), 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)
    csvfile.close 

# need to have a seperate method --> array typed input.
def csv_write_sqcd(
    working: str=Path(working_dir), 
    name: str=Path(data_dir, 'MSSM_SLHA.in')):
    tmp = main()
    norder = tmp['norder']
    sigmaLO = tmp['sigmaLO']
    errorLO = tmp['errorLO']
    sigmaNLO = tmp['sigmaNLO']
    errorNLO = tmp['errorNLO']

    #If(Micheal == False):
    #params = Namespace(**load_dict(Path(model_dir, 'models.json')))
    params = getdata(Path(data_dir, 'CSQCD.csv'))
    
    Ab = params.Ab
    At = params.At
    
    msb1 = params.msb1
    msb2 = params.msb2
    mst1 = params.mst1
    mst2 = params.mst2

    mA = params.mA
    tanb = params.tanb
    mu = params.mu

    #sthetat = params.sthetat
    #cthetat = params.cthetat
    #sthetab = params.sthetab
    #cthetab = params.cthetab

    #gc = params.gc
    #gt = params.gt

    # higgs-sbottom couplings mixing
    #gbh11 = params.gbh11
    #gbh12 = params.gbh12
    #gbh21 = params.gbh21
    #gbh22 = params.gbh22

    # higgs-stop couplings mixing
    #gth11 = params.gth11
    #gth12 = params.gth12
    #gth21 = params.gth21
    #gth22 = params.gth22

    #mb = params.mb
    #mt = params.mt

    filename = 'xs_'+str(name).strip('.in') + '.csv'
    fields = ['norder', 'sigmaLO', 'errorLO', 'sigmaNLO', 'errorNLO', 'mu', 'mA', 'tanb', 'msb1', 'msb2', 'mst1', 'mst2', 'Ab', 'At']
    #rows = [[str(norder), str(sigmaLO), str(errorLO), str(sigmaNLO), str(errorNLO), str(mA), str(tanb), str(mG), str(msb1), str(msb2), str(mst1), str(mst2), str(Ab), str(At)]]
    with open(Path(working_dir, str(filename)), 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        for idx in range(len(params.mA)):
            if norder == 0:
                rows = [[norder, str(sigmaLO[0]), str(errorLO[0]), str(sigmaNLO[0]), str(errorNLO[0]), str(params.mu[idx]), str(params.tanb[idx]), str(params.mA[idx]), str(params.Ab[idx]), str(params.At[idx]), str(params.mst1[idx]), str(params.mst2[idx]), str(params.msb1[idx]), str(params.msb2[idx])]]
            elif norder == 1:
                rows = [[norder, str(sigmaLO[0]), str(errorLO[0]), str(sigmaNLO[idx]), str(errorNLO[idx]), str(params.mu[idx]) ,str(params.tanb[idx]), str(params.mA[idx]), str(params.Ab[idx]), str(params.At[idx]), str(params.mst1[idx]), str(params.mst2[idx]), str(params.msb1[idx]), str(params.msb2[idx])]]
            csvwriter.writerows(rows)   
    csvfile.close 

if __name__ == "__main__":
    if(params.sqcd == True):
        csv_write_sqcd(working=working_dir,name='MSSM_SLHA.in')
    else:
        csvwrite(working=working_dir, name='MSSM_SLHA.in')