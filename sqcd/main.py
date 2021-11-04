from run import *
from auxiliary import *
from cluster import *
from config import *
from utils import *

import csv
import os, sys
import pandas as pd
from argparse import Namespace

param_fp = Path(config_dir, "params.json")
params = Namespace(**load_dict(param_fp))

dia_index = [2,3,5,7,9,11,12,13,14,15]
filename = 'D'+str(dia_index[params.index_in])+'_'+params.name+'.csv'
Richardson = params.richardson

if (Richardson==False):
    #csvwrite(name, index_in, lam=lam1)
    cluster_run(
        index_in = params.index_in,
        lam = params.lam,
        k=params.k,
        input_file = params.name,
        output_file = Path(logs_dir, params.name),
        amp_dir = amp_dir,
        core_dir = core_dir,
        verbose = False,
        richardson=False
    )
else:
    cluster_run(
        index_in = params.index_in,
        lam = params.lam,
        k=params.k,
        input_file =  params.name,
        output_file = Path(logs_dir, params.name),
        amp_dir = amp_dir,
        core_dir = core_dir,
        verbose = False,
        richardson = True)
    
    while True:
        filenamer='D'+str(dia_index[params.index_in]) +'_'+ params.name +'_'+'richardson' + '.csv'
        os.chdir(working_dir)
        if os.path.isfile(filenamer):
            csvwrite_richard(
                working=str(working_dir), 
                name=params.name, 
                index_in=params.index_in, 
                k=params.k)
            break
        else:
            time.sleep(60)
    print(f"Richardson extrapolation on {params.name} is done")
    print(f'The result is stored at {working_dir}/{filename}')