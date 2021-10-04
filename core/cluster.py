import subprocess as sps
import os, sys, time, csv
import pandas as pd
import time 
from datetime import date
from argparse import Namespace

from config import *
from utils import *

param_fp = Path(config_dir, "params.json")
params = Namespace(**load_dict(param_fp))

def sbatch(
    job_name = params.name,
    mem = 8,
    input_file =params.name,
    output_file = Path(logs_dir, params.name),
    amp_dir=amp_dir,
    core_dir=core_dir,
    index_in=params.index_in,
    lam=params.lam):

    os.chdir(input_dir)
    wrap = f"python3 {core_dir}/run.py -id {index_in} -lam {lam} -in {input_file} -amp {amp_dir}"
    print(wrap)
    sub = [
        'sbatch',
        '-p albatros',
        f'--job-name="{job_name}"',
        f'--output="{output_file}"',
        f'--mem={mem}000',
        f'--wrap="{wrap.strip()}"']
    
    process = sps.Popen(" ".join(sub), shell=True, stdout=sps.PIPE)
    stdout = process.communicate()[0].decode('utf-8')

def cluster_run(
    index_in=params.index_in, 
    lam = params.lam,
    k = params.k,
    input_file = params.name,
    output_file=Path(logs_dir, params.name),
    amp_dir=amp_dir,
    core_dir=base_dir, 
    verbose=False, 
    richardson=False):
    diagrams = [2,3,5,7,9,11,12,13,14,15]
    filename = 'D'+str(diagrams[index_in]) +'_'+ str(input_file) + '.csv'
    if (richardson==True):
        #lams = [lam, 2*lam, 4*lam, 8*lam, 16*lam]
        lams = [2**ik*lam for ik in range(k)]
        # save all of the different lam files into D3_richardson directory.
        # save combined lam csv into working directory
        if not os.path.isfile(filename):
            for ilam in lams:
                print(ilam)
                job_name = 'D' + str(diagrams[index_in])+'_'+str(ilam)
                sbatch(
                    job_name=job_name, 
                    mem=8, 
                    input_file=input_file,
                    output_file=output_file, 
                    amp_dir=amp_dir, 
                    core_dir=core_dir, 
                    index_in=index_in, 
                    lam=ilam)
        os.chdir(working_dir)
        while True:
            allfiles = ['D'+str(diagrams[index_in]) + '_' + str(ilam)+'.csv' for ilam in lams]
            filenamer = 'D'+str(diagrams[index_in]) +'_'+ str(input_file)+'_'+'richardson' + '.csv'
            
            if all([os.path.isfile(f) for f in allfiles]):
                combined = pd.concat([pd.read_csv(f) for f in allfiles])
                result = combined.sort_values(by=['lam'], ascending=False)
                result.to_csv(filenamer, index=False,  encoding='utf-8-sig')
                break
            else:
                if(verbose==True):
                    print('waiting for all runs to arrive!!!')
                time.sleep(60)
        print(f'All jobs in {input_file} of D {diagrams[index_in]} are done.')
        # clean csv file
        for f in allfiles:
            os.remove(f)
    
    else:
        job_name = filename.strip('.csv')
        sbatch(
            job_name=job_name,
            mem=8,
            input_file=input_file,
            output_file=output_file,
            amp_dir=amp_dir,
            core_dir=core_dir,
            index_in=index_in,
            lam=lam
        )
        print(input_file)
        print(f'Job {job_name} is submitted.')



def romberg_table(k, func, err):
    Rtable = [[None for i in range(k)] for j in range(k)]
    Etable = [[None for i in range(k)] for j in range(k)]
    for i in range(k):
        Rtable[i][0] = func[i]
        Etable[i][0] = err[i]
    for ik in range(1, k):
        for j in range(k-ik):
            Rtable[j][ik] = (2**(2*ik)*Rtable[j+1][ik-1] - Rtable[j][ik-1])/(2**(2*ik)-1)
            Etable[j][ik] = (4**(k-1)*(Etable[j+1][ik-1] - Etable[j][ik-1])/(4**(k+1) - 1))        
    
    return [Rtable, Etable]
    

def romberg(filename, k):
    """
    Log out the func['refinavg', 'refinerr', 'imfinavg', 'imfinerr']
    Arguments:
        filename: contains numerical values of target function for Richardson extrapolation
        k: order of Richardson extrapolation
        func: k-dataframe of func used in k-order romberg table

    return:
        k order on the romberg table.
        the refinavg, imfinavg
        refinerr --> sqrt((from ibp error (error propagation))**2 + richardson_error**2)
        same for imfinerr
    """
    func = pd.read_csv(filename)
    re_fin_temp = romberg_table(k, func['ReFinAvg'], func['ReFinErr'])
    im_fin_temp = romberg_table(k, func['ImFinAvg'], func['ImFinErr'])
    re_div_temp = romberg_table(k, func['ReDivAvg'], func['ReDivErr'])
    im_div_temp = romberg_table(k, func['ImDivAvg'], func['ImDivErr'])

    ReFinAvg = re_fin_temp[0][0][k-1]
    ImFinAvg = im_fin_temp[0][0][k-1]
    ReDivAvg = re_div_temp[0][0][k-1]
    ImDivAvg = re_div_temp[0][0][k-1]

    err_refin_rb = re_fin_temp[1][0][k-1] 
    err_imfin_rb = im_fin_temp[1][0][k-1]
    err_rediv_rb = re_div_temp[1][0][k-1]
    err_imdiv_rb = im_div_temp[1][0][k-1]

    ReFinErr = err_refin_rb
    ReDivErr = err_imfin_rb
    ImFinErr = err_rediv_rb
    ImDivErr = err_imdiv_rb

    print(f"result of {k} order of Richardson Extrapolation:")
    print(f"ReFinAvg: {ReFinAvg}, ReFinErr: {ReFinErr}, ImFinAvg: {ImFinAvg}, ImFinErr: {ImFinErr}")

    return ReFinAvg, ReDivAvg, ReFinErr, ReDivErr, ImFinAvg, ImDivAvg, ImFinErr, ImDivErr


def csvwrite_richard(working, name, index_in, k):
    dia_index = [2,3,5,7,9,11,12,13,14,15]
    filename = 'D'+str(dia_index[index_in])+'_'+ str(name)+'_'+'richardson'+ '.csv'
    temp = pd.read_csv(filename)
    msq1 = temp['msq1'][0]
    msq2 = temp['msq2'][0]
    mA =  temp['mA'][0]
    tanb = temp['tanb'][0]
    Aq = temp['Aq'][0]
    lam = temp['lam'][0]
    eps = temp['eps'][0]
    startdiv = temp['startdiv'][0]
    enddiv =  temp['enddiv'][0]
    startfin = temp['startfin'][0]
    endfin = temp['endfin'][0]
    TIME = temp['TIME'][0]
    squark = temp['squark'][0]

    ReFinAvg, ReDivAvg, ReFinErr, ReDivErr, ImFinAvg, ImDivAvg, ImFinErr, ImDivErr = romberg(filename, k)
    fields = ['diagrams', 'squark', 'who', 'ReFinAvg', 'ReFinErr', 'ReDivAvg', 'ReDivErr', 'ImFinAvg', 'ImFinErr', 'ImDivAvg', 'ImDivErr','msq1', 'msq2', 'mA', 'tanb', 'Aq', 'lam', 'eps','startdiv', 'enddiv', 'startfin', 'endfin', 'TIME']
    diagrams ='D'+str(dia_index[index_in])
    who = 'KIT'

    rows = [[str(dia_index[index_in]), str(squark), who, str(ReFinAvg), str(ReFinErr), str(ReDivAvg), str(ReDivErr), str(ImFinAvg), str(ImFinErr), str(ImDivAvg), str(ImDivErr), str(msq1.real), str(msq2.real), str(mA.real), str(tanb), str(Aq.real), str(lam), str(eps), str(startdiv), str(enddiv), str(startfin), str(endfin), str(time)]]
    filename = 'D'+str(dia_index[index_in]) +'_'+ str(name) + '.csv'
    with open(working+'/'+filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)
    csvfile.close
    
