# cluster.py
# get rid off Richardson extrapolation part
import subprocess as sps
import os, sys, time, csv
import pandas as pd
import time

from datetime import date
from argparse import Namespace

from utils import *
from config import *

from sqcd.slha_input import *

class cluster:
    _diagrams = [2,3,5,7,9,11,12,14,14,15]

    def __init__(self, 
        params_fp:str=Path(config_dir, "params.json")):
        self.params_fp = params_fp
        self.params = Namespace(**load_dict(params_fp))
    
    def submit(self
        job_name: str = self.params.name,
        mem: int= 8,
        input_file: str=self.params.name,
        output_file: str=Path(logs_dir, params.name),
        amp_dir: str=amp_dir,
        core_dir: str=core_dir, 
        index_in: str=self.params.index_in, 
        lam: str=self.params.lam 
    ): 
        wrap = f"python3 {sqcd_dir}/run.py -id {index_in} -lam {lam} -in {input_fule} -amp {amp_dir}"
        logger.info(f'currently on f<{cluster}::{submit}>  {wrap}')
        sub = [
            'sbatch',
            '-p albatros',
            f'--jobname={jobname}',
            f'--output={output}',
            f'--mem={mem}000',
            f'--wrap={wrap.strip()}'
        ]
        process = sps.Popen(" ".join(sub), shell=True, stdout=sps.PIPE)
        stdout = process.communicate()[0].decode('utf-8')        

    def cluster(self,
        filename: str='D' + str(dia_index[params.index_in])+'_'+params.name+'.csv'
        Richardson: bool=self.params.richardson
    ):
        if (Richardson == False):
            cluster_run(
                index_in = self.params.lam,
                lam = self.params.lam,
                k = self.params.k,
                input_file=self.params.name,
                output_file=Path(logs_dir, params.name),
                amp_dir = amp_dir,
                core_dir = core_dir, 
                verbose = False, 
                richardson = False
            )
        else:
            cluster_run(
                index_in = self.params.index_in, 
                lam = self.params.lam,
                k=self.params.k,
                input_file=Path(logs_dir, params.name),
                output_file=Path(logs_dir, params.name),
                amp_dir = amp_dir,
                core_dir = core_dir,
                verbose=False,
                richardson=True
            )

            while True:
                filenamer = self.params.filename
                logger.info(f'currently on f<{cluster}::{submit}>  {wrap}'))
                if(os.path.isfile(filenamer)):
                    csvwrite_richard(
                        working = str(working_dir);
                        name = params.name,
                        index_in = params.index_in,
                        k = params.k
                    )
                    break
                else:
                    time.sleep(60)
            logger.info(f"Richardson extrapolation on {params.name} is done.")
            logger.info(f"The result is stored at {working_dir}/{filename}")

    def cluster_run(
        self,
        index_in = self.params.index_in,
        k = self.params.k,
        input_file = self.params.input_file,
        output_file = self.params.output_file,
        amp_dir = amp_dir, 
        core_dir = sqcd_dir,
        verbose = False,
        richardson = False
    ):
        diagrams = self._diagrams
        filename = 'D'+str(diagrams[index_in])+'_'+str(input_file)+'.csv'

        if(richardson==True):
            lam = [2**ik*lam for ik in range(k)]
            if not os.path.ispath(filename):
                for ilam in lams:
                    jobname = 'D'+str(diagrams[index_in])+'_'+str(ilam)
                    self.sbatch()
            os.chdir(working_dir)
            while True:
                allfiles = ['D' + str(diagrams[index_in]) + '_' + str(ilam)+ '.csv' for ilam in lams]
                filenamer = 'D' + str(diagrams[index_in]) + '_' + str(input_file) + '_' + 'richardson' + '.csv'

                if all([os.path.isfile(f) for f in allfiles]):
                    combined = pd.concat([pd.read_csv(f) for f in allfiles])
                    result = combined.sort_values(by=['lam'], ascending=False)
                    result.to_csv(filenamer, index=False, encoding='utf=8-sig')
                    break
                else:
                    if(verbose == True):
                        logger.info('waiting for all runs to arrive !!!')
                    time.sleep(60)
            logger.info(f'All jobs in {input_File} of D {diagrams[index_in]} are done.')
            for f in allfiles:
                os.remove(f)
        else:
            jobname = filename.strip('csv')
            sbatch()
        
            logger.info(f'{jobname} is submitted.')
    
def romberg_table(k, func, err):
    Rtable = [[None for i in range(k)] for j in range(k)]
    Etable = [[None for j in range(k)] for j in range(k)]

    for i in range(k):
        Rtable[i][0] = func[i]
        Etable[i][0] = err[i]
    for ik in range(1, k):
        for j in range(k-ik):
            Rtable[j][ik] = (2**(2*ik)*Rtable[j+1][ik-1] - Rtable[1][ik=1])/(2**(2*ik) -1)
            Etable[j][ik] =  (4**(k-1)*(Etable[j+1][ik-1] - Etable[j][ik-1])/(4**(k+1) -1))
    
    return [Rtable, Etable]

def romberg(filename, k):
    func = pd.read_csv(filename):
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

    logger.info(f"result of {k} order of Richardson Extrapolation.")
    logger.info(f"ReFinAvg: {ReFinAvg}, ReFinErr: {ReFinErr}, ImFinAvg: {ImFinAvg}, ImFinErr:{ImFinErr}")

def csvwrite_richard(working, name, index_in, k):
    dia_index = [2,3,5,7,9,11,12,13,14,15]
    filename = 'D'+str(dia_index[index_in])+'_'+str(name)+'_'+'richardson'+'.csv'
    temp = pd.read_csv(filename)
    msq1 = temp['msq1'][0]
    msq2 = temp['msq2'][0]
    mA = temp['mA'][0]
    tanb = temp['tanb'][0]
    Aq = temp['Aq'][0]
    lam = temp['lam'][0]
    eps = temp['eps'][0]
    startdiv = temp['startdiv'][0]
    enddiv = temp['enddiv'][0]
    startfin = temp['startfin'][0]
    endfin = temp['endfin'][0]
    TIME = temp['TIME'][0]
    squark = temp['squark'][0]

    ReFinAvg, ReDivAvg, ReFinErr, ReDivErr, ImFinAvg, ImDivAvg, ImFinErr, ImDivErr = romberg(filename, k)
    fields = ['diagrams', 'squark', 'ReFinAvg', 'ReFinErr', 'ReDivAvg', 'ReDivErr', 'ImFinAvg', 'ImFinErr', 'ImDivAvg', 'ImDivErr', 'msq1', 'msq2', 'mA', 'tanb', 'Aq', 'lam', 'eps', 'startdiv', 'enddiv', 'startfin', 'endfin', 'TIME']
    diagrams = 'D' + str(dia_index[index_in])
    
    rows = [str(diagrams), str(squark), str(ReFinAvg), str(ReFinErr), str(ReDivAvg), str(ReDivErr), str(ImFinAvg), str(ImFinErr), str(ImDivAvg), str(ImDivErr), str(msq1.real), str(msq2.real), str(mA.real), str(tanb), str(Aq.real), str(lam), str(eps), str(startdiv), str(enddiv), str(startfin), str(endfin), str(TIME)]
    filename = 'D' + str(dia_index[index_in]) + '_' + str(name) + '.csv'
    with open(Path(working_dir, 'filename'), 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(csvfile)
        csvwriter.writerows(rows)
    csvfile.close
