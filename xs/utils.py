import csv
import pandas as pd
from getdata import *

def writefile(filename, params, sigmalo, sigmanlo, error, norder, tanbresum, tantresum, mG):
    fields = ['who', 'norder', 'tanbresum', 'tantresum', 'sigmaLO', 'sigmaNLO', 'errorNLO', 'mu', 'tanb', 'mA', 'mG', 'Ab', 'At', 'mst1', 'mst2', 'msb1', 'msb2', 'alphaqsb', 'alphaqst']
    who = 'KIT'
    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        for idx in range(len(params.mA)):
            if norder == 0:
                rows = [[who, norder, tanbresum, tantresum, str(sigmalo[idx]), 0.0, 0.0, 0.0, str(params.mu[idx]), str(params.tanb[idx]), str(params.mA[idx]), str(mG), str(param.Ab[idx]), str(params.At[idx]), str(params.mst1[idx]), str(params.mst2[idx]), str(params.msb1[idx]), str(params.msb2[idx]), str(params.alphaqsb[idx]), str(params.alphaqst[idx])]]
            elif norder == 1:
                rows = [[who, norder, tanbresum, tantresum, str(sigmalo[idx]), 0.0, str(sigmanlo[idx]), str(error[idx]), str(params.mu[idx]), str(params.tanb[idx]), str(params.mA[idx]), str(mG), str(params.Ab[idx]), str(params.At[idx]), str(params.mst1[idx]), str(params.mst2[idx]), str(params.msb1[idx]), str(params.msb2[idx]), str(params.alphaqsb[idx]), str(params.alphaqst[idx])]]
            csvwriter.writerows(rows)   
    csvfile.close