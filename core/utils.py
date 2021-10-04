import os, sys
import json, argparse
from typing import Dict, List

import numpy as np
import pandas as pd
from argparse import Namespace
from utils import *
from config import *

def load_dict(filepath):
    with open(filepath) as fp:
        d = json.load(fp)
    return d

"""
params_fp = Path(config_dir, "params.json")
params = Namespace(**load_dict(params_fp))
print(params)
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description='gga_mssm_sqcd')
    parser.add_argument('-id', '--index_in', type=lambda x: int(float(x)), metavar='', help='diagram index',)
    parser.add_argument('-lam', '--lamb', type=lambda x: float(x), metavar='', help='lambda squark mass shifted',)

    parser.add_argument('-in', '--input', type=argparse.FileType('r'), metavar='INPUT_FILE', help='input directory',)
    parser.add_argument('-amp', '--amp_dir', type=dir_path, metavar='AMP_DIR', help='amplitude directory',)

    return parser.parse_args()

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentDefaultsHelpFormatter(f'readable_dir:{path} is not a valid path')


def combined_csv(allfiles, filename):
    combined = pd.concat([pd.read_csv(f) for f in allfiles])
    result = combined.sort_values(by=['lam'], ascending=True)
    result.to_csv(filename, index=False, encoding='utf-8-sig')
    return 0