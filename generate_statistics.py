#!/usr/bin/env python3 
from tqdm import tqdm
import joblib
import h5py
import glob
import numpy as np
import argparse
import sys
import time
import pandas as pd

def extractF5Data(f, rn):
    sequence = str(f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'][...]).split("\\n")[1]
    trace = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Trace'][...]
    move = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Move'][...]
    # return trace where move == 1, so we only get basecalls present in the seq
    return sequence.replace("U","T"), trace[move==1]







def GetSignalByFilename(file_name: str, target_base="A"):
    """
    Wrapper function that extracts kmer table from fast5 file
    """
    f = h5py.File(file_name, 'r')
    big_list = []
    for rn in f.keys():
        seq_all, sig_all = extractF5Data(f, rn)
        seq_idx = [target_base == b for b in seq_all]
        iw = np.where(seq_idx)[0]
        ## vector of individual kmer seq/sig pairs
        seq = [seq_all[i-2:i+3] for i in iw]
        sig = [sig_all[i-2:i+3] for i in iw]
        # remove truncated
        iw2 = np.where([len(s) == 5 for s in seq])[0]
        # check
        if len(seq) != len(sig):
            print("AVAST AND CURSES, ME VECTORS' LENGTHS ARE ASKEW")
        # make serieses
        big_list.append(pd.DataFrame([pd.Series(data={'signal':sig[i].flatten(),'kmer':seq[i],'read_name':rn,'seqloc':str(iw[i])}) for i in iw2]))
    kmer_table = pd.concat(big_list).reset_index()
    del big_list, kmer_table['index']
    return kmer_table



base_directory = sys.argv[1]

st = time.time()
all_fast5_files = glob.glob(base_directory + "*.fast5")
kmer_tables_list = joblib.Parallel(n_jobs=50, verbose=10)(joblib.delayed(GetSignalByFilename)(file_name) for file_name in all_fast5_files)
kmer_table = pd.concat(kmer_tables_list).reset_index()
print("Finished generating all separate kmer tables")

X = np.array([x for x in kmer_table['signal']])

summary = kmer_table.groupby("kmer").aggregate({"signal": lambda x: np.mean(np.array([x_ for x_ in x]))})
fn = time.time()

delta = fn - st

print('Elapsed time:', delta, 'seconds')
