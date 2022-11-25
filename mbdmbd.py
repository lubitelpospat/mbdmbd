#!/usr/bin/env python3

import sys
import os
import glob
from argparse import ArgumentParser
# dataframes
import pandas as pd
import numpy as np
from dask import dataframe as dd
from dask import compute
# bars
from tqdm import tqdm
from dask.diagnostics import ProgressBar
# stats
from scipy import stats
from sklearn.mixture import GaussianMixture
# specialized datatypes
import h5py
# plots
from matplotlib import pyplot as plt
import seaborn as sns

def get_args(argv):
    parser = ArgumentParser()
    optional_args = parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    ### req
    required_args.add_argument('--dirA', dest='f5dirA', required=True, help='directory of fast5 files from sample A')
    required_args.add_argument('--dirB', dest='f5dirB', required=True, help='directory of fast5 files from sample B')
    required_args.add_argument('--output', dest='outputDir', required=True, help='where the output files are to be stored')
    # optional
    optional_args.add_argument('--verbose', dest='verbose', action='store_true', help='call this to make program loud')
    ##########
    parser._action_groups.append(optional_args)
    parser.set_defaults()
    return parser.parse_args(argv[1:])

def extract_f5_data(f, rn):
    sequence = str(f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'][...]).split("\\n")[1]
    trace = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Trace'][...]
    move = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Move'][...]
    # return trace where move == 1, so we only get basecalls present in the seq
    return sequence.replace("U","T"), trace[move==1]

# https://towardsdatascience.com/gaussian-mixture-models-explained-6986aaf5a95
# https://github.com/ocontreras309/ML_Notebooks/blob/master/GMM_Implementation.ipynb
### go all the way to the bottom
# https://cmdlinetips.com/2021/03/gaussian-mixture-models-with-scikit-learn-in-python/amp/

#def main(argv=sys.argv):
if True:
    argv=sys.argv
    ###
    args = get_args(argv)
    # set bars verbosity 
    tqdm.pandas(disable = not args.verbose)
    ### output storage
    if not os.path.exists(args.outputDir):
        os.mkdir(args.outputDir)
    ########################################
    ############ iterate f5 file ###########
    ########################################
    big_sig = []
    for f5file in tqdm(glob.glob(os.path.join(args.f5dirA,"*"))):
        f = h5py.File(f5file, 'r')
        for rn in f.keys():
            seq_all, sig_all = extract_f5_data(f, rn)
            seq_idx = ['A' == b for b in seq_all]
            seq = [seq_all[i-2:i+3] for i in np.where(seq_idx)[0]]
            sig = [sig_all[i-2:i+3] for i in np.where(seq_idx)[0]]
            # removing truncated kmers
            sig = [sig[i] for i in np.where([len(s) == 5 for s in sig])[0]]
            big_sig.append(sig) 
    # with single bases:
    # meaning signal axis 1 gives mean per base, axis 0 gives 5 mean voltages for all base
    # with kmer:
    # mean signal axis 1 gives 8 voltage channel mean per kmer
    # mean signal axis 0 gives 5x8 mean for whole set
    # meaning again axis 1 on 1 gives one mean value per kmer
    X = np.mean(np.mean(np.concatenate(big_sig), axis=1), axis=1)
    plot = sns.histplot(X, bins=50)
    fig = plot.get_figure()
    fig.savefig("test.png")

    exit()
    n_clusters = 2
    gmm = GaussianMixture(n_components=n_clusters, max_iter=100).fit(X)
    gmm_scores = gmm.score_samples(X)
    print("mean: {}\n score: {}\n".format(gmm.means_, gmm_scores))

#if __name__=="__main__":
#    main(sys.argv)
