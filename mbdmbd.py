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
    required_args.add_argument('--output', dest='outDir', required=True, help='where the output files are to be stored')
    # optional
    optional_args.add_argument('--verbose', dest='verbose', action='store_true', help='call this to make program loud')
    optional_args.add_argument('--targetbase', dest='targetBase', default='A', help='you can change which base is in the middle of NN.NN kmer motifs that we compare, default A')
    optional_args.add_argument('--test', dest='test', action='store_true', help='run on only one fast5 as a test')
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

#def main(argv=sys.argv):
if True:
    argv=sys.argv
    ###
    args = get_args(argv)
    # set bars verbosity 
    tqdm.pandas(disable = not args.verbose)
    ### output storage
    if not os.path.exists(args.outDir):
        os.mkdir(args.outDir)
    ### decide f5
    f5str = "*_100.fast5" if args.test else "*.fast5"
    ########################################
    ############ iterate f5 file ###########
    ########################################
    big_list = []
    for f5file in tqdm(glob.glob(os.path.join(args.f5dirA,f5str))):
        f = h5py.File(f5file, 'r')
        for rn in f.keys():
            seq_all, sig_all = extract_f5_data(f, rn)
            seq_idx = [args.targetBase == b for b in seq_all]
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

    kmer_subset = kmer_table[kmer_table['kmer']=='CGACG']

    X = np.array([x for x in kmer_subset['signal']])
    ### gmm ###
    gmm = GaussianMixture(2, covariance_type='full', random_state=0).fit(X)
    #kmer_subset['predicted_clust'] = gmm.predict(X)
    #
    ### plot for sanity checking ###
    plotdata = pd.DataFrame(data={'signal_mean': np.mean(X, axis=1), 'cluster': gmm.predict(X)})
    plot = sns.violinplot(data = plotdata, x='cluster', y='signal_mean')
    plot.set(ylim=(31.725, 32.0))
    fig = plot.get_figure()
    fig.savefig(os.path.join(args.outDir, "test_CGACG.png"))
    exit()

#if __name__=="__main__":
#    main(sys.argv)
