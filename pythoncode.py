#!/usr/bin/env python3

import sys
import os
import math
from argparse import ArgumentParser
# dataframes
import pandas as pd
import numpy as np
from dask import dataframe as dd
from dask import compute
# bars
from tqdm import tqdm
from dask.diagnostics import ProgressBar
# plots
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
# stats
from scipy import stats
from Levenshtein import distance
# specialized datatypes
import h5py
import pysam
from cigar import Cigar

############### args
def getArgs(argv):
    parser = ArgumentParser()
    optional_args = parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    ### required
    required_args.add_argument('--npa', dest='npA', required=True, help='nanopolish eventalign table A')
    required_args.add_argument('--npas', dest='npAsum', required=True, help='nanopolish summary A')
    required_args.add_argument('--bama', dest='bamPathA', required=True, help='bamfile A for getting mapping position')
    required_args.add_argument('--npb', dest='npB', required=True, help='nanopolish eventalign table B')    
    required_args.add_argument('--npbs', dest='npBsum', required=True, help='nanopolish summary B')
    required_args.add_argument('--bamb', dest='bamPathB', required=True, help='bamfile B for getting mapping position')
    required_args.add_argument('--output', dest='outputDir', required=True, help='where the output files are to be stored')
    ### optional
    optional_args.add_argument('--verbose', dest='verbose', action='store_true', help='call this to make program loud')
    ################
    parser._action_groups.append(optional_args)
    parser.set_defaults()
    return parser.parse_args(argv[1:])

def load_nanopolish(path, verbose):
    ###### nanopolish #######
    ##### cols as follows ###
    # contig	position	reference_kmer	read_index	strand	
    # event_index	event_level_mean	event_stdv	event_length	
    # model_kmer	model_mean  model_stdv  standardized_level  
    # start_idx   end_idx
    #########################
    col_list = [
        "contig",
        "position",
        "strand",
        "reference_kmer",
        "model_kmer", 
        "read_index",
        "event_level_mean", 
        "event_stdv", 
        "event_length",
        "model_mean", 
        "model_stdv", 
        "start_idx",
        "end_idx"]
    col_typedict = {
        'contig': str,
        'position': int,
        'strand': str,
        'reference_kmer': str,
        'model_kmer': str, 
        'read_index': int,
        'event_level_mean': np.float32, 
        'event_stdv': np.float32, 
        'model_mean': np.float32,
        'model_stdv': np.float32,
        'event_length': np.float32,
        'start_idx': int,
        'end_idx': int}
    npdd = dd.read_table(path, 
        delimiter='\t', 
        usecols = col_list, 
        dtype=col_typedict, 
        #low_memory=True,
        blocksize=25e6)
    if verbose:
        print('loading nanopolish')
    with ProgressBar():
        npAny = npdd.compute()
    return npAny

def load_nanopolish_summary(path, verbose):
    ###### nanopolish summary #######
    # read_index	read_name	fast5_path	model_name	
    # strand	num_eventsnum_steps	num_skips	num_stays	
    # total_duration	shift	scale	drift	var
    col_list = ["read_index", "read_name", "fast5_path"]
    col_typedict = {'read_index': int, 'read_name': str, 'fast5_path': str}
    nps = dd.read_csv(path, delimiter='\t', usecols = col_list, dtype=col_typedict, low_memory=True)
    if verbose:
        print('loading nanopolish summary')
    with ProgressBar():
        npsOut = nps.compute()
    return npsOut

def np_stat(x):
    x_m = x['event_level_mean']
    x_s = x['event_stdv']
    x_l = x['end_idx']-x['start_idx']
    return x_m, x_s, x_l

def test_np_row_model(x):
    ### two sample t-test
    # need means and stdev and samplesize
    x_m, x_s, x_l = np_stat(x)
    # model
    y_m = x['model_mean']
    y_s = x['model_stdv']
    y_l = x_l
    ### PERFORM TEST
    stat, pval = stats.ttest_ind_from_stats(x_m, x_s, x_l, y_m, y_s, y_l,
        equal_var = False)
    return pd.Series(data={
        't-stat': stat, 
        'pval': pval})

def test_np_row_row(x, yDF):
    y = yDF[x.name == yDF.index]
    ### two sample t-test
    # need means and stdev and samplesize
    x_m, x_s, x_l = np_stat(x)
    y_m, y_s, y_l = np_stat(y)
    ### PERFORM TEST
    stat, pval = stats.ttest_ind_from_stats(x_m, x_s, x_l, y_m, y_s, y_l,
        equal_var = False)
    return pd.Series(data={
        't-stat': stat[0], 
        'pval': pval[0]})

def find_indices_of_substring(fs, ss):
    return [index for index in range(len(fs)) if fs.startswith(ss, index)]

def is_in_pair(x, npAny):
    return any(x.name == npAny.index)

def extract_bam_data(bam):
    # making pd df for this
    seriesList = []
    for read in tqdm(bam.fetch(until_eof=True)):
        # https://github.com/pysam-developers/pysam/blob/cb3443959ca0a4d93f646c078f31d5966c0b82eb/pysam/libcalignedsegment.pyx#L1845
        read_name = str(read.qname)
        ref_name = str(read.reference_name)
        qend = read.qstart + len(read.query_alignment_sequence)
        series = pd.Series(data={'read_name': read_name, 
                                'contig': ref_name, 
                                'query_start': read.qstart, 
                                'query_end': qend, 
                                'ref_start': read.reference_start,
                                'query_seq': read.query_sequence,
                                'cigar': read.cigarstring})
        seriesList.append(series)
    return pd.concat(seriesList, axis=1).T
    
def extract_f5_data(x):
    rn = "read_" + x['read_name']
    f = h5py.File(x['fast5_path'], 'r')
    sequence = str(f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'][...]).split("\\n")[1]
    trace = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Trace'][...]
    move = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Move'][...]
    # return trace where move == 1, so we only get basecalls present in the seq
    return pd.Series(data={'sequence': sequence.replace("U","T"), 'signal': trace[move==1]})

def interpret_cigar(x):
    # COLUMNS: ["read_name","contig","query_start","query_end","ref_start",
    # "query_seq","cigar","fast5_path","sequence","signal","position","np_seq"]
    cigar = x['cigar']
    strindex = []
    cursor = 0
    for c in Cigar(cigar).items():
        N, Op = c
        #### handle op
        if Op == "M":
            # since this is a match we append +1 repeatedly
            for i in range(N):
                strindex.append(cursor)
                cursor += 1
        elif Op == 'D':
            # deletion in ref, ignore
            continue
        elif Op == 'I':
            # instertion in read, dupe prev
            # might need to N-1
            for i in range(N):
                strindex.append(cursor)
        elif Op == 'S':
            # soft clipping
            cursor += N
        elif Op == 'H':
            # hard clipping, we ignore
            continue
        elif Op == 'N':
            # skip ref region, same as S?
            cursor += N
        else:
            print("yikes we got an unrecognized cigar Op: " + Op)
    # return strindex
    # seq    CCGA
    # np_seq CCGGA
    seq = ''
    for i in range(len(x['np_seq'])):
        ref_i = strindex[i]
        refbase = x['sequence'][ref_i]
        npbase = x['np_seq'][i]
        if refbase == npbase:
            seq += refbase
        else:
            if x['sequence'][strindex[i-1]] == npbase:
                strindex.insert(i-1, strindex[i-1])
            else:
                print("help")
    print(seq)
    print(strindex)
    return strindex

def apply_strindex(x):
    seqH = []
    sig = []
    for i in x['cigind']:
        seqH.append(x['sequence'][i])
        sig.append(x['signal'][i])
    seq = ''.join(seqH)
    return pd.Series(data={ 'sequence': seq, 'signal': sig })

def get_kmer_subset(x):
    # COLUMNS: ["read_name","contig","position","reference_kmer",
    # "query_start","query_end","ref_start","cigar",
    # "fast5_path","sequence","signal"]
    pos = int(x['position']-x['ref_start'])
    kmer = x['sequence'][pos:pos+5]
    sig = x['signal'][pos:pos+5]
    if kmer != x['reference_kmer']:
        print("kmer missmatch, {}, {}, {},  ref: {} target: {}".format(
            x['read_name'],
            x['contig'],
            x['position'],
            x['reference_kmer'],
            kmer))
    return pd.Series(data={ 'kmer': kmer,'sig': sig})

def build_np_seq(x):
    lastpos = max(x['position'])
    # combine first base of kmers
    seq = ''.join([b[0] for b in x['reference_kmer']])
    # add last 4
    lastkmer = str(x[x['position']==lastpos]['reference_kmer'].values[0][1:5])
    return seq + lastkmer

def generate_nanopolish_data(args):
    npA = load_nanopolish(args.npA, args.verbose)[0:2500]
    npB = load_nanopolish(args.npB, args.verbose)[0:2500]
    # define index and find shared entries
    npA.index = pd.MultiIndex.from_arrays(npA[['reference_kmer','position']].values.T)
    npB.index = pd.MultiIndex.from_arrays(npB[['reference_kmer','position']].values.T)
    if args.verbose:
        print('filtering for shared genome locations')
    pcolA = npA.progress_apply(lambda x: is_in_pair(x, npB), axis=1)
    pcolB = npB.progress_apply(lambda x: is_in_pair(x, npA), axis=1)
    # subset for only mutually shared entries
    npA = npA[pcolA]
    npB = npB[pcolB]
    ## compute statistic tests on the eventalign-row kmer means
    # t stat is pos when A > B and neg when B > A
    if args.verbose:
        print('doing row-wise stat test event against model')
    npA[['tstat_model','pval_model']] = npA.progress_apply(lambda x: test_np_row_model(x), axis=1)
    npB[['tstat_model','pval_model']] = npB.progress_apply(lambda x: test_np_row_model(x), axis=1)
    if args.verbose:
        print('doing row-wise stat test event against the other dataframes event')
    npA[['tstat_comp','pval_comp']] = npA.progress_apply(lambda x: test_np_row_row(x, npB), axis=1)
    npB[['tstat_comp','pval_comp']] = npB.progress_apply(lambda x: test_np_row_row(x, npA), axis=1)  
    return npA, npB

#def main(argv=sys.argv):
if True:
    argv=sys.argv
    ###
    args = getArgs(argv)
    # set bars verbosity 
    tqdm.pandas(disable = not args.verbose)
    ### output storage
    if not os.path.exists(args.outputDir):
        os.mkdir(args.outputDir)

    ############### load data
    npAu, npBu = generate_nanopolish_data(args)
    # COLUMNS: ["contig","position","strand","reference_kmer","model_kmer","read_index","event_level_mean",
    # "event_stdv","event_length","model_mean","model_stdv","start_idx","end_idx"]
    ## filter out unsuccessful alignment 
    npAu = npAu[npAu['reference_kmer'] == npAu['model_kmer']]
    npBu = npAu[npAu['reference_kmer'] == npAu['model_kmer']]

    #### filter for significantly different from the model
    pcut = 0.01
    npA = npAu.loc[(npAu['pval_model'] < pcut)]
    npB = npBu.loc[(npBu['pval_model'] < pcut)]
    ################################################
    ########### get sequence metadata ##############
    ################################################
    # get readnames
    npAs = load_nanopolish_summary(args.npAsum, args.verbose)
    npBs = load_nanopolish_summary(args.npBsum, args.verbose)
    npA = npAu.merge(npAs, on='read_index', how='outer').dropna()
    npB = npBu.merge(npBs, on='read_index', how='outer').dropna()
    # COLUMNS: ["contig","position","strand","reference_kmer","model_kmer","read_index","event_level_mean",
    # "event_stdv","event_length","model_mean","model_stdv","start_idx","end_idx", 
    # NEW: "read_name","fast5_path"]
    del npAu, npBu, npAs, npBs
    ################################################
    ############### get alignment data #############
    ################################################
    if args.verbose:
        print('extracting bam alignment data')
    # dicts w/ key: (readname, refname), val: read to ref mapping pos
    bamA = pysam.AlignmentFile(args.bamPathA, 'rb')
    bamB = pysam.AlignmentFile(args.bamPathB, 'rb')
    bamdfA = extract_bam_data(bamA)
    bamdfB = extract_bam_data(bamB)
    # COLUMNS: ["read_name","contig","query_start","query_end","ref_start","query_seq","cigar"]
    ################################################
    ############# get voltage data #################
    ################################################
    # full size and cig only needs to be done once per read
    if args.verbose:
        print('getting sequence-level signal and bases')
    readsetA = npA[['read_name', 'fast5_path']].drop_duplicates()
    readsetB = npB[['read_name', 'fast5_path']].drop_duplicates()
    readsetA[['sequence', 'signal']] = readsetA.progress_apply(lambda x: extract_f5_data(x), axis=1)
    readsetB[['sequence', 'signal']] = readsetB.progress_apply(lambda x: extract_f5_data(x), axis=1)
    # COLUMNS: ["read_name","fast5_path","sequence","signal"]
    ########## merge ##########
    rsA = bamdfA.merge(readsetA, on='read_name', how='outer').dropna()
    #rsB = bamdfB.merge(readsetB, on='read_name', how='outer').dropna()
    # COLUMNS: ["read_name","contig","query_start","query_end","ref_start","query_seq",
    # "cigar","fast5_path","sequence","signal"]
    #del bamA, bamB, bamdfA, bamdfB
    ###############################################
    ############## interpreting cigar #############
    ###############################################
    ########### first get the nanopolish sequence ############
    # rename cols using this one weird trick, doctors hate him
    def f(x):
        try:
            int(x)
            return 'np_seq'
        except ValueError:
            return x
    # afaik all we need is read_name, contig, position, reference_kmer
    seqssA = npA[['read_name','contig','position','reference_kmer']]
    #seqssB = npB[['read_name','contig','position','reference_kmer']]
    npseqA = seqssA.drop_duplicates().groupby(['read_name', 'contig']).progress_apply(lambda x: build_np_seq(x)).reset_index().rename(columns=f)
    #npseqB = seqssB.drop_duplicates().groupby(['read_name', 'contig']).progress_apply(lambda x: build_np_seq(x)).reset_index().rename(columns=f)
    #begub
    rsA = rsA[rsA['read_name']=='2fdb8afd-e3ba-4732-a2d7-9e7ef388e43f']
    rsA = rsA.merge(npseqA, on=['read_name','contig'], how='outer').dropna()

    # add np seqs
    exit()
    cigindA = rsA.progress_apply(lambda x: interpret_cigar(x), axis=1)
    #cigindB = rsB.merge(npseqB, on=['read_name','contig'], how='outer').progress_apply(lambda x: interpret_cigar(x), axis=1)
    rsA.insert(loc=len(rsA.columns), column='cigind', value=cigindA)
    #rsB.insert(loc=len(rsB.columns), column='cigind', value=cigindB)
    rsA[['sequence', 'signal']] = rsA.apply(lambda x: apply_strindex(x), axis=1)
    #rsB[['sequence', 'signal']] = rsB.apply(lambda x: apply_strindex(x), axis=1)
    ########################################################
    ############## merge with np and get kmers #############
    ########################################################
    exit()
    # build kmer seq
    ksA = npA[['read_name','contig','position','reference_kmer']].merge(rsA, on=['read_name','contig'],how='outer').dropna()
    ksA[['kmer', 'sig']] = ksA.progress_apply(lambda x: get_kmer_subset(x), axis=1)
    #ksB = npB[['read_name','contig','position','reference_kmer']].merge(rsB, on=['read_name','contig'],how='outer').dropna()
    #ksB[['kmer', 'sig']] = ksB.progress_apply(lambda x: get_kmer_subset(x), axis=1)
    exit()
    ###################################################################
    ######################### WRITE OUT ###############################
    ###################################################################
    if args.verbose:
        print('writing to file')
    npA.to_csv(os.path(args.outputDir, "Atable.csv"))
    npB.to_csv(os.path(args.outputDir, "Btable.csv"))

#if __name__=="__main__":
#    main(sys.argv)





