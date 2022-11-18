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
tqdm.pandas()
from dask.diagnostics import ProgressBar
# plots
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
# stats
from scipy import stats
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
    required_args.add_argument('--gtf', dest='gtfPath', required=True, help='the gtf which maps transcriptome to genome')
    required_args.add_argument('--output', dest='outputDir', required=True, help='where the output files are to be stored')
    ### optional
    optional_args.add_argument('--verbose', dest='verbose', default=0, type=int, help='verbosity, 0 is silent, 1 is loud')
    ################
    parser._action_groups.append(optional_args)
    parser.set_defaults()
    return parser.parse_args(argv[1:])

def importGtf(path):
    ###### gtf #######
    num_lines = sum(1 for line in open(path, 'r'))
    rows = []
    with open(path, 'r') as f:
        for idx, line in enumerate(tqdm(f, total=num_lines)):
            # skip header and only get exons
            if idx > 4:
                if line.split('\t')[2] == "exon":
                    start = line.split('\t')[3]
                    end = line.split('\t')[4]
                    strand = line.split('\t')[6]
                    fields = line.split('\t')[8]
                    for f in fields.split("; "):
                        if f.split(' ')[0] == "gene_id":
                            gene = f.split(' ')[1].strip('\"')
                        if f.split(' ')[0] == "transcript_id":
                            transcript = f.split(' ')[1].strip('\"')
                        if f.split(' ')[0] == "exon_number":
                            exon_num = f.split(' ')[1].strip('\"')
                    ## build row
                    row = pd.Series(data={
                        'geneID': gene,
                        'transcriptID': transcript,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'exon_number': exon_num
                    })
                    ## send to list
                    rows.append(row)
    return pd.concat(rows)

def load_nanopolish(path):
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
        'event_length': np.float32,
        'start_idx': int,
        'end_idx': int}
    npdd = dd.read_table(path, 
        delimiter='\t', 
        usecols = col_list, 
        dtype=col_typedict, 
        #low_memory=True,
        blocksize=25e6)
    with ProgressBar():
        npAny = npdd.compute()
    return npAny

def load_nanopolish_summary(path):
    ###### nanopolish summary #######
    # read_index	read_name	fast5_path	model_name	
    # strand	num_eventsnum_steps	num_skips	num_stays	
    # total_duration	shift	scale	drift	var
    col_list = ["read_index", "read_name", "fast5_path"]
    col_typedict = {'read_index': int, 'read_name': str, 'fast5_path': str}
    nps = dd.read_csv(path, delimiter='\t', usecols = col_list, dtype=col_typedict, low_memory=True)
    with ProgressBar():
        npsOut = nps.compute()
    return npsOut

def check_na(x):
    return math.isnan(x[0]) or x[0] in ['nan', np.nan, ... ]

def test_np_row_model(x, model):
    model = model[model['model_kmer'] == x['model_kmer']]
    ### two sample t-test
    # need means and stdev and samplesize
    x_m = x['event_level_mean']
    x_s = x['event_stdv']
    x_l = x['end_idx']-x['start_idx']
    # model
    y_m = model['model_mean']
    y_s = model['model_stdv']
    y_l = x_l
    ### PERFORM TEST
    stat, pval = stats.ttest_ind_from_stats(x_m, x_s, x_l, y_m, y_s, y_l,
        equal_var = False)
    return pd.Series(data={
        't-stat': stat[0], 
        'pval': pval[0]})

def np_stat(x):
    x_m = x['event_level_mean']
    x_s = x['event_stdv']
    x_l = x['end_idx']-x['start_idx']
    return x_m, x_s, x_l

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

def isInPair(x, npAny):
    return any(x.name == npAny.index)

def find_indices_of_substring(fs, ss):
    return [index for index in range(len(fs)) if fs.startswith(ss, index)]

def extractMappingPos(bam, npAny):
    mappingHolder = dict()
    for read in tqdm(bam.fetch(until_eof=True)):
        mappingHolder[str(read.qname)] = str(read.cigarstring)
    return npAny['read_name'].map(mappingHolder)

def extractFast5data(x):
    rn = "read_" + x['read_name']
    f = h5py.File(x['fast5_path'], 'r')
    sequence = str(f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'][...]).split("\\n")[1]
    trace = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Trace'][...]
    move = f[rn]['Analyses']['Basecall_1D_001']['BaseCalled_template']['Move'][...]
    return pd.Series(data={'sequence': sequence, 'f5signal': trace, 'f5segmentation': move})

def parseCIGARsubsetF5(x):
    # 6. CIGAR: CIGAR string. The CIGAR operations are given in the following table (set ‘*’ if unavailable):
    # Op BAM  Description                              Consumes-query  Consumes-reference
    # M 0   alignment match (can be a sequence match or mismatch)   yes yes
    # I 1   insertion to the reference                              yes  no
    # D 2   deletion from the reference                             no  yes
    # N 3   skipped region from the reference                       no  yes
    # S 4   soft clipping (clipped sequences present in SEQ)        yes no
    # H 5   hard clipping (clipped sequences NOT present in SEQ)    no no
    # P 6   padding (silent deletion from padded reference)         no  no
    # = 7   sequence match yes yes
    # X 8   sequence mismatch yes yes
    ########################################################################
    ### Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
    f5data = extractFast5data(x)
    # loop over each cigar op, apply it to the full seq string, produce the alignment string
    strindex = []
    cursor = 0
    for c in Cigar(x['cigar']).items():
        N, Op = c
        if Op == 'M':
            strindex.append(int(cursor))
            # since this is match we just append +1 repeatedly
            for i in range(N):
                last = strindex[-1]
                strindex.append(int(last+1))
        elif Op == 'D':
            # gap in the ref
            # I can ignore this because I just want to know 
            # which parts of the read map to which parts of the ref? 
            # nanopolish already has genomic position so I think I'm ok
            continue
        elif Op == 'I':
            # gap in the read
            cursor += 1
        elif Op == 'S':
            # used in soft clipping, so the whole-ass seq is in bam but we ignore N of it
            cursor = N
        elif Op == 'H':
            # for Hard clipping, so the truncated seq appears in bam
            print("I hope this doesn't print because then it means I don't understand something")
        else: 
            print("oh no, unrecognized cigar op: " + Op)
    ####### apply strindex
    sequence = ''
    signal = []
    for i in strindex:
        sequence = sequence + f5data.sequence[i]
        signal.append(f5data.f5signal[f5data.f5segmentation == 1][i])
    return sequence, np.array(signal)

def generateData():
    # are we using this?
    #gtf = importGtf(args.gtfPath)
    npA = load_nanopolish(args.npA)[0:1000]
    npB = load_nanopolish(args.npB)[0:1000]
    model = pd.read_csv('../data-files/model_kmer.csv')    
    # define index and find shared entries
    npA.index = pd.MultiIndex.from_arrays(npA[['reference_kmer','position']].values.T)
    npB.index = pd.MultiIndex.from_arrays(npB[['reference_kmer','position']].values.T)
    pcolA = npA.progress_apply(lambda x: isInPair(x, npB), axis=1)
    pcolB = npB.progress_apply(lambda x: isInPair(x, npA), axis=1)
    # subset for only mutually shared entries
    npA = npA[pcolA]
    npB = npB[pcolB]
    ## compute statistic tests on the eventalign-row kmer means
    # t stat is pos when A > B and neg when B > A
    npA[['tstat_model','pval_model']] = npA.progress_apply(lambda x: test_np_row_model(x, model), axis=1)
    npB[['tstat_model','pval_model']] = npB.progress_apply(lambda x: test_np_row_model(x, model), axis=1)
    npA[['tstat_comp','pval_comp']] = npA.progress_apply(lambda x: test_np_row_row(x, npB), axis=1)
    npB[['tstat_comp','pval_comp']] = npB.progress_apply(lambda x: test_np_row_row(x, npA), axis=1)  
    return npA, npB

def processRow(x):
    ### 
    data = parseCIGARsubsetF5(x)
    kmer = x.reference_kmer.replace("U","T")
    pairs = {"A":"T", "C":"G", "G":"C","U":"A" }
    seq = ''.join([pairs[b] for b in data.sequence])
    matches = find_indices_of_substring(seq, kmer)

#def main(argv=sys.argv):
if True:
    argv=sys.argv
    ###
    args = getArgs(argv)

    if not os.path.exists(args.outputDir):
        os.mkdir(args.outputDir)

    ############### load data
    npAu, npBu = generateData()
    ## filter out unsuccessful alignment 
    npAu = npAu[npAu['reference_kmer'] == npAu['model_kmer']]
    npBu = npAu[npAu['reference_kmer'] == npAu['model_kmer']]

    #### filter for significantly different from the model AND from each other
    #pcut = 0.01
    #npA = npAu.loc[(npAu['pval_model'] < pcut) & (npAu['pval_comp'] < pcut)]
    #npB = npBu.loc[(npBu['pval_model'] < pcut) & (npBu['pval_comp'] < pcut)]

    ################################################
    ##################### get f5 ###################
    ################################################
    # get readnames
    npAs = load_nanopolish_summary(args.npAsum)
    npBs = load_nanopolish_summary(args.npBsum)
    npA = npAu.merge(npAs, on='read_index', how='outer').dropna()
    npB = npBu.merge(npBs, on='read_index', how='outer').dropna()
    # get mapping pos
    bamA = pysam.AlignmentFile(args.bamPathA, 'rb')
    bamB = pysam.AlignmentFile(args.bamPathB, 'rb')
    mpA = extractMappingPos(bamA, npA)
    mpB = extractMappingPos(bamB, npB)
    npA.insert(loc=len(npA.columns), column='cigar', value=mpA)
    npB.insert(loc=len(npB.columns), column='cigar', value=mpB)


    row = npA.iloc[98]
    seq, sig = parseCIGARsubsetF5(row)
    kmer = row['model_kmer'].replace("T", "U")
    hits = find_indices_of_substring(seq, kmer)

#if __name__=="__main__":
#    main(sys.argv)




