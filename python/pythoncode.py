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

############### args
def getArgs(argv):
    parser = ArgumentParser()
    optional_args = parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    ### required
    required_args.add_argument('--npa', dest='npA', required=True, help='nanopolish table A')
    required_args.add_argument('--npb', dest='npB', required=True, help='nanopolish table B')    
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
        "model_kmer", 
        "event_level_mean", 
        "event_stdv", 
        "event_length"]
    col_typedict = {
        'model_kmer': str, 
        'event_level_mean': np.float32, 
        'event_stdv': np.float32, 
        'event_length': np.float32}
    npdd = dd.read_table(path, 
        delimiter='\t', 
        usecols = col_list, 
        dtype=col_typedict, 
        #low_memory=True,
        blocksize=25e6)
    with ProgressBar():
        npAny = npdd.compute()
    return npAny

def kmer_event_dist(k, npAny, rand):
    vals = npAny[npAny['model_kmer'] == k]
    s = np.mean(vals['event_stdv'])
    x = np.mean(vals['event_length'])
    if math.isnan(x) or x in ['nan', np.nan, ... ]:
        x = .01
    # sample frequency is the number hardcoded here
    l = round(x*3000)
    m = np.mean(vals['event_level_mean'])
    dist = s * rand.standard_normal(l) + m
    return dist
    
def kmer_model_dist(k, npAny, rand):
    vals = npAny[npAny['model_kmer'] == k]
    s = np.mean(vals['model_stdv'])
    l = 30
    m = np.mean(vals['model_mean'])
    dist = s * rand.standard_normal(l) + m
    return dist

def check_na(x):
    return math.isnan(x[0]) or x[0] in ['nan', np.nan, ... ]

def main(argv=sys.argv):
#if True:
#    argv=sys.argv
    ###
    args = getArgs(argv)

    if not os.path.exists(args.outputDir):
        os.mkdir(args.outputDir)
    # are we using this?
    #gtf = importGtf(args.gtfPath)

    npA = load_nanopolish(args.npA)
    npB = load_nanopolish(args.npB)
    
    ## load in model
    model = pd.read_csv('model_kmer.csv')    
    # Create a random number generator with a fixed seed for reproducibility
    rng = np.random.default_rng(19680801)
    # histogram bins
    n_bins = 10
    for k in pd.unique(model['model_kmer']):
        distA = kmer_event_dist(k, npA, rng)
        distB = kmer_event_dist(k, npB, rng)
        model_dist = kmer_model_dist(k, model, rng)
        if not check_na(distA) | check_na(distB):
            fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
            axs[0].hist(distA, bins=n_bins)
            axs[1].hist(distB, bins=n_bins)
            #axs[2].hist(model_dist, bins=n_bins)
            fn = os.path.join(args.outputDir, k+".png")
            print(fn)
            plt.savefig(fn)
            plt.clf()
    
if __name__=="__main__":
    main(sys.argv)
