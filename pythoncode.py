#!/usr/bin/env python3

import sys
import os
from argparse import ArgumentParser
import pandas as pd
import numpy as np
from dask import dataframe as dd
from dask import compute
from tqdm import tqdm
tqdm.pandas()

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

#def main(argv=sys.argv):
if True:
    argv=sys.argv
    ###
    args = getArgs(argv)

#if __name__=="__main__":
#    main(sys.argv)
