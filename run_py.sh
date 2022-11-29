#!/usr/bin/bash

#f5a="/storage/liam/re-basecall/p48_anti-adar_guppy6.3/fast5_pass"
#f5b="/storage/liam/re-basecall/p50_Scramble_Si_guppy6.3/fast5_pass"
f5a="/storage/liam/secret-project/mbdmbd/data-files/f5_antiadar/"
f5b="/storage/liam/secret-project/mbdmbd/data-files/f5_scramble/"

python -i mbdmbd.py --dirA $f5a --dirB $f5b --output "./pyout" --verbose
