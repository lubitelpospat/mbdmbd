#!/usr/bin/bash


#Adir="/work/liam/nfruns/ADARknockdowncomparison/rebasecall/xpipor-pipe-pore/output/mouse_antiadar/"
#Bdir="/work/liam/nfruns/ADARknockdowncomparison/rebasecall/xpipor-pipe-pore/output/mouse_scramble/"
Adir="/storage/liam/secret-project/mbdmbd/pipeline/output/mouse_antiadar/"
Bdir="/storage/liam/secret-project/mbdmbd/pipeline/output/mouse_scramble/"
npa=$Adir"np/eventalign.txt"
npas=$Adir"np/summary.txt"
bama=$Adir"bam/genome_index.bam"
npb=$Bdir"np/eventalign.txt"
npbs=$Bdir"np/summary.txt"
bamb=$Bdir"bam/genome_index.bam"

samtools index $bama
samtools index $bamb
python -i pythoncode.py --npa $npa --npas $npas --bama $bama --npb $npb --npbs $npbs --bamb $bamb --output "./pyout" --verbose
