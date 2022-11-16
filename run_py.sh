#!/usr/bin/bash

npa="./data-files/mouse_antiadar_eventalign.txt"
npb="./data-files/mouse_scramble_eventalign.txt"
gtf="./data-files/Mus_musculus.GRCm39.107.gtf"



python -i pythoncode.py --npa $npa  --npb $npb --gtf $gtf --output "./pyout"
