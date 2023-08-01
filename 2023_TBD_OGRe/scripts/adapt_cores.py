#! /usr/bin/python

import sys
import numpy as np

with open('run.txt','r') as f:
    lines = f.readlines()

num_cores = 128
length = len(lines)-1

if length < num_cores:
    nparts = length
else:
    nparts = num_cores


with open('ref_run.pbs','r') as f:
    lines = f.readlines()

with open('run.pbs','w') as f:
    for line in lines:
            if line.startswith("#PBS -l nodes=1:ppn="):
                f.write("#PBS -l nodes=1:ppn={}\n".format(nparts))
            else:
                f.write(line)

