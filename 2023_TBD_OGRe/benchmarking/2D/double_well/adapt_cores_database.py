#! /usr/bin/python

import sys
import numpy as np

with open('grid_restart.txt','r') as f:
    lines = f.readlines()

num_cores = 24
length = len(lines)-1

if length < num_cores:
    nparts = length
else:
    nparts = num_cores



SPACING, KAPPA, KGF = sys.argv[1:]

with open('ref_run_database.pbs','r') as f:
    lines = f.readlines()

with open('run_database.pbs','w') as f:
    for line in lines:
            if line.startswith("#PBS -l nodes=1:ppn="):
                f.write("#PBS -l nodes=1:ppn={}\n".format(nparts))
            elif line.startswith('python database_simulate_log_database.py'):
                 f.write('python database_simulate_log_database.py $layer $nr {} {} {}\n'.format(SPACING, KAPPA, KGF))
            else:
                f.write(line)

