# Script to generate traj_X_Y.h5 files for umbrella sampling simulations of Massimo from 2023_NatComm_v14_p1008_MassimoBocus

import numpy as np
import h5py, os
import glob

# start with ai example
# select one of the hops
temp = '873K'
hop = '2-3'
run_number=1
trajectories = glob.glob('/user/gent/419/vsc41947/rdm/paper/2023_NatComm_v14_p1008_MassimoBocus/data/midterm/MLP_simulations/free_energy_profiles/classical/run{}/{}/{}/U*'.format(run_number,temp,hop))


layer_lines = []

for t in trajectories:
    number = int(t.split('/')[-1][1:]) - 1
    identity = (0,number)

    data = np.loadtxt(os.path.join(t,'COLVAR_subsampled'))
    with h5py.File('trajs/traj_{}_{}.h5'.format(*identity),mode='w') as f:
        f['trajectory/cv_values'] = data[:,1]

    with open(os.path.join(t,'plumed_new.dat'),'r') as x:
        lines = x.readlines()

    for line in lines:
        if 'AT1' in line:
            data_line = line.split()
            cv = data_line[1].split('=')[-1]
            kappa = data_line[2].split('=')[-1]
    layer_lines.append(tuple([0,number,cv,kappa,'new_node']))

layer_lines = sorted(layer_lines,key=lambda x:x[1])

with open('layer00.txt','w') as g:
    g.write('layer,nr,cvs,kappas,type\n')
    for line in layer_lines:
        g.write('{},{},{},{},{}\n'.format(*line))

    