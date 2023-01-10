#! /usr/bin/python

import numpy as np
import pickle
import numpy.ma as ma
import copy


grid = np.load('grid.npy')
with open('fes.pkl', "rb") as input_file:
    fes = pickle.load(input_file)

# we need to smooth FES at edges for MEPSA to work

print(grid.shape)

with open('fes.txt','w') as f:
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            fes_ij = copy.copy(fes[i,j])
            nn = 1
            while ma.is_masked(fes_ij):
                i_min = np.clip(i-nn,0,grid.shape[0])
                i_max = np.clip(i+nn+1,0,grid.shape[0])

                j_min = np.clip(j-nn,0,grid.shape[0])
                j_max = np.clip(j+nn+1,0,grid.shape[0])

                fes_ij_nn_values = [fij for fi in fes[i_min:i_max,j_min:j_max] for fij in fi if not ma.is_masked(fij)]
                if len(fes_ij_nn_values)>0:
                    fes_ij = np.average(fes_ij_nn_values)* (1+nn)
                    break
                nn+=1

            f.write("{}\t{}\t{}\n".format(*grid[i,j],np.clip(fes_ij,0,500)))
