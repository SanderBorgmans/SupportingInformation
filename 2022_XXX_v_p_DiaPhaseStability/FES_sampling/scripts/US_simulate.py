
import numpy as np, h5py, os, yaml, sys
import matplotlib.pyplot as pt
import dask, time

import warnings

from molmod.units import *
from ndfsampler import *

def simulate_one(grid_nr,nr,cvs,kappas,input,custom_cv):
    print("{}-{}".format(grid_nr,nr))
    sim = NDFS_Simulation(grid_nr,nr,cvs,kappas,input=input,custom_cv=custom_cv)
    sim.simulate()

def run(grid_name='run.txt'):
    # Simulate on grid - ndfsampler.sim
    grid = np.genfromtxt(grid_name,skip_header=1,delimiter=',',dtype=np.str, encoding='utf-8')
    if grid.size==0:
        print("Done - {} - {} - {}".format(CONFINEMENT_THR, OVERLAP_THR, KAPPA_GROWTH_FACTOR))
        return
    if len(grid.shape)==1:
        grid = np.array([grid])


    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    sims = []

    for n,_ in enumerate(grid):
        gridnr = int(grid[n,0])
        nr     = int(grid[n,1])
        cvs    = [float(v) for v in str(grid[n,2]).split('*')]
        kappas = [float(v) for v in str(grid[n,3]).split('*')]

        simulate = dask.delayed(simulate_one)
        sims.append(simulate(gridnr,nr,cvs,kappas,input=data,custom_cv='./custom_cv.py'))

    sims = dask.compute(*sims)  # trigger computation

############################

if __name__ == '__main__':
    from dask.distributed import Client, progress
    init=time.time()
    client = Client() # take max processes, threads and mem

    run()
    print('This took {} seconds.'.format(time.time()-init))
