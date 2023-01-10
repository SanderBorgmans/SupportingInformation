
import numpy as np, h5py, os, yaml, sys
import matplotlib.pyplot as pt

import warnings

from molmod.units import *
from ndfsampler import *


def setup_inp(CONFINEMENT_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_GRID_REFINEMENT_IT):
    # Generate grid - INPUT PARAMETERS
    # mode parameter
    mode     = 'application'

    # grid parameters
    edges    = {'min':[0,0], 'max':[22,22]} # this should adapted per material
    kappas   = [10,10]  # in kJ/mol/unit**2
    spacings = [1,1]   # in unit
    units    = ['angstrom','angstrom']

    # MD parameters
    runup    = 10000
    mdsteps  = 30000
    h5steps  = 10
    timestep = 0.5*femtosecond

    temp = 300*kelvin
    timecon_thermo = 100.*femtosecond

    press = 1*bar
    timecon_baro = 1000.*femtosecond


    # CONSTANTS
    CONFINEMENT_THR = CONFINEMENT_THR
    OVERLAP_THR = OVERLAP_THR
    KAPPA_GROWTH_FACTOR = KAPPA_GROWTH_FACTOR
    MAX_GRID_REFINEMENT_IT = int(MAX_GRID_REFINEMENT_IT)

    print("Start {} - {} - {} - {}".format(CONFINEMENT_THR, OVERLAP_THR, KAPPA_GROWTH_FACTOR, MAX_GRID_REFINEMENT_IT))

    inp = NDFS_Input(mode,kappas,spacings,edges,plot=True,units=units,
                     runup=runup,mdsteps=mdsteps,h5steps=h5steps,timestep=timestep,
                     temp=temp,timecon_thermo=timecon_thermo,press=press,timecon_baro=timecon_baro,
                     CONFINEMENT_THR=CONFINEMENT_THR,OVERLAP_THR=OVERLAP_THR,KAPPA_GROWTH_FACTOR=KAPPA_GROWTH_FACTOR,
                     MAX_GRID_REFINEMENT_IT=MAX_GRID_REFINEMENT_IT)


    # Generate grid - ndfsampler.input
    inp.make_grid() # this make the initial grid00.txt file

############################

if __name__ == '__main__':
    # Load input variables
    CONFINEMENT_THR = float(sys.argv[1])
    OVERLAP_THR = float(sys.argv[2])
    KAPPA_GROWTH_FACTOR = float(sys.argv[3])
    MAX_GRID_REFINEMENT_IT = float(sys.argv[4])

    setup_inp(CONFINEMENT_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_GRID_REFINEMENT_IT)

    from shutil import copyfile
    copyfile('grid00.txt', 'run.txt') # create initial run file
