from molmod.units import *
from ogre import *
import sys

def setup_inp(CONFINEMENT_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS,init_kappa,init_spacing):
    potential_py = """#! /usr/bin/python

import numpy as np
import numpy.ma as ma
from molmod.units import kjmol
import matplotlib.pyplot as pt

from ogre.sim.utils_analytic import AbstractPotential

# 2d simple potential from pyretis
class Potential(AbstractPotential):
    '''Potential class'''
    
    def internal_compute(self,gpos):
        gpos += (np.array([self.deriv_x(self.pos),self.deriv_y(self.pos)]))*kjmol
        return self.eval(self.pos)*kjmol
    

    @staticmethod
    def eval(x):
        return (x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]) - 10 * np.exp(-30 * (x[0]-0.2)*(x[0]-0.2) - 3*(x[1]-0.4)*(x[1]-0.4)) \
                                                           - 10 * np.exp(-30 * (x[0]+0.2)*(x[0]+0.2) - 3*(x[1]+0.4)*(x[1]+0.4))
    @staticmethod
    def deriv_x(x):
        return 2*(x[0]*x[0]+x[1]*x[1])*2*x[0] - 10 * np.exp(-30 * (x[0]-0.2)*(x[0]-0.2) - 3*(x[1]-0.4)*(x[1]-0.4)) * (-30)*2*(x[0]-0.2) \
                                              - 10 * np.exp(-30 * (x[0]+0.2)*(x[0]+0.2) - 3*(x[1]+0.4)*(x[1]+0.4)) * (-30)*2*(x[0]+0.2)

    @staticmethod
    def deriv_y(x):
        return 2*(x[0]*x[0]+x[1]*x[1])*2*x[1] - 10 * np.exp(-30 * (x[0]-0.2)*(x[0]-0.2) - 3*(x[1]-0.4)*(x[1]-0.4)) * (-3)*2*(x[1]-0.4) \
                                              - 10 * np.exp(-30 * (x[0]+0.2)*(x[0]+0.2) - 3*(x[1]+0.4)*(x[1]+0.4)) * (-3)*2*(x[1]+0.4)

"""

    with open('potential.py','w') as f:
        f.write(potential_py)

    # Generate grid - INPUT PARAMETERS
    # mode parameter
    mode     = 'analytic'

    # grid parameters
    edges    = {'min':[-0.5,-1.], 'max':[0.5,1.]}
    kappas   = [init_kappa,init_kappa]  # in kJ/mol/unit**2
    spacings = [init_spacing,init_spacing]   # in unit

    # MD parameters
    runup    = 10000
    mdsteps  = 200000
    h5steps = 20
    timestep = 0.5
    temp = 300.
    timecon_thermo = 100.

    # CONSTANTS
    CONFINEMENT_THR = CONFINEMENT_THR
    OVERLAP_THR = OVERLAP_THR
    KAPPA_GROWTH_FACTOR = KAPPA_GROWTH_FACTOR
    MAX_LAYERS = MAX_LAYERS
    HISTOGRAM_BIN_WIDTHS = HISTOGRAM_BIN_WIDTHS
    CONSISTENCY_THR = 0.96

    MAX_KAPPA = 5e+03

    print("Start {} - {} - {}".format(CONFINEMENT_THR, OVERLAP_THR, KAPPA_GROWTH_FACTOR))

    inp = OGRe_Input(mode,kappas=kappas,spacings=spacings,edges=edges,plot=False,
                     runup=runup,mdsteps=mdsteps,h5steps=h5steps,timestep=timestep, temp=temp,timecon_thermo=timecon_thermo,
                     CONFINEMENT_THR=CONFINEMENT_THR,OVERLAP_THR=OVERLAP_THR,KAPPA_GROWTH_FACTOR=KAPPA_GROWTH_FACTOR,
                     MAX_LAYERS=MAX_LAYERS,HISTOGRAM_BIN_WIDTHS=HISTOGRAM_BIN_WIDTHS,MAX_KAPPA=MAX_KAPPA,CONSISTENCY_THR=CONSISTENCY_THR)


    # Generate grid - ogre.input
    inp.make_grid() # this make the initial layer00.txt file

    # Copy the layer00.txt file to run.txt
    from shutil import copyfile
    copyfile('layer00.txt', 'run.txt') # create initial run file

    return inp

############################

if __name__ == '__main__':
    # Load input variables
    MAX_LAYERS = int(sys.argv[1])
    CONFINEMENT_THR = float(sys.argv[2])
    OVERLAP_THR = float(sys.argv[3])
    init_spacing = float(sys.argv[4])
    init_kappa = float(sys.argv[5])
    KAPPA_GROWTH_FACTOR = float(sys.argv[6])
    HISTOGRAM_BIN_WIDTHS = [0.02,0.02] # manually fix this

    inp = setup_inp(CONFINEMENT_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS,init_kappa,init_spacing)