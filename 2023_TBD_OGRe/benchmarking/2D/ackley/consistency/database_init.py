from molmod.units import *
from ogre import *
import sys

def setup_inp(CONFINEMENT_THR,OVERLAP_THR,CONSISTENCY_THR,JS_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS,init_kappa,init_spacing):
    potential_py = """#! /usr/bin/python

import numpy as np
import numpy.ma as ma
from molmod.units import kjmol
import matplotlib.pyplot as pt

from ogre.sim.utils_analytic import AbstractPotential

# 2d Ackley potential

class Potential(AbstractPotential):
    '''Potential class'''

    def internal_compute(self,gpos):
        if not self.pos[0]**2 + self.pos[1]**2 == 0:
            gpos += (np.array([2.* np.exp(-0.2 * np.sqrt(0.5*(self.pos[0]**2 + self.pos[1]**2))) * self.pos[0]/(np.sqrt(0.5*(self.pos[0]**2 + self.pos[1]**2))) +  np.exp(0.5*(np.cos(2*np.pi*self.pos[0]) + np.cos(2*np.pi*self.pos[1]))) * np.pi * np.sin(2*np.pi*self.pos[0]),
                               2.* np.exp(-0.2 * np.sqrt(0.5*(self.pos[0]**2 + self.pos[1]**2))) * self.pos[1]/(np.sqrt(0.5*(self.pos[0]**2 + self.pos[1]**2))) +  np.exp(0.5*(np.cos(2*np.pi*self.pos[0]) + np.cos(2*np.pi*self.pos[1]))) * np.pi * np.sin(2*np.pi*self.pos[1]) ]))*kjmol
        return self.eval(self.pos)*kjmol

    @staticmethod
    def eval(x):
        return -20.*np.exp(-0.2 * np.sqrt(0.5*(x[0]**2 + x[1]**2))) - np.exp(0.5*(np.cos(2*np.pi*x[0]) + np.cos(2*np.pi*x[1]))) + np.exp(1) + 20
    """

    with open('potential.py','w') as f:
        f.write(potential_py)

    # Generate grid - INPUT PARAMETERS
    # mode parameter
    mode     = 'analytic'

    # grid parameters
    edges    = {'min':[-6.,-6.], 'max':[6.,6.]}
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
    CONSISTENCY_THR = CONSISTENCY_THR
    JS_THR = JS_THR

    print("Start {} - {} - {}".format(CONFINEMENT_THR, OVERLAP_THR, KAPPA_GROWTH_FACTOR))

    inp = OGRe_Input(mode,kappas=kappas,spacings=spacings,edges=edges,plot=False,
                     runup=runup,mdsteps=mdsteps,h5steps=h5steps,timestep=timestep,temp=temp,timecon_thermo=timecon_thermo,
                     CONFINEMENT_THR=CONFINEMENT_THR,OVERLAP_THR=OVERLAP_THR,KAPPA_GROWTH_FACTOR=KAPPA_GROWTH_FACTOR,
                     MAX_LAYERS=MAX_LAYERS,HISTOGRAM_BIN_WIDTHS=HISTOGRAM_BIN_WIDTHS,CONSISTENCY_THR=CONSISTENCY_THR,JS_THR=JS_THR)


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
    HISTOGRAM_BIN_WIDTHS = [0.1,0.1] # manually fix this
    CONSISTENCY_THR = float(sys.argv[7])
    JS_THR = float(sys.argv[8])

    inp = setup_inp(CONFINEMENT_THR,OVERLAP_THR,CONSISTENCY_THR,JS_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS,init_kappa,init_spacing)