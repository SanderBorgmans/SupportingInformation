from molmod.units import *
from ogre import *
import sys

def setup_inp(CONFINEMENT_THR,CONSISTENCY_THR,JS_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS,init_kappa,init_spacing):
    potential_py = """#! /usr/bin/python

import numpy as np
import numpy.ma as ma
from molmod.units import kjmol
import matplotlib.pyplot as pt

from ogre.sim.utils_analytic import AbstractPotential

class Potential(AbstractPotential):
    '''Potential class'''

    def internal_compute(self,gpos):
        gpos += (0.4*self.pos*self.pos*self.pos - 4*self.pos + 16 * self.pos * np.exp(-self.pos*self.pos))*kjmol
        return  self.eval(self.pos)*kjmol

    @staticmethod
    def eval(x):
        return (0.1*x*x*x*x - 2*x*x - 8*np.exp(-x*x) + 10)
"""

    with open('potential.py','w') as f:
        f.write(potential_py)

    # Generate grid - INPUT PARAMETERS
    # mode parameter
    mode     = 'analytic'

    # grid parameters
    edges    = {'min':[-5.], 'max':[5.]}
    kappas   = [init_kappa]  # in kJ/mol/unit**2
    spacings = [init_spacing]   # in unit

    # MD parameters
    runup    = 10000
    mdsteps  = 40000
    h5steps = 10
    timestep = 0.5
    temp = 300.
    timecon_thermo = 100.

    # CONSTANTS
    CONFINEMENT_THR = CONFINEMENT_THR
    OVERLAP_THR = OVERLAP_THR
    KAPPA_GROWTH_FACTOR = KAPPA_GROWTH_FACTOR
    MAX_LAYERS = MAX_LAYERS
    MAX_KAPPA = 5e+03
    HISTOGRAM_BIN_WIDTHS = HISTOGRAM_BIN_WIDTHS
    JS_THR = JS_THR
    CONSISTENCY_THR = CONSISTENCY_THR

    print("Start {} - {} - {}".format(CONFINEMENT_THR, OVERLAP_THR, KAPPA_GROWTH_FACTOR))

    inp = OGRe_Input(mode,kappas=kappas,spacings=spacings,edges=edges,plot=False,
                     runup=runup,mdsteps=mdsteps,h5steps=h5steps,timestep=timestep, temp=temp,timecon_thermo=timecon_thermo,
                     CONFINEMENT_THR=CONFINEMENT_THR,CONSISTENCY_THR=CONSISTENCY_THR,JS_THR=JS_THR,OVERLAP_THR=OVERLAP_THR,KAPPA_GROWTH_FACTOR=KAPPA_GROWTH_FACTOR,
                     MAX_LAYERS=MAX_LAYERS,MAX_KAPPA=MAX_KAPPA,HISTOGRAM_BIN_WIDTHS=HISTOGRAM_BIN_WIDTHS)


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

    if CONFINEMENT_THR>0.0:
        JS_THR = 0.05 # FIXED
    else:
        JS_THR = 1.0

    CONSISTENCY_THR=0.98 # manually fix this
    HISTOGRAM_BIN_WIDTHS = [0.02] # manually fix this, twice the smallest standard deviation 

    inp = setup_inp(CONFINEMENT_THR,CONSISTENCY_THR,JS_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS,init_kappa,init_spacing)
