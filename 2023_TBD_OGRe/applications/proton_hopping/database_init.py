from molmod.units import *
from ogre import *
import sys

def setup_inp(CONFINEMENT_THR,CONSISTENCY_THR,JS_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS):
    # Generate grid - INPUT PARAMETERS
    # mode parameter
    mode     = 'application'

    # grid parameters
    edges    = {'min':[-0.85], 'max':[0.85]}
    kappas   = [1000]  # in kJ/mol/unit**2
    spacings = [0.05]   # in unit

    # MD parameters
    runup    = 10000
    mdsteps  = 200000
    h5steps = 10
    timestep = 0.5*femtosecond
    temp = 873.*kelvin
    timecon_thermo = 100.*femtosecond

    # CONSTANTS
    CONFINEMENT_THR = CONFINEMENT_THR
    OVERLAP_THR = OVERLAP_THR
    KAPPA_GROWTH_FACTOR = KAPPA_GROWTH_FACTOR
    MAX_LAYERS = MAX_LAYERS
    MAX_KAPPA = 2.0e+04
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
    MAX_LAYERS = 1
    CONFINEMENT_THR = 0.33
    OVERLAP_THR = 0.33
    KAPPA_GROWTH_FACTOR = 2

    if CONFINEMENT_THR>0.0:
        JS_THR = 0.05 # FIXED
    else:
        JS_THR = 1.0

    CONSISTENCY_THR=0.96 # manually fix this
    HISTOGRAM_BIN_WIDTHS = [0.01] # manually fix this

    inp = setup_inp(CONFINEMENT_THR,CONSISTENCY_THR,JS_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS)
