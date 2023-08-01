from molmod.units import *
from ogre import *
import sys

def setup_inp(CONFINEMENT_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS,init_kappa):
    # Generate grid - INPUT PARAMETERS
    # mode parameter
    mode     = 'application'

    # grid parameters
    cof2d    = True
    kappas   = [init_kappa,init_kappa]  # in kJ/mol/unit**2
    spacings = [1.5,1.5]   # in unit
    cv_units = ['angstrom','angstrom']

    # MD parameters
    runup    = 10000
    mdsteps  = 200000
    h5steps = 20

    timestep = 0.5*femtosecond

    temp = 300*kelvin
    timecon_thermo = 100.*femtosecond

    press = 1*bar
    timecon_baro = 1000.*femtosecond
    
    # CONSTANTS
    CONFINEMENT_THR = CONFINEMENT_THR
    OVERLAP_THR = OVERLAP_THR
    KAPPA_GROWTH_FACTOR = KAPPA_GROWTH_FACTOR
    MAX_LAYERS = MAX_LAYERS
    HISTOGRAM_BIN_WIDTHS = HISTOGRAM_BIN_WIDTHS
    CONSISTENCY_THR=0.9

    print("Start {} - {} - {}".format(CONFINEMENT_THR, OVERLAP_THR, KAPPA_GROWTH_FACTOR))

    inp = OGRe_Input(mode,kappas=kappas,spacings=spacings,cof2d=cof2d,plot=False,cv_units=cv_units,
                     runup=runup,mdsteps=mdsteps,h5steps=h5steps,timestep=timestep,temp=temp,timecon_thermo=timecon_thermo,press=press,timceon_baro=timecon_baro,
                     CONFINEMENT_THR=CONFINEMENT_THR,OVERLAP_THR=OVERLAP_THR,KAPPA_GROWTH_FACTOR=KAPPA_GROWTH_FACTOR,
                     MAX_LAYERS=MAX_LAYERS,HISTOGRAM_BIN_WIDTHS=HISTOGRAM_BIN_WIDTHS,CONSISTENCY_THR=CONSISTENCY_THR)


    # Generate grid - ogre.input
    inp.make_grid() # this make the initial layer00.txt file

    return inp

############################

if __name__ == '__main__':
    # Load input variables
    MAX_LAYERS = 2
    CONFINEMENT_THR = 0.33
    OVERLAP_THR = 0.50
    KAPPA_GROWTH_FACTOR = 2.0
    init_kappa = 2.5

    HISTOGRAM_BIN_WIDTHS = [0.05, 0.05] # manually fix this

    inp = setup_inp(CONFINEMENT_THR,OVERLAP_THR,KAPPA_GROWTH_FACTOR,MAX_LAYERS,HISTOGRAM_BIN_WIDTHS,init_kappa)
