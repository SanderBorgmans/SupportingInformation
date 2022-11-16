
import os, yaml, sys, h5py
import matplotlib.pyplot as pt
from yaff import *

from pathlib import Path
from molmod.units import *

def load_potential_file(file_loc):
    import importlib.util
    spec = importlib.util.spec_from_file_location("module",file_loc)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def simulate(name,kappas,cvs,fn_chk='init.chk',custom_cv='./custom_cv.py',temp=300*kelvin,press=1e5*pascal,mdsteps=10000,timestep=0.5*femtosecond,h5steps=5,timecon_thermo=100*femtosecond,timecon_baro=1000*femtosecond):
    log.set_level(log.medium)

    # this should be a file name
    if os.path.exists(custom_cv):
        custom_cv = load_potential_file(custom_cv)
    else:
        raise IOError('Could not locate the custom_cv module!')




    system = System.from_file(fn_chk)

    # Create a force field object
    pars=[]
    for fn in os.listdir(os.getcwd()):
        if fn.startswith('pars') and fn.endswith('.txt'):
            pars.append(fn)

    ff = ForceField.generate(system, pars, rcut=12*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True, tailcorrections=True)

    # Create CV list
    cv = custom_cv.get_cv(ff)

    # Define umbrella
    umbrella = ForcePartBias(system)
    for n,c in enumerate(cv):
        bias = HarmonicBias(kappas[n], cvs[n], c)
        umbrella.add_term(bias)

    # Add the umbrella to the force field
    ff.add_part(umbrella)

    # Initialize the input/output files:
    Path('trajs/').mkdir(parents=True, exist_ok=True)
    f=h5py.File('trajs/traj_{}.h5'.format(name),mode='w')
    hdf=HDF5Writer(f,step=h5steps)

    # Initialize the thermostat, barostat and the screen logger
    thermo = NHCThermostat(temp, timecon=timecon_thermo)
    baro = MTKBarostat(ff,temp, press, timecon=timecon_baro, vol_constraint=False)
    TBC = TBCombination(thermo, baro)
    vsl = VerletScreenLog(step=100)

    # Initialize the US simulation
    verlet = VerletIntegrator(ff, timestep, hooks=[vsl, TBC, hdf], state=[CVStateItem(cv)])

    # Run the simulation
    verlet.run(mdsteps)



############################

if __name__ == '__main__':
    name   = sys.argv[1]

    # CVs and k
    cvs    = [float(sys.argv[2]), float(sys.argv[3])]
    kappas = [float(sys.argv[4]), float(sys.argv[5])]

    # IMPORTANT
    # This script should be executed in a folder with
    # - pars*.txt files with the force field definition
    # - init.chk file for the structure definition, it can be renamed as a argument of simulate
    # - custom_cv.py script so the simulator knows what to bias, these can be found per material under FES_sampling

    simulate(name,cvs,kappas,mdsteps=100)
