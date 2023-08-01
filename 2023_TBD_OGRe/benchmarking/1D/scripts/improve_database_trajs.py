import numpy as np
import glob,os
import h5py
import scipy
from molmod.constants import *
from molmod.units import *
from yaff import log
log.set_level(log.silent)

cv_units = 1.
fes_unit = eval('kjmol')

def zscore(x,initial_kappa,kgf,kappa_factor):
    mu = x[0]
    std = np.sqrt(boltzmann*300/(initial_kappa*kgf**kappa_factor))
    return np.abs((x-mu))/std


def zscore_traj(traj,init_kappa,kgf,kappa_factor):
    with h5py.File(traj,'r') as f:
        try:
            tr = f['trajectory/cv_values'][:].ravel()
        except KeyError:
            # we did something, so rerun simulation
            max_z_score = 10000
        else:
            #deltas = np.abs(tr[1:]-tr[:-1])
            #np.round(np.max(deltas)/np.std(tr),2))
            max_z_score = np.max(zscore(tr,init_kappa,kgf,kappa_factor))

    return max_z_score

def get_cv(traj,spacing,mins,maxs):
    identity = traj.split('/')[-1]
    layer_nr = int(identity.split('.')[0].split('_')[1])
    cv_nr = int(identity.split('.')[0].split('_')[2])
    layer_spacing = spacing/(2**layer_nr)
    cv_range = np.arange(mins-spacing,
                         maxs+spacing+layer_spacing,
                         layer_spacing) # take possible refined virtual nodes into account
    return cv_range[cv_nr]

def load_potential_file(file_loc):
    import importlib.util
    spec = importlib.util.spec_from_file_location("module",file_loc)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def perform_simulation(traj,cvs,kappas,fname_potential='scripts/potential.py',mdsteps=40000,h5steps=10,timestep=0.5,temp=300.,timecon_thermo=100.,):
    from molmod.log import TimerGroup, ScreenLog
    from yaff import VerletScreenLog
    from ogre.sim.utils_analytic import SimpleHDF5Writer, SimpleVerletIntegrator, SimpleCSVRThermostat

    timer = TimerGroup()
    log = ScreenLog('OGRe', 1.0, '', '', timer)
    log.set_level(log.silent)

    # Run the simulation
    # Define umbrella
    potential = load_potential_file(fname_potential).Potential
    pot = potential(np.array([cvs])*cv_units, np.array([kappas])*fes_unit/cv_units**2)
    f=h5py.File(traj,mode='w')

    hdf=SimpleHDF5Writer(f,step=h5steps)

    # Initialize the thermostat and the screen logger
    thermo = SimpleCSVRThermostat(temp, timecon=timecon_thermo)
    vsl = VerletScreenLog(step=100)

    # Initialize the US simulation
    verlet = SimpleVerletIntegrator(pot, timestep, hooks=[vsl, thermo, hdf], state=[pot.state],timer=timer,log=log)

    # Run the simulation
    verlet.run(mdsteps)


for dir in glob.glob('scans_*_*_*/'):
    init_kappa = float(dir[:-1].split('_')[2])
    init_spacing = float(dir[:-1].split('_')[1])
    kgf = int(dir[:-1].split('_')[3])
    for traj in glob.glob(os.path.join(dir,'database/*.h5')):
        identity = traj.split('/')[-1]
        kappa_factor = int(identity.split('.')[0].split('_')[-1])
        kappas = init_kappa*kgf**kappa_factor
        cvs = get_cv(traj,init_spacing,-10.,10.)
        zsc = zscore_traj(traj,init_kappa,kgf,kappa_factor)
        while zsc> 400.:
            print(traj, zsc)
            perform_simulation(traj,cvs,kappas,mdsteps=40000,h5steps=10,timestep=0.5,temp=300.,timecon_thermo=100.,)
            zsc = zscore_traj(traj,init_kappa,kgf,kappa_factor)
            print(traj, zsc)