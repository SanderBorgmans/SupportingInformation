import os, tempfile
import numpy as np, h5py, yaml, sys
import time

from types import SimpleNamespace


from molmod.units import *
from ogre import *

def symlink(target, link_name, overwrite=False):
    '''
    Create a symbolic link named link_name pointing to target.
    If link_name exists then FileExistsError is raised, unless overwrite=True.
    When trying to overwrite a directory, IsADirectoryError is raised.
    '''

    if not overwrite:
        os.symlink(target, link_name)
        return

    # os.replace() may fail if files are on different filesystems
    link_dir = os.path.dirname(link_name)

    # Create link to target with temporary filename
    while True:
        temp_link_name = tempfile.mktemp(dir=link_dir)

        # os.* functions mimic as closely as possible system functions
        # The POSIX symlink() returns EEXIST if link_name already exists
        # https://pubs.opengroup.org/onlinepubs/9699919799/functions/symlink.html
        try:
            os.symlink(target, temp_link_name)
            break
        except FileExistsError:
            pass

    # Replace link_name with temp_link_name
    try:
        # Pre-empt os.replace on a directory with a nicer message
        if not os.path.islink(link_name) and os.path.isdir(link_name):
            raise IsADirectoryError(f"Cannot symlink over existing directory: '{link_name}'")
        os.replace(temp_link_name, link_name)
    except:
        if os.path.islink(temp_link_name):
            os.remove(temp_link_name)
        raise

class OGRe_MLPSimulation_database(OGRe_Simulation):
    def __init__(self,layer,nr,input=None,temp=300*kelvin,press=None,mdsteps=5000,timestep=0.5*femtosecond,h5steps=5,timecon_thermo=100*femtosecond,timecon_baro=1000*femtosecond):
        ''' 
            We need to override constructor function since there is no potential/custom_CV
            **Arguments**

            layer
                the layer number

            nr
                the number of the grid point in this layer

            **Optional Arguments**

            input
                yaml dict or OGRe_Input object

        '''
        # If input is yaml dict, convert it to attribute object using SimpleNamespace
        if isinstance(input,dict):
            input = SimpleNamespace(**input)

        # Input parameters
        self.input = input

        # Grid point parameters
        self.layer = layer
        self.nr = nr

        # Load the rest from the grid file
        fname = 'layer{0:0=2d}.txt'.format(self.layer)
        
        grid = np.genfromtxt(fname, delimiter=',',dtype=None,skip_header=1)
        if grid.size==1: # problem with single line files
            grid = np.array([grid])

        point = grid[self.nr]
        assert point[0]==self.layer
        assert point[1]==self.nr

        if isinstance(point[2],float):
            cvs = np.array([point[2]])
        else:
            cvs = np.array([float(p) for p in point[2].decode().split('*')])

        # Define kappas
        if isinstance(point[3],float):
            kappas = np.array([point[3]])
        else:
            kappas = np.array([float(k) for k in point[3].decode().split('*')])

        dtype = point[4].decode()

        self.cvs = cvs
        self.kappas = kappas
        self.type = dtype

        # Simulation parameters
        self.temp = input.temp if hasattr(input, 'temp') else temp
        self.press = input.press if hasattr(input, 'press') else press
        self.timestep = input.timestep if hasattr(input, 'timestep') else timestep
        self.mdsteps = input.mdsteps if hasattr(input, 'mdsteps') else mdsteps
        self.h5steps = input.h5steps if hasattr(input, 'h5steps') else h5steps
        self.timecon_thermo = input.timecon_thermo if hasattr(input, 'timecon_thermo') else timecon_thermo
        self.timecon_baro = input.timecon_baro if hasattr(input, 'timecon_baro') else timecon_baro

        # Evaluate units
        if hasattr(input, 'cv_units'):
            if not isinstance(input.cv_units, list):
                input.cv_units = [input.cv_units]
            self.cv_units = np.array([eval(unit) if unit is not None else 1.0 for unit in input.cv_units])
        else:
            self.cv_units = np.ones(len(cvs))

        if hasattr(input, 'fes_unit'):
            self.fes_unit = eval(input.fes_unit)
        else:
            self.fes_unit = 1

        # Convert all quantities to atomic units
        self.cvs    = np.asarray(self.cvs)*self.cv_units
        self.kappas = np.asarray(self.kappas)*self.fes_unit/self.cv_units**2

    def database_identity(self):
        # Determine kappa index
        kappa_identity = "{}".format(int(np.round(np.log(np.max(self.kappas)/np.max(np.array(self.input.kappas)*self.fes_unit/self.cv_units**2))/
                                         np.log(self.input.KAPPA_GROWTH_FACTOR)))) # this is the base for the logarithm

        # Determine layer index
        layer_identity = "{}".format(self.layer)

        # Determine trajectory number, as if the whole layer was indexed
        mins = [-0.95] # adapted
        maxs = [0.95] # adapted
        cv_idx = []
        for n,cv in enumerate(self.cvs):
            layer_spacing = self.input.spacings[n]/(2**self.layer)
            cv_range = np.arange(mins[n]-self.input.spacings[n],maxs[n]+self.input.spacings[n]+layer_spacing,layer_spacing) # take possible refined virtual nodes into account
            cv_idx.append(str(int(np.argmin(np.abs(cv_range-cv/self.cv_units[n])))))
        number_identity = "_".join(cv_idx)
        return "{}_{}_{}".format(layer_identity,number_identity,kappa_identity)

    def simulate(self,new_run):
        # Select simulation mode
        if self.input.mode in ['analytic']:
            new_run = self.sim_analytic(new_run)
            return new_run
        elif self.input.mode == 'application':
            new_run = self.sim_application(new_run)
            return new_run
        else:
            raise NotImplementedError('Unknown mode')
        
    def sim_application(self,new_run):
        from pathlib import Path

        # Check whether the simulation has already been performed
        traj_identity = "{}_{}".format(int(self.layer),int(self.nr))
        database_identity = self.database_identity()

        database_path = '../database'
        database_path = Path(database_path).expanduser().resolve().absolute().as_posix()
        database_traj_path = os.path.join(database_path,'traj_{}.h5'.format(database_identity))

        trajs_path = 'trajs/'
        trajs_path = Path(trajs_path).expanduser().resolve().absolute().as_posix()
        trajs_traj_path = os.path.join(trajs_path, 'traj_{}.h5'.format(traj_identity))

        # Initialize the input/output files:
        Path(database_path).mkdir(parents=True, exist_ok=True)
        Path(trajs_path).mkdir(parents=True, exist_ok=True)

        sim = False
        if not os.path.exists(database_traj_path):
            sim = True
        else:
            with open('traj_log.txt','a') as f:
                f.write(database_identity)
                f.write('\n')
            try:
                f = h5py.File(database_traj_path,'r')
            except OSError:
                # this means a bad file header, and sim should be repeated
                sim = True
            else:
                sim = False

        if sim:
            print('No simulation found for {} - {}.'.format(traj_identity, database_identity))
            #print("{} {} {} {} {}\n".format(self.layer,self.nr,"*".join(['{:.8f}'.format(p) for p in self.cvs]),"*".join(['{:.8e}'.format(kappa/(self.fes_unit/self.cv_units**2)[n]) for n,kappa in enumerate(self.kappas)]),self.type))
            new_run+="{},{},{},{},{}\n".format(self.layer,self.nr,"*".join(['{:.8f}'.format(p) for p in self.cvs]),"*".join(['{:.8e}'.format(kappa/(self.fes_unit/self.cv_units**2)[n]) for n,kappa in enumerate(self.kappas)]),self.type)
            return new_run
        
        # Create symbolic link to the data, and overwrite existing symlink to trajectory
        symlink(database_traj_path,trajs_traj_path,overwrite=True)
        return new_run


if __name__ == '__main__':
    init=time.time()

    # Remove grid_restart file if it exists
    if os.path.exists('grid_restart.txt'):
        os.remove('grid_restart.txt')

    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    grid = np.genfromtxt('run.txt',skip_header=1,delimiter=',',dtype=np.str, encoding='utf-8')
    if grid.size==0:
        sys.exit(0)
        
    if len(grid.shape)==1:
        grid = np.array([grid])

    new_run = "layer,nr,cvs,kappas,type\n"

    for n,_ in enumerate(grid):
        layer_nr = int(grid[n,0])
        nr       = int(grid[n,1])

        sim = OGRe_MLPSimulation_database(layer_nr,nr,input=data)
        new_run = sim.simulate(new_run)

    if len(new_run)>28:
        with open('grid_restart.txt','w') as f:
            f.write(new_run)
            
    #print('This took {} seconds.'.format(time.time()-init))
