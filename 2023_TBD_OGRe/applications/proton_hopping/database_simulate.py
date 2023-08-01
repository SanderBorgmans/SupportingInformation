import os, tempfile
import numpy as np, h5py, yaml, sys
import time

from types import SimpleNamespace

from molmod.units import *
from ogre import *

from yaff import *
from yaff import ForcePartPlumed

import torch
import schnetpack as spk

from ase import Atoms

from schnetpack.environment import TorchEnvironmentProvider

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


class ML_FF(ForcePart):
    def __init__(self, system, model, device = 'cuda'):
        ForcePart.__init__(self, 'ml_ff', system)
        self.system = system
        self.model = model

        self.at = Atoms(positions=self.system.pos / angstrom, numbers=self.system.numbers, cell=self.system.cell.rvecs / angstrom, pbc = True)
        self.converter = spk.data.AtomsConverter(device = device, environment_provider = TorchEnvironmentProvider(6.0, device))

    def _internal_compute(self, gpos, vtens):
        # First update ASE Atoms class
        self.at.set_positions(self.system.pos / angstrom)
        self.at.set_cell(self.system.cell.rvecs / angstrom)

        inputs = self.converter(self.at)
        pred = self.model(inputs)

        energy = pred['energy'].detach().cpu().numpy()[0,0]
        forces = pred['forces'].detach().cpu().numpy()[0]

        if not gpos is None:
            gpos[:, :] = - forces * electronvolt / angstrom

        return energy * electronvolt


def load_system(fname):
    with open('unit_cell', 'r') as f:
        f.readline()
        rvecs = np.array(f.readline().split()[2:-1], dtype=float).reshape((3, 3)) * angstrom
    s = System.from_file(fname, rvecs = rvecs)
    return s

def write_plumed_file(cvs,kappas,plumed_path,colvar_path):
    plumed_file="""RESTART
UNITS LENGTH=A TIME=fs ENERGY=kj/mol

# Definition of CVs
coord1: COORDINATION GROUPA=109 GROUPB=88 R_0=1.6 NN=4
coord2: COORDINATION GROUPA=109 GROUPB=53 R_0=1.6 NN=4
coord3: COORDINATION GROUPA=109 GROUPB=88 R_0=1.4
coord4: COORDINATION GROUPA=109 GROUPB=53 R_0=1.4
cv: MATHEVAL ARG=coord1,coord2 FUNC=x-y PERIODIC=NO
#walls
w1: MATHEVAL ARG=coord3,coord4 FUNC=x+y PERIODIC=NO
LOWER_WALLS ARG=w1 AT=0.65 KAPPA=10000.0 LABEL=lwall
#PRINT ARG=lwall.bias FILE=lwall

#Umbrella
MOVINGRESTRAINT ...
ARG=cv
STEP0=0 AT0={} KAPPA0={}
STEP1=20000 AT1={} KAPPA1={}
 ... MOVINGRESTRAINT

#Create output
PRINT STRIDE=1 ARG=cv FILE={}

FLUSH STRIDE=1
""".format(cvs[0],kappas[0],cvs[0],kappas[0],colvar_path)
    
    with open(plumed_path,'w') as f:
        f.write(plumed_file)


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

        # Determine trajectory number, as if the whole grid was indexed
        mins = [-0.95]
        maxs = [0.95]
        cv_idx = []
        for n,cv in enumerate(self.cvs):
            layer_spacing = self.input.spacings[n]/(2**self.layer)
            cv_range = np.arange(mins[n]-self.input.spacings[n],maxs[n]+self.input.spacings[n]+layer_spacing,layer_spacing) # take possible refined virtual nodes into account
            cv_idx.append(str(int(np.argmin(np.abs(cv_range-cv/self.cv_units[n])))))
        number_identity = "_".join(cv_idx)
        return "{}_{}_{}".format(layer_identity,number_identity,kappa_identity)
         
    def sim_application(self):
        from pathlib import Path
        
        # Check whether the simulation has already been performed
        traj_identity = "{}_{}".format(int(self.layer),int(self.nr))
        database_identity = self.database_identity()

        database_path = '../database'
        database_path = Path(database_path).expanduser().resolve().absolute().as_posix()
        database_traj_path = os.path.join(database_path,'traj_{}.h5'.format(database_identity))

        start_path = '../start_geometries'
        start_path = Path(start_path).expanduser().resolve().absolute().as_posix()
        start_traj_path = os.path.join(start_path,'start_{}.xyz'.format(traj_identity))

        trajs_path = 'trajs/'
        trajs_path = Path(trajs_path).expanduser().resolve().absolute().as_posix()
        trajs_traj_path = os.path.join(trajs_path, 'traj_{}.h5'.format(traj_identity))

        plumed_path = 'plumed/'
        plumed_path = Path(plumed_path).expanduser().resolve().absolute().as_posix()
        plumed_dat_path = os.path.join(plumed_path, 'plumed_{}.dat'.format(traj_identity))
        plumed_log_path = os.path.join(plumed_path, 'plumed_{}.log'.format(traj_identity))
        colvar_path     = os.path.join(plumed_path, 'colvar_{}'.format(traj_identity))


        # Initialize the input/output files:
        Path(database_path).mkdir(parents=True, exist_ok=True)
        Path(trajs_path).mkdir(parents=True, exist_ok=True)
        Path(plumed_path).mkdir(parents=True, exist_ok=True)

        sim = False
        if not os.path.exists(database_traj_path):
            sim = True
        else:                
            try:
                f = h5py.File(database_traj_path,'r')
            except OSError:
                # this means a bad file header, and sim should be repeated
                sim = True
            else:
                sim = False

        if sim:
            with open('traj_log.txt','a') as f:
                f.write(database_identity)
                f.write('\n')
                
            print('Simulating {}'.format(database_traj_path))
            # Run the simulation

            # Check if plumed files still exist, if they do, remove them
            if os.path.exists(plumed_dat_path):
                os.remove(plumed_dat_path)

            if os.path.exists(plumed_log_path):
                os.remove(plumed_log_path)

            if os.path.exists(colvar_path):
                # this is important, colvar file gets appended to
                os.remove(colvar_path)


            # Load the stored model
            device = 'cpu'
            model = torch.load('best_model', map_location=device).to(device)

            # Load the system
            if not os.path.exists(start_traj_path):
                # if the path does not exist, choose the closest one
                mins = self.input.edges['min'][0]
                maxs = self.input.edges['max'][0]
                layer_spacing = self.input.spacings[0]
                init_range = np.arange(mins,maxs+layer_spacing,layer_spacing)
                identity = "0_{}".format(np.argmin(np.abs(self.cvs[0]-init_range)))
                start_traj_path = os.path.join(start_path,'start_{}.xyz'.format(identity))

            system = load_system(start_traj_path)
            ff = ForceField(system, [ML_FF(system, model, device = device)])

            # Define umbrella
            write_plumed_file(self.cvs,self.kappas/self.fes_unit,plumed_dat_path,colvar_path)

            plumed = ForcePartPlumed(system, timestep=self.timestep, fn=plumed_dat_path, fn_log=plumed_log_path)
            ff.add_part(plumed)

            thermo = NHCThermostat(temp=self.temp)
            sl = VerletScreenLog(step=100)

            f=h5py.File(database_traj_path,mode='w')
            hdf=HDF5Writer(f,step=self.h5steps)

            verlet = VerletIntegrator(ff, self.timestep, hooks =[sl, thermo, hdf, plumed], temp0=self.temp)
            verlet.run(self.mdsteps)

            # Add the colvar data to the trajectory file for safekeeping
            colvar_data = np.loadtxt(colvar_path)
            with h5py.File(database_traj_path,mode='a') as f:
                f['trajectory/cv_values'] = colvar_data[:,1]

        # Create symbolic link to the data, and overwrite existing symlink to trajectory
        symlink(database_traj_path,trajs_traj_path,overwrite=True)


if __name__ == '__main__':
    init = time.time()
    grid_entry = sys.argv[1:]
    layer_nr = int(grid_entry[0])
    nr       = int(grid_entry[1])

    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    print("{}-{}".format(layer_nr,nr))
    sim = OGRe_MLPSimulation_database(layer_nr,nr,input=data)
    sim.simulate()
        
    print('This took {} seconds.'.format(time.time()-init))
