import os, tempfile
import numpy as np, h5py, yaml, sys
import time

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

class OGRe_Simulation_database(OGRe_Simulation):
    def database_identity(self):
        # Determine kappa index
        kappa_identity = "{}".format(int(np.round(np.log(np.max(self.kappas)/np.max(np.array(self.input.kappas)*self.fes_unit/self.cv_units**2))/
                                         np.log(self.input.KAPPA_GROWTH_FACTOR)))) # this is the base for the logarithm

        # Determine layer index
        layer_identity = "{}".format(self.layer)

        # Determine trajectory number, as if the whole grid was indexed
        mins = self.input.edges['min']
        maxs = self.input.edges['max']
        cv_idx = []
        for n,cv in enumerate(self.cvs):
            layer_spacing = self.input.spacings[n]/(2**self.layer)
            cv_range = np.arange(mins[n]-self.input.spacings[n],maxs[n]+self.input.spacings[n]+layer_spacing,layer_spacing) # take possible refined virtual nodes into account
            cv_idx.append(str(int(np.argmin(np.abs(cv_range-cv/self.cv_units[n])))))
        number_identity = "_".join(cv_idx)
        return "{}_{}_{}".format(layer_identity,number_identity,kappa_identity)

    def sim_application(self):
        from pathlib import Path
        from yaff import VerletScreenLog, ForceField, System, ForcePartBias, HarmonicBias, NHCThermostat, MTKBarostat, TBCombination, VerletIntegrator, CVStateItem, HDF5Writer

        # Check whether the simulation has already been performed
        traj_identity = "{}_{}".format(int(self.layer),int(self.nr))
        database_identity = self.database_identity()

        database_path = '../database/'
        database_path = Path(database_path).expanduser().resolve().absolute().as_posix()
        database_traj_path = os.path.join(database_path,'traj_{}.h5'.format(database_identity))

        logs_path = 'logs/'
        logs_path = Path(logs_path).expanduser().resolve().absolute().as_posix()
        logs_traj_path = os.path.join(logs_path,'log_{}.txt'.format(database_identity))

        trajs_path = 'trajs/'
        trajs_path = Path(trajs_path).expanduser().resolve().absolute().as_posix()
        trajs_traj_path = os.path.join(trajs_path, 'traj_{}.h5'.format(traj_identity))

        # Initialize the input/output files:
        Path(database_path).mkdir(parents=True, exist_ok=True)
        Path(trajs_path).mkdir(parents=True, exist_ok=True)
        Path(logs_path).mkdir(parents=True, exist_ok=True)

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
            #store every new simulation name
            with open('traj_log.txt','a') as f:
                f.write(database_identity)
                f.write('\n')

            self.log.set_level(self.log.medium)
            with open(logs_traj_path, 'w') as g:
                self.log.set_file(g)

                system = System.from_file('init.chk', log=self.log)

                # Create a force field object
                pars=[]
                for fn in os.listdir(os.getcwd()):
                    if fn.startswith('pars') and fn.endswith('.txt'):
                        pars.append(fn)

                ff = ForceField.generate(system, pars, log=self.log, timer=self.timer, rcut=12*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True, tailcorrections=True)

                # Create CV list
                cv = self.custom_cv.get_cv(ff)

                # Adapt structure
                self.custom_cv.adapt_structure(ff,self.cvs)

                # Define umbrella
                umbrella = ForcePartBias(system, log=self.log, timer=self.timer)
                for n,c in enumerate(cv):
                    bias = HarmonicBias(self.kappas[n], self.cvs[n], c)
                    umbrella.add_term(bias)

                # Add the umbrella to the force field
                ff.add_part(umbrella)

                # Initialize the input/output files:
                f=h5py.File(database_traj_path,mode='w')
                hdf=HDF5Writer(f,step=self.h5steps)

                # Initialize the thermostat, barostat and the screen logger
                thermo = NHCThermostat(self.temp, timecon=self.timecon_thermo)
                baro = MTKBarostat(ff,self.temp, self.press, timecon=self.timecon_baro, vol_constraint=False)
                TBC = TBCombination(thermo, baro)
                vsl = VerletScreenLog(step=100)

                # Initialize the US simulation
                verlet = VerletIntegrator(ff, self.timestep, hooks=[vsl, TBC, hdf], state=[CVStateItem(cv)])

                # Run the simulation
                verlet.run(self.mdsteps)
        
        # Create symbolic link to the data, and overwrite existing symlink to trajectory
        symlink(database_traj_path,trajs_traj_path,overwrite=True)
         

if __name__ == '__main__':
    init=time.time()
    grid_entry = sys.argv[1:]

    layer_nr = int(grid_entry[0])
    nr     = int(grid_entry[1])
   
    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    print("{}-{}".format(layer_nr,nr))
    sim = OGRe_Simulation_database(layer_nr,nr,input=data,custom_cv='./custom_cv.py')
    sim.simulate()
        
    print('This took {} seconds.'.format(time.time()-init))
