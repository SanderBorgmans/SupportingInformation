
import numpy as np, h5py, os, yaml, sys
import matplotlib.pyplot as pt
import time

from pathlib import Path
from yaff import *

from molmod.units import *
from ndfsampler import *
from ghostatoms import GhostForceField

def write_ghost_pos(system, ghost_indexes):
    '''
    Update M site position of TIP4P water based on position of other atoms in molecule
    Assumes that atoms are always ordered as O-H-H-M

        r_M = r_O + d_OM^rel/2 * [ (1+d02/d01)*r01 + (1+d01/d02)*r02 ]
    '''
    d_om_rel = 0.13194
    for iatom in ghost_indexes:
        # Vector pointing from O to H1
        r01 = system.pos[iatom-1] - system.pos[iatom-3]
        system.cell.mic(r01)
        d01 = np.linalg.norm(r01)
        # Vector pointing from O to H2
        r02 = system.pos[iatom-2] - system.pos[iatom-3]
        system.cell.mic(r02)
        d02 = np.linalg.norm(r02)
        # Set M position
        system.pos[iatom] = system.pos[iatom-3] + 0.5*d_om_rel*((1.0+d02/d01)*r01 + (1.0+d01/d02)*r02)

def write_ghost_gpos(gpos, vtens, system, ghost_indexes):
    d_om_rel = 0.13194
    for iatom in ghost_indexes:
        # Vector pointing from O to H1
        r01 = system.pos[iatom-1] - system.pos[iatom-3]
        system.cell.mic(r01)
        d01 = np.linalg.norm(r01)
        # Vector pointing from O to H2
        r02 = system.pos[iatom-2] - system.pos[iatom-3]
        system.cell.mic(r02)
        d02 = np.linalg.norm(r02)
        # Partial derivatives of M positions
        pdiff_01 = gpos[iatom,:]*(1.0+d02/d01) - r01*np.dot(gpos[iatom,:],(d02/d01/d01/d01*r01-r02/d01/d02))
        pdiff_02 = gpos[iatom,:]*(1.0+d01/d02) - r02*np.dot(gpos[iatom,:],(d01/d02/d02/d02*r02-r01/d02/d01))
        # Apply chain rule
        gpos[iatom-3,:] += gpos[iatom,:]
        gpos[iatom-3,:] -= 0.5*d_om_rel*pdiff_01
        gpos[iatom-3,:] -= 0.5*d_om_rel*pdiff_02
        gpos[iatom-1,:] += 0.5*d_om_rel*pdiff_01
        gpos[iatom-2,:] += 0.5*d_om_rel*pdiff_02
        if vtens is not None:
            r_mo = 0.5*d_om_rel*((1.0+d02/d01)*r01 + (1.0+d01/d02)*r02)
            vtens[:] -= np.outer(gpos[iatom],r_mo)
            vtens[:] += np.outer(0.5*d_om_rel*pdiff_01,r01)
            vtens[:] += np.outer(0.5*d_om_rel*pdiff_02,r02)

class NDFS_GhostSimulation(NDFS_Simulation):
    def sim_application(self):
        self.log.set_level(self.log.medium)
        Path('logs/').mkdir(parents=True, exist_ok=True)

        with open('logs/log_{}_{}.txt'.format(self.grid,self.nr), 'w') as g:
            self.log.set_file(g)

            system = System.from_file('init.chk', log=self.log)

            # Get the indexes of the ghost atoms (M atoms of TIP4P water)
            ghost_indexes = np.array([iatom for iatom in range(system.natom) if system.get_ffatype(iatom)=='TM'])

            # Set the masses, and adapt the ghost masses to 0
            system.set_standard_masses()
            system.masses[ghost_indexes] = 0.
            system.masses = np.array(system.masses,dtype=np.float) # cast back to float array


            # Create a force field object
            pars=[]
            for fn in os.listdir(os.getcwd()):
                if fn.startswith('pars') and fn.endswith('.txt'):
                    pars.append(fn)

            ff = ForceField.generate(system, pars, log=self.log, timer=self.timer, rcut=12*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True, tailcorrections=True)

            # Create CV list
            cv = self.custom_cv.get_cv(ff)
            
            # Define umbrella
            umbrella = ForcePartBias(ff.system, log=self.log, timer=self.timer)
            for n,c in enumerate(cv):
                bias = HarmonicBias(self.kappas[n], self.cvs[n], c)
                umbrella.add_term(bias)

            # Add the umbrella to the force field
            ff.add_part(umbrella)

            
            ghostff = GhostForceField(ff, ghost_indexes, write_ghost_pos, write_ghost_gpos)


            # Initialize the input/output files:
            Path('trajs/').mkdir(parents=True, exist_ok=True)
            f=h5py.File('trajs/traj_{}_{}.h5'.format(int(self.grid),int(self.nr)),mode='w')
            hdf=HDF5Writer(f,step=self.input.h5steps)

            # Initialize the thermostat, barostat and the screen logger
            thermo = NHCThermostat(self.temp, timecon=self.timecon_thermo)
            baro = MTKBarostat(ghostff, self.temp, self.press, timecon=self.timecon_baro, vol_constraint=False)
            TBC = TBCombination(thermo, baro)
            vsl = VerletScreenLog(step=100)

            # Initialize the US simulation
            verlet = VerletIntegrator(ghostff, self.timestep, hooks=[vsl, TBC, hdf], state=[CVStateItem(cv)])

            # Run the simulation
            verlet.run(self.input.mdsteps)

if __name__ == '__main__':
    init=time.time()
    grid_entry = sys.argv[1:]

    grid_nr = int(grid_entry[0])
    nr     = int(grid_entry[1])
    cvs    = [float(v) for v in str(grid_entry[2]).split('*')]
    kappas = [float(v) for v in str(grid_entry[3]).split('*')]

    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    print("{}-{}".format(grid_nr,nr))
    sim = NDFS_GhostSimulation(grid_nr,nr,cvs,kappas,input=data,custom_cv='./custom_cv.py')
    sim.simulate()
        
    print('This took {} seconds.'.format(time.time()-init))
