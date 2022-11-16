import numpy as np
from yaff import *
import copy, sys, os

log.set_level(0)

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


def geometric_com(system):
    return np.average(system.pos,axis=0)

def com(system):
    return np.einsum('ij,i->j',system.pos,system.masses)/np.sum(system.masses)

def rotate(system, theta, phi, psi):
    pos = copy.copy(system.pos)
    a = np.array([np.cos(theta)*np.cos(psi) ,-np.cos(phi)*np.sin(psi) + np.sin(phi)*np.sin(theta)*np.cos(psi) , np.sin(phi)*np.sin(psi) + np.cos(phi)*np.sin(theta)*np.cos(psi)])
    b = np.array([np.cos(theta)*np.sin(psi) , np.cos(phi)*np.cos(psi) + np.sin(phi)*np.sin(theta)*np.sin(psi) ,-np.sin(phi)*np.cos(psi) + np.cos(phi)*np.sin(theta)*np.sin(psi)])
    c = np.array([-np.sin(theta)            , np.sin(phi)*np.cos(theta)                                       , np.cos(phi)*np.cos(theta)])
    Q = np.array([a,b,c])

    #transform according to Q
    for i in range(0,len(pos)):
        pos[i] = Q.dot(pos[i])

    tmp = copy.copy(system)
    tmp.pos = pos
    return tmp

def deepcopy_system(system):
    return System(
            numbers = system.numbers,
            pos = copy.copy(system.pos), # deepcopy positions or they will get altered
            scopes=system.scopes,
            ffatypes=system.ffatypes,
            ffatype_ids=system.ffatype_ids,
            bonds=system.bonds,
            rvecs=system.cell.rvecs,
            charges=system.charges,
            radii=system.radii,
            valence_charges=system.valence_charges,
            dipoles=system.dipoles,
            radii2=system.radii2,
            masses=system.masses,
        )

def offset(system,vector):
    tmp = deepcopy_system(system)
    tmp.pos[:] += vector
    return tmp

def dist(system,i,j):
    return np.linalg.norm(system.pos[i]-system.pos[j])

def sane(system, num_molecules, length_molecule, verbose=False):
    dists = []
    for n in range(num_molecules):
        for m in range(n+1,num_molecules):
            dists.append([dist(system,i,j) for i in range(len(system.pos)-length_molecule*(n+1), len(system.pos)-length_molecule*n) 
                                           for j in range(len(system.pos)-length_molecule*(m+1), len(system.pos)-length_molecule*m) ])
    if verbose:
        return np.array(dists)
            
    if np.any(np.array(dists)<2.5*angstrom):
        return False
    else:
        return True

def brute_force_insertion(framework,molecule,num_molecules,max_attempts=50):
    attempt=0
    length_molecule = len(molecule.numbers)
    while True:
        new_framework = deepcopy_system(framework)
        com_framework = geometric_com(new_framework)
        molecule_cluster = None

        for n in range(num_molecules):
            molecule_tmp = copy.copy(molecule)

            molecule_tmp.pos -= geometric_com(molecule_tmp) # center molecule at origin

            # Give random orientation
            molecule_tmp = rotate(molecule_tmp,*(np.random.rand(3)*np.pi - np.pi/2))

            # Set the molecule in the center of the framework
            molecule_tmp.pos += com_framework + np.array([-2,1,0])

            # Give random offset
            max_dist= 3.25 * angstrom
            molecule_tmp = offset(molecule_tmp, (np.random.rand(3)*2 -1)*max_dist)

            if n==0:
                molecule_cluster = copy.copy(molecule_tmp)
            else:
                molecule_cluster = molecule_cluster.merge(molecule_tmp)

        # Each unit cell contains 4 pores
        new_framework = new_framework.merge(molecule_cluster)
        
        if sane(new_framework,num_molecules,length_molecule) or attempt>max_attempts:
            if not sane(new_framework,num_molecules,length_molecule):
                print('Brute force did not converge for {} molecules. Trying iterative MC'.format(num_molecules))
                return new_framework, True
            else:
                print('Found a structure after {} attempts for {} molecules.'.format(attempt, num_molecules))
                return new_framework, False
        attempt+=1

def MC_insertion(framework,molecule,num_molecules,iterations=10):
    # this is the length without ghost atoms
    length_molecule = int(len([iatom for iatom in range(molecule.natom) if not molecule.get_ffatype(iatom)=='TM']))
    # Get the indexes of the ghost atoms (M atoms of TIP4P water)
    ghost_indexes = np.array([iatom for iatom in range(framework.natom) if framework.get_ffatype(iatom)=='TM'])

    # Set the masses, and adapt the ghost masses to 0
    framework.set_standard_masses()
    framework.masses[ghost_indexes] = 0.
    framework.masses = np.array(framework.masses,dtype=float) # cast back to float array

    # Create a force field object
    pars=[]
    for fn in os.listdir(os.getcwd()):
        if fn.startswith('pars') and fn.endswith('.txt'):
            pars.append(fn)

    ff = ForceField.generate(framework, pars, rcut=12*angstrom, alpha_scale=3.2, gcut_scale=1.5, smooth_ei=True, tailcorrections=True)
    ghostff = GhostForceField(ff, ghost_indexes, write_ghost_pos, write_ghost_gpos)
    init_energy = ghostff.compute()

    for i in range(iterations):
        if (i%100)==0: print("{}/{}".format(i,iterations))
        for n in range(1,num_molecules): # dont move the first atom, use it as reference
            old_pos = copy.copy(ghostff.system.pos)
            new_pos = copy.copy(ghostff.system.pos)

            # Alter position of molecule n
            new_pos[ghostff.system.natom-(num_molecules-(n))*length_molecule:ghostff.system.natom-(num_molecules-(n+1))*length_molecule] += (np.random.rand(3)*2 -1)
            ghostff.update_pos(new_pos)
            new_energy = ghostff.compute()

            delta_e = new_energy - init_energy

            if delta_e<0 or np.random.uniform(0,1)<np.exp(-delta_e/(boltzmann*300*kelvin)):
                init_energy = new_energy
                ghostff.ff.system.to_file('update.chk')
            else:
                ghostff.update_pos(old_pos)

    return ghostff.ff.system

            

def insert_molecules(framework,molecule,num_molecules,max_attempts=50,max_iterations=10):
    filled_framework,faulty = brute_force_insertion(framework,molecule,num_molecules,max_attempts=max_attempts)
    if faulty:
        filled_framework = MC_insertion(filled_framework,molecule,num_molecules,iterations=max_iterations)

    # Extract molecular cluster
    molecule_cluster = filled_framework.subsystem(np.arange(framework.natom,framework.natom+num_molecules*molecule.natom))
    # Expand over all pores
    # Add another cluster to the right, bottom, and bottom right
    filled_framework = filled_framework.merge(offset(molecule_cluster, framework.cell.rvecs[0]/2))
    filled_framework = filled_framework.merge(offset(molecule_cluster, framework.cell.rvecs[1]/2))
    filled_framework = filled_framework.merge(offset(molecule_cluster, (framework.cell.rvecs[0]+framework.cell.rvecs[1])/2))
    filled_framework.to_file('init.chk')

######################
if __name__=='__main__':
    num_molecules = int(sys.argv[1])
    framework = System.from_file('framework.chk')
    framework.set_standard_masses()
    molecule = System.from_file('H2O.chk')
    insert_molecules(framework,molecule,num_molecules,max_attempts=50,max_iterations=10)        

