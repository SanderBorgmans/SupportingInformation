#! /usr/bin/python
import os,time,h5py,glob,warnings,yaml,sys
import numpy as np
from molmod.units import *


def load_compressed_trajs(compressed_fn,number,grids,trajs,dkappas,energies,identities):
    with h5py.File(compressed_fn,'r') as f:
        grids[number]   = np.array(f['grids'])
        trajs[number]   = np.array(f['trajs'])
        dkappas[number] = np.array(f['dkappas'])
        energies[number] = np.array(f['energies'])
        identities[number] = [tuple(identity) for identity in f['identities']]


def load_trajs(gn,data,grids,trajs,dkappas,energies,identities):
    grid = np.genfromtxt(gn, delimiter=',',dtype=None,skip_header=1)
    for point in grid:
        # Define identity
        gnr = int(point[0])
        nr = int(point[1])
        identity = (gnr,nr)

        if 'multi_load' in data:
            multi_load = data['multi_load']
            # In this case some data points may have already been loaded
            if multi_load and identity in identities[gnr]:
                continue # skip this point

        # Define point
        if isinstance(point[2],float):
            cvs = np.array([point[2]])
        else:
            cvs = np.array([float(p) for p in point[2].decode().split('*')])

        # Define kappas
        if isinstance(point[3],float):
            kappas = np.array([point[3]])
        else:
            kappas = np.array([float(k) for k in point[3].decode().split('*')])
        try:
            with h5py.File('trajs/traj_{}_{}.h5'.format(gnr,nr),'r') as f:
                if 'tdcof' in data and data['tdcof']:
                    try:
                        tr = f['trajectory/cv_values'][:].reshape((1,-1,len(kappas)))
                    except ValueError: # backwards compatibility
                        tr = f['trajectory/cv_values'][:,:,:2].reshape((1,-1,len(kappas)))
                else:
                    try:
                        tr = f['trajectory/cv_values'][:].reshape((1,-1,len(kappas)))
                    except KeyError: # backwards compatibility
                        tr = f['trajectory/cv'][:].reshape((1,-1,len(kappas)))


                energy = f['trajectory/epot'][:].reshape((1,-1))

                if not gnr in grids.keys() or not gnr in trajs.keys():
                    grids[gnr]   = cvs
                    trajs[gnr]   = tr
                    dkappas[gnr] = kappas
                    energies[gnr] = energy
                    identities[gnr] = [identity]
                else:
                    grids[gnr]   = np.vstack((grids[gnr], cvs))
                    trajs[gnr]   = np.vstack((trajs[gnr], tr))
                    dkappas[gnr] = np.vstack((dkappas[gnr], kappas))
                    energies[gnr] = np.vstack((energies[gnr], energy))
                    identities[gnr].append(identity)
                    
        except OSError:
            raise ValueError('Could not find one of the required trajectories! Exiting ...')


def compress_grid(data,index):
    """
        Create h5 file to compress trajectories to contain only required information
    """
    # Load all grids and divide in sub grids based on gridnr
    init = time.time()
    grid_names = sorted(glob.glob('grid[0-9][0-9].txt'), key = lambda gn: int(gn.split('.')[-2].split('grid')[-1]))
    numbers = [int(gn.split('.')[-2].split('grid')[-1]) for gn in grid_names]

    grids = {}
    trajs = {}
    dkappas = {}
    identities = {}
    energies = {}

    multi_load = False
    if 'multi_load' in data:
        multi_load = data['multi_load']

    number = numbers[numbers.index(index)]
    gn = grid_names[numbers.index(index)]

    compressed_fn = 'trajs/compressed_{}.h5'.format(number)
    compressed = os.path.exists(compressed_fn)
    if multi_load:
        # Load both the compressed part and new trajectories, this is convenient when expanding grid
        print("Loading both compressed file and trajectory files for grid {}.".format(number))
        assert compressed
        load_compressed_trajs(compressed_fn,number,grids,trajs,dkappas,energies,identities)
        load_trajs(gn,data,grids,trajs,dkappas,energies,identities)
    else:
        if compressed:
            print("Loading compressed file for grid {}.".format(number))
            load_compressed_trajs(compressed_fn,number,grids,trajs,dkappas,energies,identities)
        else:
            print("Loading trajectory files for grid {}.".format(number))
            load_trajs(gn,data,grids,trajs,dkappas,energies,identities)


    assert len(grids.keys())==1

    h5 = h5py.File("trajs/compressed_{}.h5".format(index),mode='w')
    h5['grids'] = grids[index]
    h5['trajs'] = trajs[index]
    h5['dkappas'] = dkappas[index]
    h5['energies'] = energies[index]
    h5['identities'] = identities[index]

    print("Compressing trajectories took {} seconds.".format(time.time()-init))

if __name__ == '__main__':
    index = int(sys.argv[1])

    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

        compress_grid(data,index)
    else:
        raise IOError('No data file found.')
