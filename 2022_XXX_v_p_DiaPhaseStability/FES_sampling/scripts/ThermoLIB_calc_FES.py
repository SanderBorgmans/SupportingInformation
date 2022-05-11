#! /usr/bin/python
import os,sys,h5py,copy,warnings,subprocess,glob,pickle,yaml,time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as pt

from molmod.constants import boltzmann
from molmod.units import *

from thermolib.thermodynamics.fep import FreeEnergySurface2D, plot_feps
from thermolib.thermodynamics.histogram import Histogram2D
from thermolib.tools import read_wham_input_2D


def get_units(data):
    if 'units' in data:
        if not isinstance(data['units'], list):
            units = [data['units']] * len(data['spacings'])
        else:
            units = data['units']
        assert len(units)==len(data['spacings'])

        units = np.array([eval(unit) if unit is not None else 1.0 for unit in units])
    else:
        units = np.ones(len(data['spacings']))
    return units


def load_compressed_trajs(compressed_fn,number,grids,trajs,dkappas):
    with h5py.File(compressed_fn,'r') as f:
        grids[number]   = np.array(f['grids'])
        trajs[number]   = np.array(f['trajs'])
        dkappas[number] = np.array(f['dkappas'])

def load_grid(path,data,verbose=True,force=False):
    """
        Create grid if it does not yet exist
    """
    # Load all grids and divide in sub grids based on gridnr
    grids = {}
    trajs = {}
    dkappas = {}

    init = time.time()
    grid_names = sorted(glob.glob(os.path.join(path,'grid[0-9][0-9].txt')), key = lambda gn: int(gn.split('.')[-2].split('grid')[-1]))
    numbers = [int(gn.split('.')[-2].split('grid')[-1]) for gn in grid_names]


    for n,gn in enumerate(grid_names):
        number = numbers[n]

        # If we are extending grids and already have some iterations, this will not take the subsequent iterations into account
        if 'MAX_GRID_REFINEMENT_IT' in data and not number < data['MAX_GRID_REFINEMENT_IT']:
            continue

        compressed_fn = os.path.join(path,'trajs/compressed_{}.h5'.format(number))
        compressed = os.path.exists(compressed_fn)
        assert compressed
        print("Loading compressed file for grid {}.".format(number))
        load_compressed_trajs(compressed_fn,number,grids,trajs,dkappas)


    if verbose:
        print("Loading trajectories took {} seconds.".format(time.time()-init))

    return grids, trajs, dkappas


def write_colvars(path,filename,rtrajs,rgrids,rkappas,verbose=False):
    if not os.path.exists(os.path.join(path,'colvars')):
        os.makedirs(os.path.join(path,'colvars'), exist_ok = True)

    with open(filename,'w') as g:
        for n,traj in enumerate(rtrajs):
            # Define point
            cvs = rgrids[n]
            # Define kappa
            kappas = rkappas[n]
            if verbose:
                print("{} | {} - {}".format(n," ".join(["{: 2.2f}".format(cv) for cv in cvs])," ".join(["{: 2.2f}".format(k) for k in kappas])))

            # Create a 1D time series of the collective variable
            t = np.arange(0,len(traj)).reshape(-1,1)

            # Save this time series in a file named colvar
            np.savetxt(os.path.join(path,'colvars/colvar_{}'.format(n)), np.hstack((t, traj)))

            # Write the value of the collective variable and the harmonic spring constant
            g.write(os.path.join(path,'colvars/colvar_{}\t'.format(n)) + "\t".join([str(cv) for cv in cvs]) + '\t' + "\t".join([str(kappa) for kappa in kappas]) + '\n')


def generate_fes(source_path,target_path,data,tol=1e-6):
    """
        Based on umbrella integration
    """

    grids, trajs, kappas = load_grid(source_path,data)
    units = get_units(data)
    temp = data['temp'] if 'temp' in data else 300.*kelvin

    # Ravel the trajectories and grids
    rtrajs = np.zeros((0,*trajs[0][:,data['runup']:].shape[1:]))
    rgrids = np.zeros((0,*grids[0].shape[1:]))
    rkappas = np.zeros((0,*kappas[0].shape[1:]))
    for key in grids.keys():
        rtrajs = np.vstack((rtrajs,trajs[key][:,data['runup']:]))
        rgrids = np.vstack((rgrids,grids[key]))
        rkappas = np.vstack((rkappas,kappas[key]))

    #print("The full trajectories shape taken into account is: ", rtrajs.shape)

    # Convert rgrids and rkappas to atomic units
    rkappas /= units**2 # kappa stays in kjmol
    rgrids *= units

    #WHAM analysis
    fn_wham = os.path.join(target_path,'metadata')
    write_colvars(target_path,fn_wham,rtrajs,rgrids,rkappas)

    temp_none, biasses, trajectories = read_wham_input_2D(fn_wham, path_template_colvar_fns='%s', stride=1, verbose=False)

    rows = [np.linspace(0,25,101)*u for u in units]
    hist = Histogram2D.from_wham_c(rows, trajectories, biasses, temp, error_estimate=None, #'mle_p',
                                   verbosity='medium', convergence=tol, Nscf=20000, cv1_output_unit='angstrom', cv2_output_unit='angstrom')

    fes = FreeEnergySurface2D.from_histogram(hist, temp)
    fes.set_ref(ref='min')
    fes.plot(os.path.join(target_path,'fes.pdf'), ncolors=20)
    fes.savetxt(os.path.join(target_path,'fes.txt'))


def analysis(source_path,target_path):
    print('Analysing {}'.format(target_path))
    # Load yaml file for post-processing
    if os.path.exists(os.path.join(source_path,'data.yml')):
        with open(os.path.join(source_path,'data.yml'),'r') as f:
            data = yaml.full_load(f)
    # Convert runup to units of h5steps
    data['runup'] = data['runup']//data['h5steps']
    generate_fes(source_path,target_path,data)


# MAIN
if __name__ == '__main__':
    #for fname in glob.glob('COF-*/*K/[0-9]*fold/'):
    for fname in glob.glob('NPN-*/*K/[0-9]*fold/'):
        source_path = fname
        target_path = os.path.join(os.getcwd(), '/'.join(fname.split('/')[-4:]))
        if not os.path.exists(os.path.join(target_path,'fes.txt')):
            try:
                analysis(source_path,target_path)
            except (FloatingPointError,AssertionError) as e:
                print('Something went wrong! Skipping ...')
                print('----------------------------------')
                continue
