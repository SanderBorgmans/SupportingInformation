
import numpy as np, h5py, os, yaml, sys
import matplotlib.pyplot as pt

import warnings

from molmod.units import *
from ndfsampler import *

from yaff import log
log.set_level(0)


def post(do_overlap=True,do_generate_FES=False):
    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    # Convert runup to units of h5steps
    data['runup'] = data['runup']//data['h5steps']

    # Make intermediary plots
    data['plot_con']     = False
    data['plot_overlap'] = False
    data['refresh'] = True

    if data['refresh']:
        if os.path.exists('trajs.pkl'):
            os.remove('trajs.pkl')
        if os.path.exists('grids.pkl'):
            os.remove('grids.pkl')
        if os.path.exists('kappas.pkl'):
            os.remove('kappas.pkl')

    if do_overlap:
        investigate_overlap(data)
    if do_generate_FES:
        generate_fes_WHAM(data)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        do_overlap = False
        do_generate_FES = False

        if sys.argv[1]=='overlap':
            do_overlap = True
        elif sys.argv[1]=='fes':
            do_generate_FES = True
        elif sys.argv[1]=='post':
            do_overlap = True
            do_generate_FES = True
        else:
            raise ValueError

        post(do_overlap=do_overlap,do_generate_FES=do_generate_FES)

    else:
        post()
