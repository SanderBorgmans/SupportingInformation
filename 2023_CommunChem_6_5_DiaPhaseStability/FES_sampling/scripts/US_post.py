
import numpy as np, h5py, os, yaml, sys
import matplotlib.pyplot as pt

import warnings

from molmod.units import *
from ndfsampler import *

from yaff import log
log.set_level(0)


def post(refresh=True):
    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    # Convert runup to units of h5steps
    data['runup'] = data['runup']//data['h5steps']

    # Make intermediary plots
    data['plot_con']     = False
    data['plot_overlap'] = False
    data['refresh'] = refresh

    if data['refresh']:
        if os.path.exists('trajs.pkl'):
            os.remove('trajs.pkl')
        if os.path.exists('grids.pkl'):
            os.remove('grids.pkl')
        if os.path.exists('kappas.pkl'):
            os.remove('kappas.pkl')

    investigate_overlap(data)


if __name__ == '__main__':
    post()
