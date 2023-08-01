import os
import yaml, sys

from molmod.units import *
from ogre import generate_fes

if __name__ == '__main__':
    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    generate_fes(data,error_estimate=None)