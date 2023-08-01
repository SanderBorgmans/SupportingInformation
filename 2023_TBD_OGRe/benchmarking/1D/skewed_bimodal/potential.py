#! /usr/bin/python

import numpy as np
import numpy.ma as ma
from molmod.units import kjmol
import matplotlib.pyplot as pt

from ogre.sim.utils_analytic import AbstractPotential

# 1d potential with two minima separated by a barrier
class Potential(AbstractPotential):
    '''Potential class'''

    def internal_compute(self,gpos):
        gpos += (3 + 1*(0.04*self.pos**3 - 2*self.pos))*kjmol
        return self.eval(self.pos)*kjmol

    @staticmethod
    def eval(x):
        return 3*x + 1*(0.01*x**4 - x**2)