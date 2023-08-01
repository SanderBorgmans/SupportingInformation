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
        gpos += (0.4*self.pos*self.pos*self.pos - 4*self.pos + 16 * self.pos * np.exp(-self.pos*self.pos))*kjmol
        return  self.eval(self.pos)*kjmol

    @staticmethod
    def eval(x):
        return (0.1*x*x*x*x - 2*x*x - 8*np.exp(-x*x) + 10)