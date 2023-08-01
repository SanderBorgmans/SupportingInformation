#! /usr/bin/python

import numpy as np
import numpy.ma as ma
from molmod.units import kjmol
import matplotlib.pyplot as pt

from ogre.sim.utils_analytic import AbstractPotential

# 2d simple potential from pyretis
class Potential(AbstractPotential):
    '''Potential class'''
    
    def internal_compute(self,gpos):
        gpos += (np.array([self.deriv_x(self.pos),self.deriv_y(self.pos)]))*kjmol
        return self.eval(self.pos)*kjmol
    

    @staticmethod
    def eval(x):
        return (x[0]*x[0]+x[1]*x[1])*(x[0]*x[0]+x[1]*x[1]) - 10 * np.exp(-30 * (x[0]-0.2)*(x[0]-0.2) - 3*(x[1]-0.4)*(x[1]-0.4)) \
                                                           - 10 * np.exp(-30 * (x[0]+0.2)*(x[0]+0.2) - 3*(x[1]+0.4)*(x[1]+0.4))
    @staticmethod
    def deriv_x(x):
        return 2*(x[0]*x[0]+x[1]*x[1])*2*x[0] - 10 * np.exp(-30 * (x[0]-0.2)*(x[0]-0.2) - 3*(x[1]-0.4)*(x[1]-0.4)) * (-30)*2*(x[0]-0.2) \
                                              - 10 * np.exp(-30 * (x[0]+0.2)*(x[0]+0.2) - 3*(x[1]+0.4)*(x[1]+0.4)) * (-30)*2*(x[0]+0.2)

    @staticmethod
    def deriv_y(x):
        return 2*(x[0]*x[0]+x[1]*x[1])*2*x[1] - 10 * np.exp(-30 * (x[0]-0.2)*(x[0]-0.2) - 3*(x[1]-0.4)*(x[1]-0.4)) * (-3)*2*(x[1]-0.4) \
                                              - 10 * np.exp(-30 * (x[0]+0.2)*(x[0]+0.2) - 3*(x[1]+0.4)*(x[1]+0.4)) * (-3)*2*(x[1]+0.4)

