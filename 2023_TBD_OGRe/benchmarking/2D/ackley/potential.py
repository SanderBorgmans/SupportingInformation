#! /usr/bin/python

import numpy as np
import numpy.ma as ma
from molmod.units import kjmol
import matplotlib.pyplot as pt

from ogre.sim.utils_analytic import AbstractPotential

# 2d Ackley potential

class Potential(AbstractPotential):
    '''Potential class'''

    def internal_compute(self,gpos):
        if not self.pos[0]**2 + self.pos[1]**2 == 0:
            gpos += (np.array([2.* np.exp(-0.2 * np.sqrt(0.5*(self.pos[0]**2 + self.pos[1]**2))) * self.pos[0]/(np.sqrt(0.5*(self.pos[0]**2 + self.pos[1]**2))) +  np.exp(0.5*(np.cos(2*np.pi*self.pos[0]) + np.cos(2*np.pi*self.pos[1]))) * np.pi * np.sin(2*np.pi*self.pos[0]),
                               2.* np.exp(-0.2 * np.sqrt(0.5*(self.pos[0]**2 + self.pos[1]**2))) * self.pos[1]/(np.sqrt(0.5*(self.pos[0]**2 + self.pos[1]**2))) +  np.exp(0.5*(np.cos(2*np.pi*self.pos[0]) + np.cos(2*np.pi*self.pos[1]))) * np.pi * np.sin(2*np.pi*self.pos[1]) ]))*kjmol
        return self.eval(self.pos)*kjmol

    @staticmethod
    def eval(x):
        return -20.*np.exp(-0.2 * np.sqrt(0.5*(x[0]**2 + x[1]**2))) - np.exp(0.5*(np.cos(2*np.pi*x[0]) + np.cos(2*np.pi*x[1]))) + np.exp(1) + 20
    