#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as pt
from mpl_toolkits.axes_grid1 import make_axes_locatable


data = np.loadtxt('fes.dat')
grid = data[:,:-1]
fes = data[:,-1]
fes -= np.nanmin(fes)

pt.plot(grid,fes)
pt.xlabel(r'$\xi$ [a.u.]')
pt.ylabel('F [kJ/mol]')

pt.savefig('us.pdf',bbox_inches='tight')
pt.close()