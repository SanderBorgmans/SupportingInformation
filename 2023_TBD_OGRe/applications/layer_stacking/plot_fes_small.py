#! /usr/bin/python

import numpy as np
import numpy.ma as ma
import itertools
import matplotlib.pyplot as pt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.path as Path
import matplotlib.pyplot as pt
import matplotlib.patches as patches

from scipy.spatial import Voronoi,voronoi_plot_2d

from molmod.io import load_chk
from molmod.units import *


def sort_vertices(vertices):
    # sort according to angle, and add first vertex as last one to ensure closure
    com_vertices = vertices - np.average(vertices,axis=0)
    angles = np.arctan2(com_vertices[:,1],com_vertices[:,0])
    sorted_v = vertices[angles.argsort()]
    return np.concatenate((sorted_v,[sorted_v[0]]),axis=0)

def wigner_seitz_cell(vecs,plot=True):
    assert vecs.shape[0]==2
    # make wigner seitz cell boundary
    images = np.array([sum(n * vec for n, vec in zip(ns, vecs)) for ns in itertools.product([-1,0,1],repeat=2)]) # 2 dimensions
    images = images[np.where(np.linalg.norm(images,axis=-1)!=np.linalg.norm(images,axis=-1).max())] # get nearest neighbors
    vor = Voronoi(images)
    if plot:
        fig = voronoi_plot_2d(vor)
        pt.show()
        pt.close()
    return vor,Path.Path(sort_vertices(vor.vertices), closed=True)  # make sure that ordening is right

def get_path(rvecs,plot):
    _,path = wigner_seitz_cell(np.array(rvecs)[0:2,0:2],plot=plot)
    return path


cmap = 'viridis_r'

fig = pt.figure()
ax = pt.gca()

chk_file = 'init.chk'
chk = load_chk(chk_file)
rvecs = chk['rvecs']
rvecs = rvecs/angstrom/4.

# Define potential
data = np.loadtxt('fes.dat')

grid = data[:,:-1]
fes = data[:,-1]
fes -= np.nanmin(fes)

path = get_path(rvecs,False)
mask = path.contains_points(grid)

xdim = len(set(grid[:,0]))
ydim = len(set(grid[:,1]))
grid = grid.reshape((xdim,ydim,2))
fes = fes.reshape((xdim,ydim))
mask = mask.reshape((xdim,ydim))
fes_masked = ma.array(fes,mask=~mask)

levs = np.linspace(0,40,21)
im = pt.contourf(grid[:,:,0],grid[:,:,1],fes_masked,levs,cmap=cmap)


pt.xlim(-5,5)
pt.ylim(-5,5)
ax.set_xlabel(r'$\xi_x$ [Å]')
ax.set_ylabel(r'$\xi_y$ [Å]')

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)
cbar = fig.colorbar(im,cax=cax,orientation="vertical")
cbar.set_label('ΔF [kJ/mol]')

ax.set_aspect('equal')

fig.tight_layout()    
pt.savefig('us_small_clipped.pdf',bbox_inches='tight')
pt.close()