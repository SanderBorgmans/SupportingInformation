#! /usr/bin/python
import numpy as np, yaml, itertools, copy, os, warnings
import matplotlib.path as Path
import matplotlib.pyplot as pt
import matplotlib.patches as patches

from scipy.spatial import KDTree
from scipy.spatial import Voronoi,voronoi_plot_2d

from molmod.units import *

# This script allows to add new grid points to existing grid_files (and writes these additional points to run.txt)
# Afterwards a plot is generated with the original grid points and the added points
# This script will also alter the edges in data.yml!

def save2yaml(options,mode=None,edges=None):
    if isinstance(options,dict):
        d = options
    else:
        d = vars(options)
    if mode is not None: d.update({'mode' : mode})
    if edges is not None: d.update({'edges' : edges})
    if os.path.exists('data.yml'):
        with open('data.yml','r') as yamlfile:
            cyaml = yaml.safe_load(yamlfile) # Note the safe_load
            cyaml.update(d)
        if cyaml:
            with open('data.yml','w') as yamlfile:
                yaml.safe_dump(cyaml, yamlfile) # Also note the safe_dump
    else:
        with open('data.yml','w') as yamlfile:
            yaml.safe_dump(d, yamlfile) # Also note the safe_dump


def sort_vertices(vertices):
    # sort according to angle, and add first vertex as last one to ensure closure
    com_vertices = vertices - np.average(vertices,axis=0)
    angles = np.arctan2(com_vertices[:,1],com_vertices[:,0])
    sorted_v = vertices[angles.argsort()]
    return np.concatenate((sorted_v,[sorted_v[0]]),axis=0)


def make_path(min,max):
    if len(min) == 1:
        vertices = np.array([[min[0],-1],[min[0],1],[max[0],-1],[max[0],1]])
    elif len(min) == 2:
        vertices = np.array([[min[0],min[1]],[min[0],max[1]],[max[0],min[1]],[max[0],max[1]]])
    else:
        return None
    return Path.Path(sort_vertices(vertices), closed=True)

def make_grid(new_edges, data):
    min = np.asarray(new_edges['min'])
    max = np.asarray(new_edges['max'])
    spacings = data['spacings']
    kappas = data['kappas']

    assert((min<max).all())
    assert len(min) == len(max)

    try:
        assert len(min) == len(spacings)
        spacing = np.array(spacings)
    except AssertionError:
        spacing = np.ones_like(min) * spacings[0]
        print('ATTENTION: USING IDENTICAL SPACING FOR EVERY DOF, SINCE NUMBER OF SPACINGS WAS FAULTY!')

    try:
        assert len(min) == len(kappas)
    except AssertionError:
        kappas = np.ones_like(min) * kappas[0]
        print('ATTENTION: USING IDENTICAL KAPPA FOR EVERY DOF, SINCE NUMBER OF KAPPAS WAS FAULTY!')

    arrs = [np.arange(min[n],max[n]+spacing[n],spacing[n]) for n,_ in enumerate(min)]
    g = np.meshgrid(*arrs,indexing='ij')
    mesh = np.vstack(map(np.ravel, g)).T
    write_grid(mesh,data,make_path(min,max))


def write_grid(points,data,path):
    grid = np.genfromtxt('grid00.txt',skip_header=1,delimiter=',',dtype=np.str, encoding='utf-8')

    old_grid = np.array([[float(v) for v in str(grid[n,2]).split('*')] for n,_ in enumerate(grid)])
    new_grid = []

    tree = KDTree(np.round(np.array([loc for loc in old_grid]),6))

    for p in points:
        idx = tree.query_ball_point(p,0.5)
        if len(idx)==0:
            new_grid.append(p)

    new_grid = np.asarray(new_grid)


    with open('new.txt','w') as f:
        f.write('grid,nr,cvs,kappas\n')
        for n,point in enumerate(new_grid):
            f.write('{},{},{},{}\n'.format(0,len(old_grid)+n, '*'.join(['{:.8f}'.format(p) for p in point]), '*'.join(['{:.8e}'.format(k) for k in data['kappas']]))) # use * as separator


    pt.figure(figsize=(8,8))
    pt.scatter(old_grid[:,0],old_grid[:,1],color='gray')
    pt.scatter(new_grid[:,0],new_grid[:,1],color='red')

    patch = patches.PathPatch(path, lw=2, alpha=0.1)
    ax = pt.gca()
    ax.add_patch(patch)
    pt.show()


if __name__ == '__main__':
    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    # new grid parameters
    edges    = {'min':[0,0], 'max':[22,22]}
    save2yaml(data,edges=edges)

    make_grid(edges, data)
