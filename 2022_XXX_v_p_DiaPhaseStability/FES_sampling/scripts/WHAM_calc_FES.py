#! /usr/bin/python
import os,sys,copy,warnings,subprocess,glob,yaml
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as pt
import pickle
import networkx as nx

import ndfsampler.post.grid_utils as grid_utils
import ndfsampler.post.sampling as sampling
from ndfsampler.input.utils import get_units

from scipy.spatial import KDTree,Delaunay
from matplotlib.ticker import MaxNLocator

from molmod.constants import boltzmann
from molmod.units import *

def get_executable(name):
    def is_tool(name):
        """Check whether `name` is on PATH and marked as executable."""

        # from whichcraft import which
        from shutil import which

        return which(name) is not None

    if is_tool(name):
        return name
    else:
        if os.path.exists('./'+name) and os.access('./'+name,os.X_OK):
            return './'+name
        else:
            raise IOError('There is no wham program available. Provide either a wham executable in the directory or have it available in the PATH.')

def generate_fes_WHAM(data,edges,suffix,tol=0.00001,periodicity_x=0,periodicity_y=0,attempt=None):
    """
        Based on umbrella integration
    """
    warnings.warn('This method requires the wham and wham-2d scripts to be available in kjmol units. Load the WHAM/...-kj_mol module when on the HPC!')

    grids, trajs, kappas, _, _ = grid_utils.load_grid(data,verbose=False)
    units = get_units(data)
    edges = { k:copy.copy(np.array(v)*units) for k,v in edges.items() }

    spacings = data['spacings'] * units
    steps = [spacing*0.25 for spacing in spacings]
    temp = data['temp'] if 'temp' in data else 300.*kelvin


    # Ravel the trajectories and grids
    rtrajs = np.zeros((0,*trajs[0][:,data['runup']:].shape[1:]))
    rgrids = np.zeros((0,*grids[0].shape[1:]))
    rkappas = np.zeros((0,*kappas[0].shape[1:]))
    for key in grids.keys():
        rtrajs = np.vstack((rtrajs,trajs[key][:,data['runup']:]))
        rgrids = np.vstack((rgrids,grids[key]))
        rkappas = np.vstack((rkappas,kappas[key]))

    #print("The full trajectories shape taken into account is: ", rtrajs.shape)

    # Convert rgrids and rkappas to atomic units
    rkappas /= units**2 # kappa stays in kjmol
    rgrids *= units

    bins = [int(np.round((edges['max'][i]-edges['min'][i])/(steps[i]))) for i,_ in enumerate(edges['min'])]
    #print(bins)

    #WHAM analysis
    filename = 'combined_fes/metadata_{}'.format(suffix)
    outputname = 'combined_fes/free_{}'.format(suffix)
    if os.path.exists(filename):
        os.remove(filename)
    if os.path.exists(outputname):
        os.remove(outputname)

    # remove trajectories which do not lie within the bounds
    mask = np.ones(rtrajs.shape[0],dtype=np.bool)
    for n,traj in enumerate(rtrajs):
        mask_n = np.ones(rtrajs.shape[1],dtype=np.bool)
        for i in range(rtrajs.shape[2]):
            mask_n *= (traj[:,i]<=edges['max'][i]) * (traj[:,i]>=edges['min'][i])
        mask[n] = np.all(mask_n)

    grid_utils.write_colvars(filename,rtrajs[mask],rgrids[mask],rkappas[mask],verbose=False)

    # Find executable
    if rtrajs.shape[-1] == 1:
        exec = get_executable('wham')
        command = exec + ' {} {} {} {} {} 0 {} {} > fes.log'.format(edges['min'][0], edges['max'][0], bins[0], tol, temp, filename, outputname)
    elif rtrajs.shape[-1] == 2:
        exec = get_executable('wham-2d')
        # ./whame-2d min_x max_x bins_x min_y max_y bins_y tolerance temp numpad metadataFile freeEnergyFile use_mask
        command = exec + ' Px={} {} {} {} Py={} {} {} {} {} {} 0 {} {} 1 > combined_fes/fes_{}.log'.format(periodicity_x, edges['min'][0], edges['max'][0], bins[0],
         periodicity_y, edges['min'][1], edges['max'][1], bins[1], tol, temp, filename, outputname, suffix)
        # use_mask will not take into account bins for which there is no data
    else:
        raise NotImplementedError('WHAM has not yet been implemented for more than 2 dimensions, use the UI method.')


    os.system(command)

    # Load FES analysis
    try:
        if rtrajs.shape[-1] == 1:
            free = np.genfromtxt(outputname,usecols=(0,1))
        if rtrajs.shape[-1] == 2:
            free = np.genfromtxt(outputname,usecols=(0,1,2))
    except OSError:
        return [np.nan], [np.nan]

    free = free.reshape(*bins,rtrajs.shape[-1]+1) # one extra dimension for free energy
    fes = free[...,-1].copy()
    fes = ma.array(fes)
    fes[(fes==np.inf)|(fes==9999999.000000)] = ma.masked

    grid = free[...,:-1]/units # convert from atomic units to defined units

    #grid_utils.write_fes(data,grid,fes,suffix)
    #grid_utils.plot_fes(data,grid,fes,suffix)

    return grid,fes


class GridLinkNode():
    def __init__(self,index,grids,fess,spacings):
        self.index = index
        self.links = self.find_overlapping_grids(grids,fess,spacings)
        self.delta_fes = 0.

        #List information
        self.previous = None # this is in direction of reference


    def find_overlapping_grids(self,grids,fess,spacings):
        overlapping_grids = {}
        for n,_ in enumerate(grids):
            if n == self.index:
                continue
            overlap,common_point = self._overlapping(grids[self.index],fess[self.index],grids[n],fess[n],spacings)
            if overlap:
                overlapping_grids[n] = common_point

        if len(overlapping_grids)==0:
            raise ValueError('Could not find overlapping grid for grid {}.'.format(self.index))

        return overlapping_grids

    @staticmethod
    def _overlapping(grid1,fes1,grid2,fes2,spacings):
        common_point = None
        # For each point in grid2, check if it exists in grid1

        tree = KDTree(grid1.reshape((-1,grid1.shape[-1])))
        idx = tree.query_ball_point(grid2.reshape((-1,grid2.shape[-1])),np.min(spacings)/2)

        # For each point that exists in both grids, check whether both FES values exist
        fes_values_1 = [(np.unravel_index(k[0],fes1.shape), fes1[np.unravel_index(k[0],fes1.shape)]) for n,k in enumerate(idx) if len(k)>0]
        fes_values_2 = [(np.unravel_index(n,fes2.shape),    fes2[np.unravel_index(n,fes2.shape)])    for n,k in enumerate(idx) if len(k)>0]
        common_point_info = [(fes_values_1[n][0], fes_values_2[n][0], fes_values_2[n][1]-fes_values_1[n][1]) for n,_ in enumerate(fes_values_1)
                             if not ma.is_masked(fes_values_1[n][1]) and not ma.is_masked(fes_values_2[n][1])]

        overlap = len(common_point_info)>0
        if overlap:
            common_point = common_point_info[-1]

        return overlap, common_point


class GridLinkList():
    def __init__(self,grid_link_nodes,reference=0):
        self.reference = reference
        self.node_dictionary = {gln.index : gln for gln in grid_link_nodes}
        self.delta_fes_dict = {(gln.index,self.node_dictionary[linked_node].index) : linked_node_info[2] \
                                 for gln in self.node_dictionary.values() \
                                 for linked_node, linked_node_info in gln.links.items()}

        # Create directed graph
        self.graph = self.create_graph()

        # Calculate energy difference with respect to reference
        self.calculate_delta_fes()

    def create_graph(self, plot=True):
        # Construct unique node connections
        node_connections = [(gln.index, self.node_dictionary[linked_node].index) for gln in self.node_dictionary.values() for linked_node in gln.links.keys()]
        node_connections = [tuple(sorted(nc)) for nc in node_connections]
        node_connections = list(set(node_connections))

        # Assign delta fes to each node connection
        node_connections = [(*nc,{'delta_fes':-self.delta_fes_dict[nc]} if nc in self.delta_fes_dict \
                            else {'delta_fes':+self.delta_fes_dict[nc[::-1]]}) for nc in node_connections]

        # For each directed edge, add the inverse as well
        node_connections += [(nc_d[1],nc_d[0],{'delta_fes':-nc_d[2]['delta_fes']}) for nc_d in node_connections]

        # Add node connections to graph with edge attributes
        graph = nx.DiGraph()
        graph.add_edges_from(node_connections)
        if plot:
            pt.figure()
            nx.draw(graph, with_labels=True, font_weight='bold')
            pt.show()

        return graph


    def calculate_delta_fes(self):
        for k,v in self.node_dictionary.items():
            if k==self.reference:
                continue

            # Generate path to reference
            path = nx.shortest_path(self.graph, k, self.reference)
            # Accumulate delta fes contributions
            delta_fes = 0
            for n in range(len(path)-1):
                delta_fes += self.graph.get_edge_data(path[n],path[n+1])['delta_fes']
            v.delta_fes = delta_fes


def daisy_chain_FES(data,grids,fess,edges,plot_sub_fes=False):
    # If finding a common point does not work, try to construct a fully connected network of minimally pairwise overlapping FESs
    spacings = [s*0.25 for s in data['spacings']]  # in units

    reference = 0

    # Create GridLinkNodes
    grid_link_nodes = []
    for n,_ in enumerate(grids):
        grid_link_nodes.append(GridLinkNode(n,grids,fess,spacings))

    # Create GridLinkList
    grind_link_list = GridLinkList(grid_link_nodes,reference=reference)

    # Normalize each fes
    for n,fes in enumerate(fess):
        fes -= grind_link_list.node_dictionary[n].delta_fes

    # Combine fes
    #   Construct full grid
    minima = [np.min([np.min(grid[:,:,i]) for grid in grids]) for i in range(grids[0].shape[-1])]
    maxima = [np.max([np.max(grid[:,:,i]) for grid in grids]) for i in range(grids[0].shape[-1])]

    rows = [np.arange(minima[n],maxima[n]+spacings[n],spacings[n]) for n,_ in enumerate(minima)]
    xx,yy = np.meshgrid(*rows,indexing='ij') # only works for 2D!

    full_grid = np.zeros((*xx.shape,grids[0].shape[-1]))
    full_fes = ma.zeros(xx.shape)

    trees = [KDTree(grid.reshape((-1,grids[0].shape[-1]))) for grid in grids]

    # STEP 3: Merge all the grids
    for i in range(full_grid.shape[0]):
        for j in range(full_grid.shape[1]):
            # Find the grid point in the grids and get the corresponding fes value
            grid_point = np.array([xx[i,j],yy[i,j]])

            idx = [tree.query_ball_point(grid_point,np.min(spacings)/2) for tree in trees]
            fes_values = [fess[n][np.unravel_index(k,fess[n].shape)] for n,k in enumerate(idx) if len(k)>0]
            if len(fes_values)==0:
                full_fes[i,j] = ma.masked
            else:
                if np.all([ma.is_masked(f) for f in fes_values]):
                    full_fes[i,j] = ma.masked
                else:
                    full_fes[i,j] = np.average([f for f in fes_values if not ma.is_masked(f)])
            full_grid[i,j] = grid_point

    # STEP 4: Set reference energy of full grid to 0
    full_fes_min = np.nanmin(full_fes)
    full_fes -= full_fes_min

    # STEP 4bis : Plot all FES components with correct reference energy
    if plot_sub_fes:
        for n,fes in enumerate(fess):
            plot_FES(grids[n],fes-full_fes_min,edges,max_fes=np.nanmax(full_fes),name="combined_fes/fes_{}".format(n))

    return full_grid, full_fes


def combine_FES(data,grids,fess,edges,plot_sub_fes=False):

    # STEP 1: Find common point
    #   Construct full grid
    minima = [np.min([np.min(grid[:,:,i]) for grid in grids]) for i in range(grids[0].shape[-1])]
    maxima = [np.max([np.max(grid[:,:,i]) for grid in grids]) for i in range(grids[0].shape[-1])]

    spacings = [s*0.25 for s in data['spacings']]  # in units

    rows = [np.arange(minima[n],maxima[n]+spacings[n],spacings[n]) for n,_ in enumerate(minima)]
    xx,yy = np.meshgrid(*rows,indexing='ij') # only works for 2D!

    full_grid = np.zeros((*xx.shape,grids[0].shape[-1]))
    full_overlap = ma.zeros(xx.shape,dtype=np.bool)
    full_fes = ma.zeros(xx.shape)

    trees = [KDTree(grid.reshape((-1,grids[0].shape[-1]))) for grid in grids]

    #   Find points for which all grids have a FES value
    for i in range(full_grid.shape[0]):
        for j in range(full_grid.shape[1]):
            # Find the grid point in the grids and get the corresponding fes value
            grid_point = np.array([xx[i,j],yy[i,j]])
            idx = [tree.query_ball_point(grid_point,np.min(spacings)/2) for tree in trees]
            fes_values = [fess[n][np.unravel_index(k,fess[n].shape)] for n,k in enumerate(idx) if len(k)>0]
            if len(fes_values)==len(fess): # all grids should have a value at this point
                full_overlap[i,j] = all([not ma.is_masked(f) for f in fes_values])
            full_grid[i,j] = grid_point

    common_points = full_grid[np.where(full_overlap)]
    if len(common_points)==0:
        print('Could not find common reference, daisy chaining ...')
        full_grid, full_fes = daisy_chain_FES(data,grids,fess,edges,plot_sub_fes=plot_sub_fes)
        return full_grid, full_fes
    common_point = common_points[-1] # this is normally the point with maximal CVs
    print('Found a common reference point at ', common_point)

    # STEP 2: Use the common point as reference
    common_point_idx = [tree.query_ball_point(common_point,np.min(spacings)/2) for tree in trees]
    for n,fes in enumerate(fess):
        assert(len(common_point_idx[n])==1)
        k = common_point_idx[n][0]
        fes -= fes[np.unravel_index(k,fes.shape)]

    # STEP 3: Merge all the grids
    for i in range(full_grid.shape[0]):
        for j in range(full_grid.shape[1]):
            # Find the grid point in the grids and get the corresponding fes value
            grid_point = full_grid[i,j]
            idx = [tree.query_ball_point(grid_point,np.min(spacings)/2) for tree in trees]
            fes_values = [fess[n][np.unravel_index(k,fess[n].shape)] for n,k in enumerate(idx) if len(k)>0]
            if len(fes_values)==0:
                full_fes[i,j] = ma.masked
            else:
                if np.all([ma.is_masked(f) for f in fes_values]):
                    full_fes[i,j] = ma.masked
                else:
                    full_fes[i,j] = np.average([f for f in fes_values if not ma.is_masked(f)])
            full_grid[i,j] = grid_point

    # STEP 4: Set reference energy of full grid to 0
    full_fes_min = np.nanmin(full_fes)
    full_fes -= full_fes_min

    # STEP 4bis : Plot all FES components with correct reference energy (and do not hide any information)
    if plot_sub_fes:
        for n,fes in enumerate(fess):
            plot_FES(grids[n],fes-full_fes_min,edges,max_fes=np.nanmax(full_fes),name="combined_fes/fes_{}".format(n))

    return full_grid, full_fes


def plot_FES(grid,fes,edges,max_fes=400,name="combined_fes"):
    # Grid is expressed in units and does not need to be converted
    fig = pt.figure('fes',figsize = (6,6))
    ax = fig.gca()
    cmap = 'viridis_r'

    fes[fes>max_fes] = ma.masked
    levels = np.arange(0,max_fes+1,10)

    im = ax.contourf(grid[:,:,0],grid[:,:,1],fes,levels=levels,cmap=cmap,vmin=0,vmax=max_fes)
    fig.colorbar(im,label="kJ/mol")

    pt.xlim(edges['min'][0],edges['max'][0])
    pt.ylim(edges['min'][1],edges['max'][1])

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))


    fig.savefig('{}.png'.format(name),bbox_inches='tight')
    pt.close('fes')


def validate_combination(data, edges, name, verbose=False):
    if verbose:
        print('Trying edges ',edges['min'],edges['max'])
    grid_test,fes_test = generate_fes_WHAM(data,edges,name)
    refined = not (np.any(np.isnan(fes_test))) and not fes_test.mask.all() # also check whether all elements are masked
    return refined, grid_test, fes_test


def edge_to_vertices(edges):
    # assume format is ((min_x,min_y),(max_x,max_y)) or {'min':(min_x,min_y), 'max':(max_x,max_y)}
    if isinstance(edges,dict):
        min_x,min_y = edges['min']
        max_x,max_y = edges['max']
    elif isinstance(edges,tuple):
        min_x,min_y = edges[0]
        max_x,max_y = edges[1]
    else:
        raise ValueError('Expected tuple or dict, instead got {}.'.format(type(edges)))
    return [(min_x,min_y),(min_x,max_y),(max_x,min_y),(max_x,max_y)]


def vertices_in_convex_hull(test_edges, succesfull_edges):
    # Create delaunay triangulation of all vertices
    vertices = [e for s in succesfull_edges for e in edge_to_vertices(s)]
    hull = Delaunay(vertices)

    # If any of the vertices of the new edges does not lie inside one of the created triangles, the edges are not redundant
    # function returns simplex index, or -1 if no simplex could be found
    # if all indices are larger than or equal to 0, the vertices all lie inside the hull
    return np.all(hull.find_simplex([edge_to_vertices(test_edges)])>=0)

def get_y_ranges(edges,left,right):
    # STEP 2: Splitting all rects on X coordinates
    relevant_edges = []
    for e in edges:
        if e[0][0]<=left and e[1][0]>=right:
            relevant_edges.append(e)


    def _overlapping(edge,range):
        return not (edge[0][1]>range[1] or edge[1][1]<range[0]) # bottom>max_range or top<min_range

    def _merge(edge,range):
        return (min(edge[0][1],range[0]), max(edge[1][1],range[1]))

    #Step 3: Combining rects on the Y boundaries
    yrs = []
    for e in relevant_edges:
        if len(yrs)==0:
            yrs.append((e[0][1],e[1][1])) # (bottom, top)
        else:
            try:
                # If an overlapping combination exists merge the y limits
                index = next(i for i,yr in enumerate(yrs) if _overlapping(e,yr))
            except StopIteration:
                # else append as independent range
                yrs.append((e[0][1],e[1][1])) # (bottom, top)
            else:
                yrs[index] = _merge(e,yrs[index])

    # Calculate heights
    yrs = [yr[1]-yr[0] for yr in yrs]
    return yrs

def area_of_overlapping_rectangles(edges):
    # algorithm can be found here http://codercareer.blogspot.com/2011/12/no-27-area-of-rectangles.html
    # edges: (left,bottom),(right,top)
    area=0

    # STEP 1: Determine and sort the unique X values (left and right) for all rectangles
    unique_x_vals = sorted(set(e[i][0] for e in edges for i in [0,1])) # get all left and right values and create a set

    # Loop over unique x-ranges (we don't need the final x)
    for n,x in enumerate(unique_x_vals[:-1]):
        left = unique_x_vals[n]
        right = unique_x_vals[n+1]
        xr = right - left # this will be the width of each rectangle in this range
        area += sum([xr*yr for yr in get_y_ranges(edges, left, right)])

    return area

def coverage_unaffected(test_edges, succesfull_edges):
    return np.isclose(area_of_overlapping_rectangles(succesfull_edges),
                      area_of_overlapping_rectangles(succesfull_edges + [test_edges]))

def is_redundant(data, edges, successful_edges):
    # Check whether this combination is redundant given the successful_edges
    # we need to check whether the edge captures some area that was not captured before
    # assume edge format in successful_edges is ((min1,min2),(max1,max2))

    # If there are no succesfull_edges yet, all combinations are relevant
    if len(successful_edges)==0:
        return False


    # We will take a reduction of spacings/2 into account to avoid edge effects in the FES evaluation
    reduced_edges = [(tuple(np.array(s[0]) + np.array(data['spacings'])/2 ),
                      tuple(np.array(s[1]) - np.array(data['spacings'])/2 )) for s in successful_edges]

    reduced_test_edges = (tuple(edges['min']+ np.array(data['spacings'])/2), tuple(edges['max'] - np.array(data['spacings'])/2))

    # APPROACH 1: check whether at least one of the vertices of the test rectangle lies outside the convex hull of the existing rectangles
    # THIS FAILS IN A SITUATION WHERE THE EXISTING RECTANGLES LEAVE A GAP IN THE MIDDLE
    # AND THE REMAINING GAP IS SMALLER THAN THE REQUIRED SPACING (i.e. at least one data['spacing'])
    #return vertices_in_convex_hull(reduced_test_edges, reduced_edges)

    # APPROACH 2: explicitely calculate whether the covered area increases
    return coverage_unaffected(reduced_test_edges, reduced_edges)

def base_refine_lower_bound(data,edges,grids,fess,successful_edges,refine_inverse=False,refinement=[],num_spacings=1,square_edges=False):
    # Refinement defines which parameter should be changed (cv_index,step_size,limit)
    edges_test = edges.copy()

    if square_edges:
        # both min bounds should be refined
        assert len(refinement)==2

        # the min bounds should be equal
        assert np.isclose(edges_test['min'][0],edges_test['min'][1])
    else:
        # only a single parameter should be refined for rectangular edges
        assert len(refinement)==1

        # the min bounds should not be equal, if they are add a step
        if np.isclose(edges_test['min'][0],edges_test['min'][1], atol=refinement[0][1]/2):
            print('Skipping edges ',edges['min'],edges['max'])
            for r in refinement:

                # if one of the parameters has reached its limit, go to the next iteration
                if np.isclose(locals()['edges_test']['min'][r[0]],r[2], atol=refinement[0][1]/2):
                    print("-"*15)
                    return

                # if not, increase parameter with step
                locals()['edges_test']['min'][r[0]] += r[1]*num_spacings

                # finally check whether this combination would already be redundant
                if is_redundant(data, edges_test, successful_edges):
                    #print("-"*15) # do not print these lines
                    return # go the the next iteration (decrease max_edge)

    refined = False
    while not refined:
        if is_redundant(data, edges_test, successful_edges):
            print("-"*15)
            return # go the the next iteration (decrease max_edge)

        # we want the lower bound to be at least one spacing smaller than the upper bond
        # tolerance is spacing/2
        if  ( any([np.isclose(edges_test['max'][i] - edges_test['min'][i], data['spacings'][i], atol=data['spacings'][i]/2) for i,_ in enumerate(edges_test['min'])]) or \
              any([edges_test['max'][i] - edges_test['min'][i] <= data['spacings'][i] for i,_ in enumerate(edges_test['min'])])):
            print('Edge min-max difference too small, ', edges_test['min'], edges_test['max'])
            return # go the the next iteration (decrease max_edge and start from lowest min_edge

        refined, grid_test, fes_test = validate_combination(data, edges_test, 'test', verbose=True)
        if refined:
            print('Found a successful combination ', edges_test['min'], edges_test['max'])
            successful_edges.append((tuple(edges_test['min']), tuple(edges_test['max'])))
            grids.append(grid_test)
            fess.append(fes_test)

            # If the min bounds are different, we have to check the inverse combination as well
            if not square_edges and not refine_inverse and not np.isclose(edges_test['min'][0],edges_test['min'][1]):
                print('Validating inverse combination ...')
                edges_inverse = edges_test.copy()
                edges_inverse['min'] = edges_inverse['min'][::-1]
                edges_inverse['max'] = edges_inverse['max'][::-1]

                r = refinement[0]
                refinement_inverse = [((r[0]+1)%2,r[1],r[2])] # index 0 or 1
                refined_inverse = base_refine_lower_bound(data,edges_inverse,grids,fess,successful_edges,refine_inverse=True,refinement=refinement_inverse,num_spacings=num_spacings)
                if refined_inverse:
                    print('Found the corresponding inverse combination ', tuple(edges_inverse['min']), tuple(edges_inverse['max']))
                else:
                    print('Could not find inverse edges for the following edges: ', edges_test, ' Skipping ...')
        else:
            for r in refinement:
                # if one of the parameters has reached its limit, go to the next iteration
                if np.isclose(locals()['edges_test']['min'][r[0]],r[2]):
                    print("-"*15)
                    return

                # if not, increase parameter with step
                locals()['edges_test']['min'][r[0]] += r[1]*num_spacings
    return refined

def get_edges(data,spacings,initial_minimum_min_edge,initial_maximum_min_edge,initial_minimum_max_edge,initial_maximum_max_edge,grids,fess,successful_edges,rectangular_max_edges,num_spacings=1):

    # STEP 1: Get square edges

    print('---- SQUARE EDGES ----')
    max_bound = initial_maximum_max_edge  # current upper bound

    # Refine upper bound
    while max_bound>=initial_minimum_max_edge:
        # For each upper bound we have to initialize the lower bounds
        min_bound_1 = initial_minimum_min_edge
        min_bound_2 = min_bound_1

        edges = {
        'min': [min_bound_1,min_bound_2],
        'max': [max_bound,max_bound]
        }

        # Refine min_bounds
        refinement = [(0,spacings[0],initial_maximum_min_edge), (1,spacings[0],initial_maximum_min_edge)]
        base_refine_lower_bound(data,edges,grids,fess,successful_edges,refinement=refinement,num_spacings=num_spacings,square_edges=True)
        max_bound-=spacings[0]*num_spacings

    # STEP 2: Get rectangular edges
    print('---- RECTANGULAR EDGES ----')
    max_bound_1 = initial_maximum_max_edge # current upper bound
    max_bound_2 = initial_maximum_max_edge # current upper bound

    # Refine upper bound 2
    while max_bound_2>=initial_minimum_max_edge:
        if rectangular_max_edges:
            max_bound_1 = initial_maximum_max_edge

        # Refine upper bound 1
        while max_bound_1>=initial_minimum_max_edge:
            # For each upper bound we have to initialize the lower bounds
            min_bound_1 = initial_minimum_min_edge  # current lower bound for CV1, this will be optimized in base_refine
            min_bound_2 = initial_minimum_min_edge  # current lower bound for CV2

            # Refine lower bound min 2
            while min_bound_2<=initial_maximum_min_edge:
                edges = {
                'min': [min_bound_1,min_bound_2],
                'max': [max_bound_1,max_bound_2]
                }

                # Refine min_bound_1
                refinement = [(0,spacings[0],initial_maximum_min_edge)]
                base_refine_lower_bound(data,edges,grids,fess,successful_edges,refinement=refinement,num_spacings=num_spacings,square_edges=False)
                min_bound_2+=spacings[0]*num_spacings
            max_bound_1-=spacings[0]*num_spacings

            # If the max bound edges can not be rectangular, keep them identical
            if not rectangular_max_edges:
                max_bound_2 = max_bound_1

        # If the max bound edges can be rectangular alter them
        if rectangular_max_edges:
            max_bound_2-=spacings[0]*num_spacings

def refine_edges(data,initial_minimum_min_edge=0,initial_maximum_min_edge=10,initial_minimum_max_edge=15,initial_maximum_max_edge=20,rectangular_max_edges=False):
    # Get spacings from data
    spacings = [s*0.25 for s in data['spacings']]  # in units
    assert spacings[0]==spacings[1] # this makes it easier

    num_spacings = 4 # determines how many spacings are between each attempt

    # Define data containers
    grids = []
    fess = []
    successful_edges = []

    # Execute get_edges function
    get_edges(data,spacings,initial_minimum_min_edge,initial_maximum_min_edge,initial_minimum_max_edge,
                            initial_maximum_max_edge,grids,fess,successful_edges,
                            num_spacings=num_spacings,rectangular_max_edges=rectangular_max_edges)

    if len(successful_edges) == 0:
        raise ValueError('Could not find edges for the given parameters.')
    else:
        print('Found the following edges:')
        for s in successful_edges:
            print(s)

    return grids,fess

if __name__=='__main__':
    # Load yaml file for post-processing
    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)

    with open('wham.yml','r') as f:
        wham_data = yaml.full_load(f)

    initial_minimum_min_edge = float(wham_data['initial_minimum_min_edge'])
    initial_maximum_min_edge = float(wham_data['initial_maximum_min_edge'])
    initial_minimum_max_edge = float(wham_data['initial_minimum_max_edge'])
    initial_maximum_max_edge = float(wham_data['initial_maximum_max_edge'])

    rectangular_max_edges = True

    save_intermediate = True
    plot_sub_fes = True

    edges = {
    'min':[initial_minimum_min_edge,initial_minimum_min_edge],
    'max':[initial_maximum_max_edge,initial_maximum_max_edge]
    }

    # Convert runup to units of h5steps
    data['runup'] = data['runup']//data['h5steps']

    # save the grid and fes for later
    if not os.path.exists('combined_fes/'):
        os.makedirs('combined_fes')

    if os.path.exists('combined_fes/fes.pkl') and os.path.exists('combined_fes/grid.npy'):
        grid = np.load('combined_fes/grid.npy')
        with open('combined_fes/fes.pkl', "rb") as input_file:
            fes = pickle.load(input_file)
    else:
        if save_intermediate:
            if os.path.exists('combined_fes/fes_0.pkl') and os.path.exists('combined_fes/grid_0.npy'):
                grids = []
                fess = []
                for g in sorted(glob.glob('combined_fes/grid_*.npy'), key=lambda x: int(x.split('/')[-1].split('.')[0].split('_')[-1])):
                    grids.append(np.load(g))
                for f in sorted(glob.glob('combined_fes/fes_*.pkl'), key=lambda x: int(x.split('/')[-1].split('.')[0].split('_')[-1])):
                    with open(f, "rb") as input_file:
                        fess.append(pickle.load(input_file))
            else:
                grids, fess = refine_edges(data,initial_minimum_min_edge=initial_minimum_min_edge,initial_maximum_min_edge=initial_maximum_min_edge,
                                                initial_minimum_max_edge=initial_minimum_max_edge,initial_maximum_max_edge=initial_maximum_max_edge,
                                                rectangular_max_edges=rectangular_max_edges)
                for n,g in enumerate(grids):
                    np.save('combined_fes/grid_{}.npy'.format(n),g)
                for n,f in enumerate(fess):
                    f.dump('combined_fes/fes_{}.pkl'.format(n))
        else:
            grids, fess = refine_edges(data,initial_minimum_min_edge=initial_minimum_min_edge,initial_maximum_min_edge=initial_maximum_min_edge,
                                            initial_minimum_max_edge=initial_minimum_max_edge,initial_maximum_max_edge=initial_maximum_max_edge,
                                            rectangular_max_edges=rectangular_max_edges)

        # combine_FES
        grid,fes = combine_FES(data, grids, fess, edges, plot_sub_fes=plot_sub_fes)
        #grid,fes = daisy_chain_FES(data, grids, fess, edges, plot_sub_fes=plot_sub_fes)


        np.save('combined_fes/grid.npy',grid)
        fes.dump('combined_fes/fes.pkl')

    plot_FES(grid,fes,edges,max_fes=500)
