#! /usr/bin/python

import os, glob, sys, copy
import posixpath
import numpy as np
import pickle
import numpy.ma as ma
import copy
import warnings
import mepsa
import sklearn
import matplotlib.pyplot as pt
import matplotlib.patches as mpatches
import matplotlib as mpl
import yaml

from matplotlib.ticker import MaxNLocator
from itertools import combinations
from scipy.spatial.distance import pdist
from scipy import interpolate, signal

colors = [
'#e41a1c',
'#377eb8',
'#4daf4a',
'#984ea3',
'#ff7f00',
'#ffff33',
'#a65628',
'#f781bf'
]

overflow_markers = ['o','s','v']

path_colors = {
('SQ','sq') : '#4daf4a',
('SQ','vrect') : '#984ea3',
('sq','vrect') : '#ff7f00',
('hrect','vrect') : '#e41a1c'
}


# General plot settings
edges = {'min':(0,0), 'max':(25,25)}
level_step = 25
max_fes = 500
cmap = 'viridis_r'


# Define markers for phases
node_markers = {
'SQ' : [(-1,-1),(1,-1),(1,1),(-1,1),(-1,-1),(1,-1)],
'sq' :  [(-1,-1),(1,-1),(1,1),(-1,1),(-1,-1),(1,-1)],
'hrect' :  [(-10,-5),(10,-5),(10,5),(-10,5),(-10,-5),(10,-5)],
'vrect' :  [(-5,-10),(5,-10),(5,10),(-5,10),(-5,-10),(5,-10)],
}

node_colors = {
'SQ' : '#3579b1',
'sq' : '#fec44f',
'hrect' : '#ce1819',
'vrect' : '#ce1819',
}


node_size = {
'SQ' : 100,
'sq' :  30,
'hrect' :  100,
'vrect' :  100,
}

def dia_phase_node_path_generator(target,node_points_to_plot, node_array, points_in_minima_groups, maps_list, fes_data=None, verbose=False):
    """
    This is a custom function for the dia phase stability paper to generate the node indices
    between which the MFEP will be calculated with mepsa based on the provided MIN nodes
    and their energy.
    """

    # Create compound node data
    energy_array = maps_list.flatten()
    nodes = np.array([[*node_points_to_plot[n,:2], int(node_points_to_plot[n,2]),int(points_in_minima_groups[n][0]),energy_array[points_in_minima_groups[n][0]]] \
                                                                                                         for n,_ in enumerate(node_points_to_plot)], dtype=object)
    # Try density based clustering algorithm
    from sklearn.cluster import DBSCAN
    # eps max distance within a cluster
    clustering = DBSCAN(eps=2, min_samples=1).fit(nodes[:,:2])

    nodes_data = {}
    for label in set(clustering.labels_):
        nodes_data[label] = np.array([node for n,node in enumerate(nodes) if clustering.labels_[n]==label])

    #Identify the minimum
    ref_node = nodes[np.argmin(nodes[:,-1])]
    ref_node_label = clustering.labels_[np.argmin(nodes[:,-1])]

    if verbose:
        print('Format: (cv1, cv2, node_index, grid_index, energy)')
        #print('FES minimum identified: ', ref_node, ' with label ', ref_node_label)
        for label,nls in nodes_data.items():
            print('Label {}:'.format(label))
            for n in nls:
                print('\t', n, "==> REF"*(np.all(n==ref_node)))

    # Reduce number of nodes by taking only minimal energy node
    reduced_nodes = {}
    for label,nls in nodes_data.items():
        reduced_nodes[label] = nls[np.argmin([n[-1] for n in nls])]

    # Identify all phases and create paths between all combinations
    fig = pt.figure()
    if fes_data is not None:
        pt.contourf(fes_data[:,:,0],fes_data[:,:,1],fes_data[:,:,2],levels=25,cmap=cmap)

    for label,nls in nodes_data.items():
        pt.scatter(nls[:,0],nls[:,1],c=colors[label%len(colors)],marker=overflow_markers[label//len(colors)],alpha=0.3,edgecolors='k',s=mpl.rcParams['lines.markersize'] ** 2 * 0.5)
        pt.scatter([reduced_nodes[label][0]], [reduced_nodes[label][1]],marker=overflow_markers[label//len(colors)], c=colors[label%len(colors)], edgecolors='k',label='Node index {}'.format(label))


    pt.xlabel(r'CV$_1$ [Å]')
    pt.ylabel(r'CV$_2$ [Å]')
    pt.legend(bbox_to_anchor=(1.04,0.5), loc="center left")
    pt.tight_layout()

    if os.path.exists(os.path.join(target,'phase_idx.yml')):
        with open(os.path.join(target,'phase_idx.yml'),'r') as f:
            phase_idx = yaml.safe_load(f)
        print('Loaded the following phase_idx: ', phase_idx)
        #pt.show()
        pt.close()
    else:
        pt.draw()
        pt.pause(0.1)
        large_square_index = input("Provide node index for large square (or press enter if None): ")
        small_square_index = input("Provide node index for small square (or press enter if None): ")
        horizontal_rectangle_index = input("Provide node index for horizontal rectangle (or press enter if None): ")
        vertical_rectangle_index = input("Provide node index for vertical rectangle (or press enter if None): ")
        pt.close()

        phase_idx = {
        'SQ': int(large_square_index) if len(large_square_index)>0 else None,
        'sq': int(small_square_index) if len(small_square_index)>0 else None,
        'hrect': int(horizontal_rectangle_index) if len(horizontal_rectangle_index)>0 else None,
        'vrect': int(vertical_rectangle_index) if len(vertical_rectangle_index)>0 else None,
        }

        with open(os.path.join(target,'phase_idx.yml'), 'w') as outfile:
            yaml.dump(phase_idx, outfile, default_flow_style=False)


    phase_nodes = {
    k : reduced_nodes[v] for k,v in phase_idx.items() if v is not None
    }

    phase_idx_combinations =  combinations([k for k,v in phase_nodes.items()],2)
    node_paths = {pic : (phase_nodes[pic[0]][2],phase_nodes[pic[1]][2]) for pic in phase_idx_combinations}
    return phase_nodes,node_paths


def adapt_energy_map(key, phase_nodes, x, y, energy_map):
    """
    Adapt the energy such that the other phases in the energy profiel become unavailable, and only the the intended phases are visited
    """
    # Identify which grid_idx need to be adapted
    phases_to_adapt = {k : v for k,v in phase_nodes.items() if not k in key}

    if len(phases_to_adapt)==0:
        return energy_map

    # General settings
    distance_matrix = pdist(np.array([phase[:2] for phase in phase_nodes.values()]))
    radius = np.min(distance_matrix)/2

    # MAYBE WE SHOULD TAKE A UNIQUE RADIUS PER PHASE?
    radius=0

    # For each of the to adapt phases, adapt the energy map
    for _,phase in phases_to_adapt.items():
        xi = (np.abs(x-phase[0]) < radius)
        yi = (np.abs(y-phase[1]) < radius)
        grid_idx = np.outer(xi,yi)
        energy_map[0,grid_idx] = np.nan

    return energy_map

def save_nodes(target,phase_nodes):
    data_nodes = {}
    for k,v in phase_nodes.items():
        data_nodes[k] = {
        'cv1': float(v[0]),
        'cv2': float(v[1]),
        'node_index': int(v[2]),
        'grid_index': int(v[3]),
        'energy': float(v[4]),
        }

    with open(os.path.join(target,'path_nodes.yml'), 'w') as outfile:
        yaml.dump(data_nodes, outfile, default_flow_style=False)

def symmetrize_fes(data):
    #full_data = np.array([data[:,:,2], data[:,:,2].T])
    #data[:,:,2] = np.nanmean(full_data,axis=0)
    data[:,:,2] = 0.5*(data[:,:,2] + data[:,:,2].T)
    return data

def remove_numerical_instability(phase_nodes, forward_barriers, backward_barriers, threshold=2.):
    phases = {}
    for k,v in forward_barriers.items():
        if k[0] not in phases:
            phases[k[0]] = [v]
        else:
            phases[k[0]].append(v)

    for k,v in backward_barriers.items():
        if k[1] not in phases:
            phases[k[1]] = [v]
        else:
            phases[k[1]].append(v)

    for k,v in phases.items():
        if np.all(np.array(v)<threshold):
            phase_nodes = {ki:vi for ki,vi in phase_nodes.items() if not k==ki}
            forward_barriers = {ki:vi for ki,vi in forward_barriers.items() if not k in list(ki)}
            backward_barriers = {ki:vi for ki,vi in backward_barriers.items() if not k in list(ki)}

    return phase_nodes, forward_barriers, backward_barriers

def derive_mfep(source,target,node_path_generator,verbose=False):
    """
    This is a general MFEP derivation function calling the appropriate MEPSA functions (assuming the mepsa python script is importable)

    *** Args

    source
        a filename corresponding to the free energy surface data in three columns (x,y,energy)
    target
        a location to put the calculated data
    node_path_generator
         a function that generates the required origin->target node combinations based on input data containing the MIN nodes and the corresponding energy
         receives the following input data:

         - node_points_to_plot: array with x,y,node index (starting from 0), e.g. [(5.4, 2.1, 0), (24.1,5.0,1)]
         - node_array: flattened array with -2 everywhere and node index at indices that correspond to nodes, e.g. [-2,..,-2,1,-2,...,-2,0,-2,..]
         - points_in_minima_groups: array with flattened grid indices where nodes areg e.g. [[5124],[1069]]
         - energy_array: flattened energy array
         - (optional) fes_data: to make plots
         - (optional) verbose: whether to print information
    """
    # Create target location
    os.makedirs(target, exist_ok=True)

    # Load fes data
    data = np.loadtxt(source)
    bins = (len(set(data[:,0])), len(set(data[:,1])))
    fes_data = data[:,:3].reshape((*bins,3))
    fes_data = symmetrize_fes(fes_data)
    flattened_data = fes_data.reshape(-1,3)
    temp_x, temp_y, temp_z = flattened_data[:,0], flattened_data[:,1], flattened_data[:,2]


    # Format data for mepsa
    x, y, maps_list, error, indices = mepsa.three_dimension_sorter(temp_x, temp_y, temp_z)

    # Find nodes - MIN ONLY
    find_flats = False
    is_4_axes_bool = False # this is apparently a deprecated option
    nodes_not_found, node_points_to_plot, node_array, points_in_minima_groups, node_matrix = mepsa.node_finder(maps_list[np.shape(maps_list)[0]-1],
                                                                                                                is_4_axes_bool, x, y, find_flats)
    # node_points_to_plot: array with x,y,node index (starting from 0)
    # node_array: flattened array with -2 everywhere and node index at indices that correspond to nodes
    # points_in_minima_groups: array with flattened grid indices where nodes are
    # node_matrix: connectivity matrix

    if nodes_not_found:
        warnings.warn('Could not find nodes for this FES. Skipping ...')
        return

    # Some kind of clustering algorithm to assign nodes to the different phases
    phase_nodes,node_paths = node_path_generator(target,node_points_to_plot, node_array, points_in_minima_groups, maps_list, fes_data=fes_data, verbose=verbose)

    if len(phase_nodes)==0:
        warnings.warn('No phases could be identified. Skipping ...')
        return

    # Find paths - based on provided nodes
    paths = {}
    for k,node_path in node_paths.items():
        origin_center = [0,0]
        origin_increment = [0,0]
        target_center = [0,0]
        target_increment = [0,0]

        origin_selection_method = node_path[0] # index
        target_selection_method = node_path[1] # index

        OT_nodes_not_found, origin_point_to_plot, target_point_to_plot, origin_point, target_point = mepsa.find_OT(maps_list[np.shape(maps_list)[0]-1],
                  node_array, x, y,
                  origin_center[0], origin_center[1], origin_increment[0], origin_increment[1],
                  target_center[0], target_center[1], target_increment[0], target_increment[1],
                  origin_selection_method, target_selection_method, points_in_minima_groups)

        # origin_point_to_plot: (x,y,energy)
        # target_point_to_plot: (x,y,energy)
        # origin_point: flattened grid index of origin
        # target_point: flattened grid index of target

        if OT_nodes_not_found:
            warnings.warn('Could not find origin/target path for this FES. Skipping ...')
            return

        # Calculate path energy between origin and target
        mode = "GLOBAL"
        if mode == "GLOBAL":
            if_node_by_node = 0
        elif mode == "NODE BY NODE":
            if_node_by_node = 1


        energy_map = np.empty(maps_list.shape)
        energy_map[:] = maps_list[:]
        energy_map = adapt_energy_map(k,phase_nodes,x,y,energy_map)
        path_not_found, path_points, visited_points = mepsa.find_path(energy_map[np.shape(energy_map)[0]-1], origin_point, target_point, node_array,
                  points_in_minima_groups, node_matrix, is_4_axes_bool, x, y, if_node_by_node, indices, error)

        # path_points: array (x,y,energy,grid index) from origin to target
        # visited_points: array (x,y,energy), I guess this provides a list of all the visited points

        if path_not_found:
            warnings.warn('Could not generate path for this FES. Skipping ...')
            return

        # Save path information
        np.savetxt(os.path.join(target,'path_{}.txt'.format("_".join([ki for ki in k]))), path_points, delimiter='\t')
        paths[k] = path_points

    # Create plots
    #   Make FES plot
    fig, axs = pt.subplot_mosaic([['left', 'middle_top', 'right'],
                                   ['left', 'middle_bottom', 'right']],
                              figsize=(15,7.5), constrained_layout=True,
                              gridspec_kw={'width_ratios': [2,1,1]})

    #fig,axs = pt.subplots(ncols=3,figsize=(15,7.5),gridspec_kw={'width_ratios': [2,1,1]}) # scale needs to be square
    levels = np.arange(0,max_fes+1,level_step)
    im = axs['left'].contourf(fes_data[:,:,0],fes_data[:,:,1],fes_data[:,:,2],levels=levels,cmap=cmap,vmin=0,vmax=max_fes)

    # Add node data to plots
    #   Add markers on FES
    for k,phase in phase_nodes.items():
        if phase is not None:
            axs['left'].scatter([phase[0]],[phase[1]],marker=node_markers[k],c=node_colors[k],edgecolors='k',s=node_size[k],zorder=100)

    #   Add path plots on FES and seperately plot MFEPs as reference
    path_handles = []
    counter = 0
    forward_barriers = {}
    backward_barriers = {}
    for k,path_points in paths.items():
        if 'rect' in k[0] or 'rect' in k[1]:
            if 'rect' in k[0] and 'rect' in k[1]:
                pass
            else:
                if 'hrect' in k[0] or 'hrect' in k[1]:
                    continue

                # all hrects become vrects for plotting except when both show up
                k = [ki.replace('hrect','vrect') for ki in k]


        color = path_colors[tuple(k)] if tuple(k) in path_colors.keys() else path_colors[tuple(k[::-1])]

        #create spline function
        f, u = interpolate.splprep([path_points[:,0], path_points[:,1]], s=0)
        #create interpolated lists of points
        xint, yint = interpolate.splev(np.linspace(0, 1, len(path_points)//2), f)
        axs['left'].plot(xint, yint, lw=2, c=color, solid_capstyle='round')
        #axs[0].plot(path_points[:,0], path_points[:,1], lw=2, c=colors[n], solid_capstyle='round')


        # Calculate local maxima (using smoothed path)
        path_energy = signal.savgol_filter(path_points[:,2], 13, 3)
        axs['right'].plot(np.linspace(0,1,path_points.shape[0]), path_energy, c=color,label="_".join([ki for ki in k]))

        # Use expanded path for better min/max handling
        overlap = 5
        expanded_path_energy = list(path_energy)
        expanded_path_energy = np.array(expanded_path_energy[overlap:0:-1] + expanded_path_energy + expanded_path_energy[-2:-2-overlap:-1])
        max_peaks, _ = signal.find_peaks(expanded_path_energy, height=0, distance=5) # maxima at least 5 indices apart
        min_peaks, _ = signal.find_peaks(-expanded_path_energy, distance=5) # maxima at least 5 indices apart
        min_peaks = sorted(list(set([mp if mp>=overlap and mp<=(len(expanded_path_energy)-1-overlap) else abs(mp-overlap)+overlap if (mp-overlap)<0 else 2*(len(path_energy)-1+overlap)-mp for mp in min_peaks])))
        max_peaks = sorted(list(set([mp if mp>=overlap and mp<=(len(expanded_path_energy)-1-overlap) else abs(mp-overlap)+overlap if (mp-overlap)<0 else 2*(len(path_energy)-1+overlap)-mp for mp in max_peaks])))
        max_peaks = [mp for mp in max_peaks if mp>min_peaks[0] and mp<min_peaks[-1]]

        # min_peaks should at least contain the considered path minima
        if not all([i in min_peaks for i in [overlap, len(expanded_path_energy)-1-overlap]]):
            # remove all indices within the overlap range of the newly added minima
            min_peaks = [mp for mp in min_peaks if np.abs(overlap-mp)>overlap and np.abs(len(expanded_path_energy)-1-overlap-mp)>overlap]

            min_peaks.append(overlap)
            min_peaks.append(len(expanded_path_energy)-1-overlap)
            min_peaks = sorted(list(set((min_peaks))))

        #print(k, max_peaks, min_peaks)

        if len(max_peaks)>0:
            left_indices = [i for i in max_peaks if i < min_peaks[1]]
            left_index = left_indices[0] if len(left_indices)>0 else overlap # first max index that is smaller than the second minimum

            right_indices = [i for i in max_peaks if i > min_peaks[-2]]
            right_index = right_indices[-1] if len(right_indices)>0 else len(path_energy)-1+overlap # last max index that is larger than the second last minimum
        else:
            argmax = np.where(np.isclose(path_energy,np.max(path_energy)))[0]
            left_index, right_index = argmax[0]+overlap, argmax[-1]+overlap

        axs['right'].scatter([(left_index-overlap)/len(path_energy)], [expanded_path_energy[left_index]], marker='<' ,c=color)
        axs['right'].scatter([(right_index-overlap)/len(path_energy)], [expanded_path_energy[right_index]], marker='>' ,c=color)

        forward_barrier = expanded_path_energy[left_index]-path_points[0,2]
        backward_barrier = expanded_path_energy[right_index]-path_points[-1,2]

        forward_barriers[tuple(k)] = np.clip(forward_barrier,0,None)
        backward_barriers[tuple(k)] = np.clip(backward_barrier,0,None)

        axs['middle_top'].bar([counter], [forward_barrier], color=color, width=0.5, align='center', edgecolor='k')
        axs['middle_bottom'].bar([counter], [backward_barrier], color=color, width=0.5, align='center', edgecolor='k')

        counter+=1

    for n,(k,_) in enumerate(forward_barriers.items()):
        axs['middle_top'].scatter([n-0.15],[forward_barriers[k] + np.nanmax(list(forward_barriers.values()))/20],
                                   marker=node_markers[k[0]],c=node_colors[k[0]],edgecolors='k',s=node_size[k[0]],zorder=100)
        axs['middle_top'].scatter([n+0.15],[forward_barriers[k] + np.nanmax(list(forward_barriers.values()))/20],
                                   marker=node_markers[k[1]],c=node_colors[k[1]],edgecolors='k',s=node_size[k[1]],zorder=100)

        axs['middle_top'].scatter([n],[forward_barriers[k] + np.nanmax(list(forward_barriers.values()))/20],
                                   marker='$->$',c='k',s=25,zorder=100)


        axs['middle_bottom'].scatter([n+0.15],[backward_barriers[k] + np.nanmax(list(backward_barriers.values()))/20],
                                   marker=node_markers[k[0]],c=node_colors[k[0]],edgecolors='k',s=node_size[k[0]],zorder=100)
        axs['middle_bottom'].scatter([n-0.15],[backward_barriers[k] + np.nanmax(list(backward_barriers.values()))/20],
                                   marker=node_markers[k[1]],c=node_colors[k[1]],edgecolors='k',s=node_size[k[1]],zorder=100)

        axs['middle_bottom'].scatter([n],[backward_barriers[k] + np.nanmax(list(backward_barriers.values()))/20],
                                   marker='$->$',c='k',s=25,zorder=100)


    # Remove numerical instability stable phases
    phase_nodes, forward_barriers, backward_barriers = remove_numerical_instability(phase_nodes, forward_barriers, backward_barriers)

    # Save the phase nodes for later (e.g. plotting relative energies)
    save_nodes(target,phase_nodes)

    #   Save transition barriers
    transition_barriers = {}
    transition_barriers.update({k : (float(forward_barriers[k]), float(backward_barriers[k])) for k in forward_barriers.keys()})
    with open(os.path.join(target,'path_transitions.yml'), 'w') as outfile:
        yaml.dump(transition_barriers, outfile, default_flow_style=False)

    #   Format plots
    axs['left'].xaxis.set_major_locator(MaxNLocator(integer=True))
    axs['left'].yaxis.set_major_locator(MaxNLocator(integer=True))
    axs['left'].set_xlim(edges['min'][0],edges['max'][0])
    axs['left'].set_ylim(edges['min'][1],edges['max'][1])
    axs['left'].set_xticks(np.arange(edges['min'][0],edges['max'][0]+1,5,dtype=int))
    axs['left'].set_yticks(np.arange(edges['min'][1],edges['max'][1]+1,5,dtype=int))
    axs['left'].set_xlabel(r'CV$_1$ [Å]')
    axs['left'].set_ylabel(r'CV$_2$ [Å]')

    axs['middle_top'].set_ylabel(r'$\Delta F^{‡}_{\rightarrow}$ [kJ/mol]')
    axs['middle_bottom'].set_ylabel(r'$\Delta F^{‡}_{\leftarrow}$ [kJ/mol]')

    axs['middle_top'].set_xticks([])
    axs['middle_bottom'].set_xticks([])

    axs['middle_bottom'].invert_yaxis()
    if counter>0:
        axs['middle_bottom'].set_xlim(-0.5,counter-0.5)
        axs['middle_top'].set_xlim(-0.5,counter-0.5)

    axs['right'].set_xlabel('MFEP path [-]')
    axs['right'].set_ylabel(r'$\Delta$F [kJ/mol]')
    axs['right'].set_xticks([])

    pt.savefig(os.path.join(target,'MFEP.pdf'),bbox_inches='tight')
    pt.close()


if __name__=='__main__':
    COFs = ['COF-300','COF-320','NPN-1','NPN-3']
    path_to_fes = '2DFES_thermolib/'
    for c in COFs:
        c_source = posixpath.join(path_to_fes,c)
        for source in glob.glob(posixpath.join(c_source,'*K','[0-9]*fold','fes.txt')):
            target = posixpath.join(*source.split('/')[-4:-1])
            print("Analyzing ", target)
            derive_mfep(source,target,dia_phase_node_path_generator)
