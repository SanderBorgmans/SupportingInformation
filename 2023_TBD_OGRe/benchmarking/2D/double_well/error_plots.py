# this script should calculate the sum of the squared errors for calc vs analytic to determine optimal hyperparameters

import glob, sys, os
import numpy as np
import matplotlib.pyplot as pt
import matplotlib.colors as colors
from matplotlib.markers import MarkerStyle
from pathlib import Path
import pickle

import matplotlib as mpl

from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import AxesGrid

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern

def load_potential_file(file_loc):
    import importlib.util
    spec = importlib.util.spec_from_file_location("module",file_loc)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def load_analytic_fes(x):
    potential = load_potential_file('scripts/potential.py').Potential
    return potential.eval(x)

def get_num_sims(fn):
    data = np.genfromtxt(fn,dtype=str)
    return len(list(set([i for i in data])))

def error_metric(fes,ref):
    return np.sqrt(np.average((fes-ref)**2))
    #return residuals_rp(fes,ref)



#Path('error_plots/minimal').mkdir(parents=True, exist_ok=True)
#Path('error_plots/limited_range').mkdir(parents=True, exist_ok=True)
Path('error_plots/').mkdir(parents=True, exist_ok=True)
Path('error_data/').mkdir(parents=True, exist_ok=True)


# check if the database was already generated
if not os.path.exists('error_data/grid_data.pkl'):
    grids = {}
    grids_error = {}
    grids_info = {}

    scan_dirs = glob.glob('scans_*_*_*')
    combinations = [tuple(fn.split('_')[1:]) for fn in scan_dirs]

    for index,sd in enumerate(scan_dirs):
        print(sd)
        working_directories = glob.glob('{}/*/fes.dat'.format(sd))
        if len(working_directories)==0:
            print('Skipping {}: no fes files'.format(sd))
            continue

        data = {}
        errors = {}
        num_sims = {}
        non_converged = {}
        extreme_kappa = {}
        nan_length = {}

        CONFINEMENT = []
        OVERLAP = []
        MAX_LAYERS = []
        

        for wd in working_directories:
            parameters = wd.split('/')[1].split('_')
            try:
                c,o,k,m = parameters
            except ValueError:
                continue # this is backup or test directory

            data[tuple([c,o,m])] = np.loadtxt(wd)
            non_converged_fn = os.path.join(Path(wd).parent,'non_converged')
            extreme_kappa_fn = os.path.join(Path(wd).parent,'extreme_kappas')
            traj_log_fn = os.path.join(Path(wd).parent,'traj_log.txt')
            non_converged[tuple([c,o,m])] = os.path.exists(non_converged_fn)
            extreme_kappa[tuple([c,o,m])] = os.path.exists(extreme_kappa_fn)
            num_sims[tuple([c,o,m])] = get_num_sims(traj_log_fn) #sum(1 for line in open(traj_log_fn))
            CONFINEMENT.append(c)
            OVERLAP.append(o)
            MAX_LAYERS.append(m)


        for k,v in data.items():
            idx = ((v[:,0] >= -0.5) & (v[:,0] <= 0.5) & (v[:,1] >= -1.0) & (v[:,1] <= 1.0))
            idx_finite = (idx & np.isfinite(v[:,-1]))

            # set zero offset for both analytic and calculated data
            calc_fes = v[idx_finite,-1] - np.nanmin(v[idx_finite,-1])
            analytic_fes = load_analytic_fes(v[idx_finite,0:-1].T)
            analytic_fes -= np.nanmin(analytic_fes)

            errors[k] = error_metric(calc_fes,analytic_fes) 

            # Look at nans
            fes_array = v[idx,-1]

            nans = (~np.isfinite(fes_array))
            total_nans = np.sum(nans)
            nan_length[k] = total_nans

        CONFINEMENT = sorted(list(set(CONFINEMENT)))
        OVERLAP = sorted(list(set(OVERLAP)))
        MAX_LAYERS = sorted(list(set(MAX_LAYERS)))

        grid = np.zeros((len(CONFINEMENT),len(OVERLAP),len(MAX_LAYERS),4)) 
        grid_error = np.zeros((len(CONFINEMENT),len(OVERLAP),len(MAX_LAYERS),1))

        for n,c in enumerate(CONFINEMENT):
            for m,o in enumerate(OVERLAP):
                for j,mi in enumerate(MAX_LAYERS):
                    if (c,o,mi) in errors:
                        grid_error[n,m,j,0] = errors[(c,o,mi)]
                        grid[n,m,j,0] = non_converged[(c,o,mi)]
                        grid[n,m,j,1] = num_sims[(c,o,mi)]
                        grid[n,m,j,2] = extreme_kappa[(c,o,mi)]
                        grid[n,m,j,3] = nan_length[(c,o,mi)]
                    else:
                        grid_error[n,m,j,:] = np.nan
                        grid[n,m,j,:] = np.nan

        if np.isnan(grid).all():
            print('Skipping {}: all nan'.format(sd))
            continue

        combination = list(combinations[index]) # spacing,kappa,kg(,consistency)
        if len(combination)==3: combination += ['0.00']
        combination = tuple(combination)

        grids[combination] = grid
        grids_error[combination] = grid_error
        grids_info[combination] = [CONFINEMENT, OVERLAP, MAX_LAYERS]

    with open('error_data/grid_data.pkl','wb') as fp:
        pickle.dump(grids,fp)

    with open('error_data/grid_info.pkl','wb') as fp:
        pickle.dump(grids_info,fp)

    with open('error_data/grid_error_data.pkl','wb') as fp:
        pickle.dump(grids_error,fp)

else:
    print('Loading pickled data.')
    with open('error_data/grid_data.pkl','rb') as fp:
            grids = pickle.load(fp)

    with open('error_data/grid_info.pkl','rb') as fp:
            grids_info = pickle.load(fp)

    with open('error_data/grid_error_data.pkl','rb') as fp:
            grids_error = pickle.load(fp)

# Determine max and min error
max_error = np.max([np.nanmax(ge[:,:,:,0]) for ge in grids_error.values()])
min_error = np.min([np.nanmin(ge[:,:,:,0]) for ge in grids_error.values()])

# Determine max and min number of simulations
max_sims = np.max([np.max(g[:,:,:,1]) for g in grids.values()])
min_sims = np.min([np.min(g[:,:,:,1]) for g in grids.values()])


def make_error_figure(key,grid,grid_error,error_index=0):
    num_plots_per_row = 5
    nrows = int(np.ceil(float(grid.shape[-2])/num_plots_per_row))
    ncols = num_plots_per_row

    def get_indices(idx):
        return idx//num_plots_per_row, idx%num_plots_per_row

    fig = pt.figure(figsize=(ncols*3,nrows*3))
    axgrid = AxesGrid(fig, 111,
                nrows_ncols=(nrows, ncols),
                axes_pad=0.35,
                cbar_mode='single',
                cbar_location='right',
                cbar_pad=0.15,
                share_all=True
                )

    # Initialize im
    im = None

    CONFINEMENT = grids_info[key][0]
    OVERLAP = grids_info[key][1]
    MAX_LAYERS = grids_info[key][2]

    for n,mi in enumerate(MAX_LAYERS):
        im = axgrid[n].imshow(grid_error[:,:,n,error_index].T, origin='lower', cmap='RdYlGn_r', vmin=0., vmax=5.) # norm=colors.LogNorm(vmin=min_error, vmax=max_error)
        # Plot hatches where the protocol did not converge
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                if grid[i,j,n,0]:
                    axgrid[n].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch='///', fill=False, snap=False))
                if grid[i,j,n,2]:
                    axgrid[n].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch="xxx", fill=False, snap=False))
                if grid[i,j,n,3]>500: # more than 10% nans
                    axgrid[n].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch="...", fill=False, snap=False))
        
        
        # We want to show all ticks... 
        axgrid[n].set_xticks(np.arange(grid.shape[0]))
        axgrid[n].set_yticks(np.arange(grid.shape[1]))
        axgrid[n].title.set_text('MAX_LAYERS = {}'.format(mi))

        i,j = get_indices(n)

        # Create xy labels and edges
        if j==0:
            axgrid[n].set_yticklabels(OVERLAP)
            axgrid[n].set_ylabel('OVERLAP_THR')
        if i==nrows-1:
            axgrid[n].set_xticklabels(CONFINEMENT,rotation = 90)
            axgrid[n].set_xlabel('CONFINEMENT_THR')

    if im is not None:
        cbar = axgrid.cbar_axes[0].colorbar(im)#, ticks=[min_error, max_error])
        #cbar.ax.set_yticklabels(['lowest error', 'highest error']) 

    pt.tight_layout()
    if error_index==0:
        pt.savefig('error_plots/minimal/hyperparameter_scan_{}.pdf'.format("_".join(key)),bbox_inches='tight')
    elif error_index==1:
        pt.savefig('error_plots/full/hyperparameter_scan_{}.pdf'.format("_".join(key)),bbox_inches='tight')
    elif error_index==2:
        pt.savefig('error_plots/limited_range/hyperparameter_scan_{}.pdf'.format("_".join(key)),bbox_inches='tight')
    pt.close()


def make_simulations_figure(key,grid):
    num_plots_per_row = 5
    nrows = int(np.ceil(float(grid.shape[-2])/num_plots_per_row))
    ncols = num_plots_per_row

    def get_indices(idx):
        return idx//num_plots_per_row, idx%num_plots_per_row

    fig = pt.figure(figsize=(ncols*3,nrows*3))
    axgrid = AxesGrid(fig, 111,
                nrows_ncols=(nrows, ncols),
                axes_pad=0.35,
                cbar_mode='single',
                cbar_location='right',
                cbar_pad=0.15,
                share_all=True
                )

    # Initialize im
    im = None

    CONFINEMENT = grids_info[key][0]
    OVERLAP = grids_info[key][1]
    MAX_LAYERS = grids_info[key][2]

    for n,mi in enumerate(MAX_LAYERS):
        im = axgrid[n].imshow(grid[:,:,n,1].T, origin='lower', cmap='RdYlGn_r', vmin=0, vmax=max_sims) # norm=colors.LogNorm(vmin=min_error, vmax=max_error)
        # Plot hatches where the protocol did not converge
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                if grid[i,j,n,0]:
                    axgrid[n].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch='///', fill=False, snap=False))
                if grid[i,j,n,2]:
                    axgrid[n].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch="xxx", fill=False, snap=False))
                if grid[i,j,n,3]>500: # more than 10% nans
                    axgrid[n].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch="...", fill=False, snap=False))
        
        # We want to show all ticks... 
        axgrid[n].set_xticks(np.arange(grid.shape[0]))
        axgrid[n].set_yticks(np.arange(grid.shape[1]))
        axgrid[n].title.set_text('MAX_LAYERS = {}'.format(mi))


        i,j = get_indices(n)

        # Create xy labels and edges
        if j==0:
            axgrid[n].set_yticklabels(OVERLAP)
            axgrid[n].set_ylabel('OVERLAP_THR')
        if i==nrows-1:
            axgrid[n].set_xticklabels(CONFINEMENT,rotation = 90)
            axgrid[n].set_xlabel('CONFINEMENT_THR')

    if im is not None:
        cbar = axgrid.cbar_axes[0].colorbar(im)#, ticks=[min_error, max_error])
        #cbar.ax.set_yticklabels(['lowest error', 'highest error']) 

    pt.tight_layout()
    pt.savefig('error_plots/num_sims_{}.pdf'.format("_".join(key)),bbox_inches='tight')
    pt.close()


def make_error_simulations_figure(key,grid,grid_error,error_index=0):
    num_plots_per_row = 4
    nrows = int(np.ceil(float(grid.shape[-2])/num_plots_per_row)) * 2 # double the number for the simulations
    ncols = num_plots_per_row

    def get_indices(idx):
        return idx//num_plots_per_row, idx%num_plots_per_row
    
    def indices2index(i,j):
        return i*num_plots_per_row + j

    fig = pt.figure(figsize=(ncols*3,nrows*3))
    axgrid = AxesGrid(fig, 111,
                nrows_ncols=(nrows, ncols),
                axes_pad=0.35,
                cbar_mode='edge',
                cbar_location='right',
                cbar_pad=0.15,
                share_all=True
                )

    # Initialize im
    im = None

    CONFINEMENT = grids_info[key][0]
    OVERLAP = grids_info[key][1]
    MAX_LAYERS = grids_info[key][2]

    for n,mi in enumerate(MAX_LAYERS):
        #print(grid[:,:,n,1])
        im0 = axgrid[indices2index(0,n)].imshow(grid_error[:,:,n,error_index].T, origin='lower', cmap='RdYlGn_r', vmin=0., vmax=5.) # norm=colors.LogNorm(vmin=min_error, vmax=max_error)
        im1 = axgrid[indices2index(1,n)].imshow(grid[:,:,n,1].T, origin='lower', cmap='RdYlGn_r', norm=colors.LogNorm(vmin=min_sims, vmax=max_sims))
        
        

        # Plot hatches where the protocol did not converge
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                if grid[i,j,n,0]:
                    axgrid[indices2index(0,n)].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch='///', fill=False, snap=False))
                    axgrid[indices2index(1,n)].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch='///', fill=False, snap=False))
                if grid[i,j,n,2]:
                    axgrid[indices2index(0,n)].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch="xxx", fill=False, snap=False))
                    axgrid[indices2index(1,n)].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch="xxx", fill=False, snap=False))
                if grid[i,j,n,3]>500: # more than 10% nans
                    axgrid[indices2index(0,n)].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch="...", fill=False, snap=False))
                    axgrid[indices2index(1,n)].add_patch(mpl.patches.Rectangle((i-.5, j-.5), 1, 1, hatch="...", fill=False, snap=False))

        
        # We want to show all ticks... 
        axgrid[indices2index(0,n)].set_xticks(np.arange(grid.shape[0]))
        axgrid[indices2index(0,n)].set_yticks(np.arange(grid.shape[1]))
        axgrid[indices2index(1,n)].set_xticks(np.arange(grid.shape[0]))
        axgrid[indices2index(1,n)].set_yticks(np.arange(grid.shape[1]))

        i,j = get_indices(n)

        # Create xy labels and edges
        axgrid[n].title.set_text('MAX_LAYERS = {}'.format(mi))
        axgrid[n+num_plots_per_row].set_xticklabels(CONFINEMENT,rotation = 90)
        axgrid[n+num_plots_per_row].set_xlabel('CONFINEMENT_THR')
        

    axgrid[0].set_yticklabels(OVERLAP)
    axgrid[0].set_ylabel('OVERLAP_THR')
    axgrid[num_plots_per_row].set_yticklabels(OVERLAP)
    axgrid[num_plots_per_row].set_ylabel('OVERLAP_THR')
        

    if im0 is not None:
        cbar = axgrid.cbar_axes[0].colorbar(im0)#, ticks=[min_error, max_error])
        #cbar.ax.set_yticklabels(['lowest error', 'highest error']) 
        cbar.set_label("RMSE [kJ/mol]")

    if im1 is not None:
        cbar = axgrid.cbar_axes[1].colorbar(im1)#, ticks=[min_error, max_error])
        #cbar.ax.set_yticklabels(['lowest error', 'highest error']) 
        cbar.set_label("number of simulations")

    pt.tight_layout()
    pt.savefig('error_plots/combined_sims_error_{}.pdf'.format("_".join(key)),bbox_inches='tight')
    pt.close()


def plot_error_vs_simulations(grids,grids_error,grids_info,error_index=0):
    # single plot over all parameter combinations, for a fixed initial spacing, initial kappa, and kappa growth factor
    pt.figure()

    sizes = {'2':1.0,'10':3.0} # for KGF
    markers = []
    import matplotlib
    cmap_gr = matplotlib.cm.get_cmap('viridis')
    cmap_or = matplotlib.cm.get_cmap('cividis')

    identity_points = []
    identity2_points = []
    error_simulation_points = np.empty((0,2))

    for key,grid in grids.items():
        grid_error = grids_error[key]

        idx = ((grid[:,:,:,0]==0)&(grid[:,:,:,2]==0)&(grid[:,:,:,3]<=500)) # more than 10% nans
        indices = zip(*np.where(idx))
        for i,j,k in indices:
            identity_points.append((*key,))
            identity2_points.append((grids_info[key][0][i],grids_info[key][1][j],grids_info[key][2][k]))

      
        pt.scatter(grid[idx][...,1].ravel(),grid_error[idx][...,error_index].ravel(),s=0.25,label=key,c='#ff7f0e')

        costs = np.array([grid[idx][...,1].ravel(), grid_error[idx][...,error_index].ravel()]).T
        error_simulation_points = np.vstack((error_simulation_points, costs))

    # find the pareto optimal set, and plot it
    optimal_set = is_pareto_efficient(error_simulation_points,return_mask=False)
    optimal_set = optimal_set[error_simulation_points[optimal_set, 0].argsort()]
    pt.plot(error_simulation_points[optimal_set,0],error_simulation_points[optimal_set,1],'k-',lw=0.5,label='Pareto optimal set')


    with open('pareto_set.txt','w') as f:
        f.write('{: <15} {: <15} {: <8} {}\n'.format('intial_params', 'hyper_params', 'num_sims', 'rmse'))
        for o in optimal_set:
            f.write("{: <15} {: <15} {: <8} {}\n".format(",".join(identity_points[o]),",".join(identity2_points[o]),error_simulation_points[o,0],error_simulation_points[o,1]))

        #pt.scatter(grid[:,:,:,1].ravel(),grid_error[:,:,:,error_index].ravel(),s=0.25,label=key)
    pt.xlabel('# simulations')
    pt.ylabel('RMSE')
    pt.ylim([0,5])
    pt.legend()
    pt.tight_layout()

    pt.savefig('error_plots/error_vs_sims.pdf',bbox_inches='tight')
    pt.close()


def plot_densities(grids,grids_error,grids_info):
    from scipy.stats import gaussian_kde
    from scipy.interpolate import griddata
    # single plot over all parameter combinations, for a fixed initial spacing, initial kappa, and kappa growth factor
    fig = pt.figure(figsize=(6*3,1*3))
    axgrid = AxesGrid(fig, 111,
                nrows_ncols=(1, 6),
                axes_pad=0.45,
                cbar_mode='each',
                cbar_location='right',
                cbar_pad=0.15,
                share_all=True
                )

    data = []
    for key,grid in grids.items():
        grid_error = grids_error[key]

        idx = ((grid[:,:,:,0]==0)&(grid[:,:,:,2]==0)&(grid[:,:,:,3]<=500))
        #idx = grid[:,:,:,4]<=100
        indices = zip(*np.where(idx))
        for i,j,k in indices:
            data_point = np.array([grid[i,j,k,1],grid_error[i,j,k,0],float(key[0]), float(key[1]), float(key[2]), float(grids_info[key][0][i]),float(grids_info[key][1][j]), float(key[3])])
            data.append(data_point)

    data_array = np.array(data)
    titles = ['initial spacing', 'initial kappa', 'KAPPA_GROWTH_FACTOR', 'CONFINEMENT_THR', 'OVERLAP_THR', 'CONSISTENCY_THR'] 
    vmins = [0.25,0.1,2,0,0,0.00]
    vmaxs = [0.5,10.0,10,1,1,1.00]

    for n in range(6): # spacing, kappa, kgf, confinement, overlap, confinement   
        # Define grid
        idx = np.array([0,1,2+n])
        data = data_array[:,idx]
        #print(data[:,-1])

        # Define grid parameters
        grid_size = 35
        x_min, x_max = np.min(data[:,0]), np.max(data[:,0])
        y_min, y_max = np.min(data[:,1]), np.max(data[:,1])

        # Create a regular 2D grid
        x_range = np.linspace(x_min, x_max, grid_size)
        y_range = np.linspace(y_min, y_max, grid_size)
        xx, yy = np.meshgrid(x_range, y_range)

        if n==1:
            data[:,2] = np.log(data[:,2])

        # Interpolate the values on the regular 2D grid
        z_interp = griddata((data[:,0], data[:,1]), data[:,2], (xx, yy), method='linear')

        if n==1:
            z_interp = np.exp(z_interp)

        # Plot the density plot
        if n==1:
            im = axgrid[n].imshow(z_interp, cmap='viridis', extent=[x_min, x_max,y_min, y_max], aspect=(x_max-x_min)/(y_max-y_min), interpolation='gaussian', origin='lower', norm=colors.LogNorm(vmin=vmins[n], vmax=vmaxs[n]))
        else:
            im = axgrid[n].imshow(z_interp, cmap='viridis', extent=[x_min, x_max,y_min, y_max], aspect=(x_max-x_min)/(y_max-y_min), interpolation='gaussian', origin='lower', vmin=vmins[n], vmax=vmaxs[n])

        if n==2:
            cbar = axgrid.cbar_axes[n].colorbar(im, ticks=[vmins[n], vmaxs[n]])
        else:
            cbar = axgrid.cbar_axes[n].colorbar(im)
        #cbar.set_label(titles[n])
        axgrid[n].set_xlabel('# simulations')
        axgrid[n].set_title(titles[n])

    axgrid[0].set_ylabel('RMSE')
    #axgrid[0].set_ylim([0,3])


    pt.savefig('error_plots/density_plots.pdf',bbox_inches='tight')
    pt.close()

# Faster than is_pareto_efficient_simple, but less readable.
def is_pareto_efficient(costs, return_mask = True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
        # Look for any dimension where the given point is strictly better than all others 
        nondominated_point_mask = np.any(costs<costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient
    

def consistency_variations(grids,grids_error):
    """
        Consider the average trend of the error, and the number of simulations, as a function of the consistency threshold
    """

    # Get all combinations
    key_combinations = {}
    for key in grids.keys():
        if not tuple([key[0],key[1],key[2]]) in key_combinations.keys():
            key_combinations[tuple([key[0],key[1],key[2]])] = [key[3]]
        else:
            key_combinations[tuple([key[0],key[1],key[2]])] += [key[3]]

    # Remove those combinations without consistency variation
    reduced_key_combinations = {}
    for key,comb in key_combinations.items():
        if len(comb)>1:
            reduced_key_combinations[key] = comb


    data = {}

    for key,comb in reduced_key_combinations.items():
        comb = sorted(comb)
        print(key,comb)
        avg_rmse_list = []
        avg_num_sims_list = []

        # use first entry as the reference
        grid_ref = grids[tuple([*key,comb[0]])][:,:,-1,:]
        grid_error_ref = grids_error[tuple([*key,comb[0]])][:,:,-1,:]

        ref_idx = ((grid_ref[:,:,0]==0)&(grid_ref[:,:,2]==0)&(grid_ref[:,:,3]<=500))

        for n in range(1,len(comb)):
            grid = grids[tuple([*key,comb[n]])][:,:,-1,:]
            grid_error = grids_error[tuple([*key,comb[n]])][:,:,-1,:]

            idx = ((grid[:,:,0]==0)&(grid[:,:,2]==0)&(grid[:,:,3]<=500) & ref_idx)

            avg_num_sims_difference = np.average((grid[idx][...,1]-grid_ref[idx][...,1])/grid_ref[idx][...,1] * 100 )
            avg_rmse_difference = np.average((grid_error[idx][...,0]-grid_error_ref[idx][...,0])/grid_error_ref[idx][...,0] * 100)

            avg_rmse_list.append(avg_rmse_difference)
            avg_num_sims_list.append(avg_num_sims_difference)

        data[tuple([*key,])] = [avg_rmse_list,avg_num_sims_list]
            
   # plot the data as bars
    for key,d in data.items():
        fig,axs = pt.subplots(2,1)

        axs[0].bar(np.arange(len(d[1])), d[1], width=0.25, color='#e41a1c', edgecolor='k')
        axs[1].bar(np.arange(len(d[0])), d[0], width=0.25, color='#4daf4a', edgecolor='k')

        if '0.25' in key[0]:
            axs[0].set_ylim([0,210])
            axs[1].set_ylim([-7.5,7.5])

        if '0.5' in key[0]:
            axs[0].set_ylim([0,350])
            axs[1].set_ylim([-35,0])
        
        #axs[0].set_ylabel(r'$\Delta$ number of simulations (%)')
        #axs[1].set_ylabel(r'$\Delta$ RMSE (%)')
        

        pt.savefig('error_plots/consistency_variation_{}.pdf'.format("_".join(key)),bbox_inches='tight')
        pt.close()

    

"""
# Plot the errors
for key,grid in grids.items():
    #for i in range(3):
    i=0
    make_error_figure(key,grid,grids_error[key],error_index=i)


# Plot the simulations
for key,grid in grids.items():
    make_simulations_figure(key,grid)
"""

"""# make error(num_simulations) plots
#for i in range(3):
i=0
plot_error_vs_simulations(grids,grids_error,grids_info,error_index=i)
plot_densities(grids,grids_error,grids_info)

"""
#Combined error and num simulations plots, with num simulations in log scale (specific number is not important, order of magnitude is)
for key,grid in grids.items():
    make_error_simulations_figure(key,grid,grids_error[key],error_index=0)


#consistency_variations(grids,grids_error)