import matplotlib.pyplot as pt
import h5py
import numpy as np
import yaml
import os
import matplotlib.patches as mpatches
from molmod.units import *
from molmod.constants import *
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.spatial.distance import jensenshannon
from scipy.interpolate import interpn

def compute_histogram_bins(bin_centers):
    """
    Given an array of bin centers, returns an array of bin edges that define
    the histogram bins. Assumes uniform bin widths.
    """
    bin_width = bin_centers[1] - bin_centers[0] # Assumes uniform bin widths
    bin_edges = np.concatenate(([bin_centers[0] - bin_width/2],
                                 bin_centers + bin_width/2))
    return bin_edges

# Load FES
fes = np.loadtxt('fes.dat')
fes[:,1] -= np.nanmin(fes[:,1])

beta = 1/(boltzmann*300*kelvin/kjmol)


# Load traj data
layer00 = np.genfromtxt('layer00.txt', delimiter=',',dtype=object,skip_header=1,encoding='utf')

pt.figure()
num_plots_per_row = 1
nrows = int(np.ceil(float(len(layer00))/num_plots_per_row))
ncols = num_plots_per_row

def get_indices(idx):
    return idx//num_plots_per_row, idx%num_plots_per_row

fig,axgrid = pt.subplots(nrows, ncols,figsize=(5, 20))

for n,node in enumerate(layer00):
    traj_data = node
    biased_fes = fes.copy()
    bias = 0.5*float(traj_data[3])*(biased_fes[:,0]-float(traj_data[2]))**2
    biased_fes[:,1] += bias

    """pt.figure()
    pt.plot(fes[:,0],fes[:,1])
    pt.plot(fes[:,0],bias)
    pt.plot(biased_fes[:,0],biased_fes[:,1])
    pt.show()
    pt.close()"""

    if os.path.exists('data.yml'):
        with open('data.yml','r') as f:
            data = yaml.full_load(f)


    binwidths = data['HISTOGRAM_BIN_WIDTHS']
    edges = data['edges']
    spacings = data['spacings']
    RUN_UP_TIME = data['runup']

    with h5py.File('trajs/traj_{}_{}.h5'.format(int(traj_data[0]), int(traj_data[1])),'r') as f:
        traj = f['trajectory/cv_values'][:]

    bin_edges = compute_histogram_bins(biased_fes[:,0])
    h1, edges = np.histogramdd(traj[RUN_UP_TIME:], bins=(bin_edges,), density=True)

    # Calc expected biased probability
    biased_prob = np.exp(-beta*biased_fes[:,1])
    bincenters = (edges[0][:-1] + edges[0][1:])/2.
    index = np.argmin(np.abs(bincenters[np.argmax(h1)]-biased_fes[:,0]))
    #index = np.argmin(np.abs(biased_fes[:,0]-float(traj_data[2])))

    norm_biased_fes = biased_fes[:,1]/np.nanmax(biased_fes[:,1])
    norm_biased_prob = biased_prob/np.nanmax(biased_prob) 

    # Interpolate expected biased probability to histogram bin centers
    #biased_prob_hist = np.round(interpn(tuple((biased_fes[np.isfinite(biased_prob),0],)), biased_prob[np.isfinite(biased_prob)], bincenters, bounds_error=False),10)
    #biased_prob_hist /= np.nansum((edges[0][1:]-edges[0][:-1])*biased_prob_hist)
    #js_div = jensenshannon(h1[np.isfinite(biased_prob_hist)],biased_prob_hist[np.isfinite(biased_prob_hist)])**2
    density_biased_prob = biased_prob/np.nansum((edges[0][1:]-edges[0][:-1])*biased_prob)

    js_div = jensenshannon(h1,density_biased_prob)

    axgrid[n].set_title('{}_{} - JS div = {}'.format(int(traj_data[0]), int(traj_data[1]), js_div))
    axgrid[n].bar(edges[0][:-1],h1,align='edge',width=binwidths,edgecolor='k',zorder=-1,label='sampled prob')
    axgrid[n].plot(biased_fes[:,0],norm_biased_fes*np.max(h1), label='biased FES', c='red')
    #axgrid[n].plot(biased_fes[:,0],norm_biased_prob*np.max(h1), label='biased prob', c='orange')
    axgrid[n].bar(edges[0][:-1], density_biased_prob, align='edge',width=binwidths, label='biased prob', alpha=0.5)

    axgrid[n].tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,        # ticks along the left edge are off
    right=False,       # ticks along the right edge are off
    labelleft=False)  # labels along the left edge are off
    
    #axgrid[n].set_xlim([float(traj_data[2])-2*spacings[0], float(traj_data[2])+2*spacings[0]])
    axgrid[n].set_ylim([-0.1,np.max(h1)+0.1])

    axgrid[n].spines['top'].set_visible(False)
    axgrid[n].spines['right'].set_visible(False)
    axgrid[n].spines['left'].set_visible(False)
    if n==0:
        axgrid[n].legend(bbox_to_anchor=(0.0,0.5))

fig.tight_layout()
pt.savefig('expected_probs.pdf',bbox_inches='tight')
pt.close()