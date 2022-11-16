import glob
import numpy as np
from molmod.units import *
import h5py
import matplotlib.pyplot as pt
from mpl_toolkits.axes_grid1 import AxesGrid

run_up=999

colors = ['#3579b1','#fec44f','#ce1819']


# Parse average CV values (just from h5)
def average_cv(h5_file):
    h5 = h5py.File(h5_file,'r')
    avg_cv = np.average(h5['trajectory/cv_values'][run_up:],axis=0)/angstrom
    return avg_cv

# Parse average cell parameter values (calculate from h5)
def average_cell_params(h5_file):
    h5 = h5py.File(h5_file,'r')
    vecs = np.array(h5['trajectory/cell'][run_up:])/angstrom
    a_x = vecs[:,0,0]
    a_y = vecs[:,0,1]
    a_z = vecs[:,0,2]

    b_x = vecs[:,1,0]
    b_y = vecs[:,1,1]
    b_z = vecs[:,1,2]

    c_x = vecs[:,2,0]
    c_y = vecs[:,2,1]
    c_z = vecs[:,2,2]

    a = np.sqrt(a_x**2 + a_y**2 + a_z**2)
    b = np.sqrt(b_x**2 + b_y**2 + b_z**2)
    c = np.sqrt(c_x**2 + c_y**2 + c_z**2)

    alpha = np.arccos((b_x*c_x+b_y*c_y+b_z*c_z)/(b*c))*180./np.pi
    beta = np.arccos((a_x*c_x+a_y*c_y+a_z*c_z)/(a*c))*180./np.pi
    gamma = np.arccos((a_x*b_x+a_y*b_y+a_z*b_z)/(a*b))*180./np.pi

    params = np.array([a,b,c,alpha,beta,gamma]).T
    avg_params = np.average(params,axis=0)
    assert len(avg_params)==6

    return avg_params

# Read pore volume from zeo++ output in trajs and average
def average_pore_volume(h5_file):
    """
        Extract accessible pore volume from:
        @ COF-300_5/trajs/traj_0_100/999.vol Unitcell_volume: 3859.76   Density: 0.992433   AV_A^3: 154.854 AV_Volume_fraction: 0.04012 AV_cm^3/g: 0.0404259 NAV_A^3: 29.6429 NAV_Volume_fraction: 0.00768 NAV_cm^3/g: 0.00773856
        Number_of_channels: 1 Channel_volume_A^3: 154.854  
        Number_of_pockets: 7 Pocket_volume_A^3: 20.2251  1.85268  0  7.56513  0  0  0  

    """
    avg_pore_volume = None
    for output_file in glob.glob(h5_file[:-3]+'/*.vol'):
        with open(output_file,'r') as f:
            lines = f.readlines()
        
        avg_pore_volume = float(lines[0].split()[7])

        # Also check number of channels and its volume
        #if not int(lines[1].split()[1]) == 1:
        #    print('Not a single channel {}'.format(h5_file) + ' ({})'.format(int(lines[1].split()[1])))
        #if not np.isclose(float(lines[1].split()[3]),avg_pore_volume):
        #    print('Channel volume ({}) deviates from pore volume ({}) for {}'.format(float(lines[1].split()[3]), avg_pore_volume, h5_file))

    return avg_pore_volume


def avg_cv_to_fivefold(cvs):
    rounded = np.array(np.round(cvs/5,0)*5,dtype=int)
    print(cvs, rounded)
    return np.array(np.round(cvs/5,0)*5,dtype=int)


data = {}
for material in glob.glob('*/'):
    data_array = []
    for traj in glob.glob(material+'trajs/*.h5'):
        avg_cv = np.round(average_cv(traj),1) # convert to closest five fold
        avg_cell_params = np.round(average_cell_params(traj),1)
        avg_pore_volume = int(np.round(average_pore_volume(traj),0))
        data_array.append([*avg_cv,*avg_cell_params,avg_pore_volume])
    data[material[:-1]] = data_array


ncols = 2
nrows = int(len(data)//2)+1
fig = pt.figure(figsize=(ncols*7,nrows*7))
axgrid = AxesGrid(fig, 111,
            nrows_ncols=(nrows, ncols),
            axes_pad=0.4,
            share_all=True,
            label_mode='L'
            )


for n,k in enumerate(sorted(data.keys())):
    v = data[k]
    i,j = n//2, n%2
    ax = axgrid[n]
    ax.set_xlim(3.5,21.5)
    ax.set_ylim(3.5,21.5)
    ax.set_xticks(np.arange(5,25,5))
    ax.set_yticks(np.arange(5,25,5))

    name = k.split('_')
    ax.set_title(name[0]+'({})'.format(name[1]))

    ax.set_ylabel(r'CV$_2$ [Å]')
    ax.set_xlabel(r'CV$_1$ [Å]')
        

    for entry in v:
        if entry[0]>21.5 or entry[1]>21.5:
            continue
        #print(name, entry[0], entry[1])
        ax.scatter([entry[0]], [entry[1]], marker='x', c='k')
        #ax.text(entry[0], entry[1], ",".join([str(e) for e in entry[2:5]]) + '\n' + ",".join([str(e) for e in entry[5:8]]) + '\n' + str(entry[8]), ha='center', va='center')
        #ax.annotate(",".join([str(e) for e in entry[2:5]]) + '\n' + ",".join([str(int(np.round(e,0))) for e in entry[5:8]]) + '\n' + str(entry[8]), (np.max([entry[0]-1.5,3.75]),entry[1]-0.15))

        # print each line in a different colour
        ax.annotate(",".join([str(e) for e in entry[2:5]]),                  (np.max([entry[0]-1.5,3.75]),entry[1]-0.15+1.0), color=colors[0])
        ax.annotate(",".join([str(int(np.round(e,0))) for e in entry[5:8]]), (np.max([entry[0]-1.5,3.75]),entry[1]-0.15+0.5), color=colors[1])
        ax.annotate(str(entry[8]),                                           (np.max([entry[0]-1.5,3.75]),entry[1]-0.15), color=colors[2])

pt.savefig('cv_link.pdf', bbox_inches='tight')
pt.close()
        




