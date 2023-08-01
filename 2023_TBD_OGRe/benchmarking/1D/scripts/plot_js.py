import matplotlib.pyplot as pt
import h5py, glob, os, yaml
import numpy as np
from scipy.spatial.distance import jensenshannon


from matplotlib.cm import get_cmap

name = "tab20"
cmap = get_cmap(name)  # type: matplotlib.colors.ListedColormap
colors = cmap.colors  # type: list


# Load yaml file for post-processing
if os.path.exists('data.yml'):
    with open('data.yml','r') as f:
        data = yaml.full_load(f)
else:
    raise IOError('No data.yml file found!')


binwidths = data['HISTOGRAM_BIN_WIDTHS']
edges = data['edges']
spacings = data['spacings']

bins = tuple([int(((edges['max'][i]+spacings[i])-(edges['min'][i]-spacings[i]))//binwidths[i]) for i,_ in enumerate(spacings)])
run_up = data['runup']

def calc_JS(simulation):
    traj = simulation['trajectory/cv_values'][run_up:]
    h1, _ = np.histogramdd(traj[:len(traj)//2], bins=bins, range=[(edges['min'][i]-spacings[i],
                                                            edges['max'][i]+spacings[i]) for i,_ in enumerate(spacings)], density=True)
    h2, _ = np.histogramdd(traj[len(traj)//2:], bins=bins, range=[(edges['min'][i]-spacings[i],
                                                            edges['max'][i]+spacings[i]) for i,_ in enumerate(spacings)], density=True)
    return jensenshannon(h1,h2)**2


js_dv = {}

for fn in glob.glob('trajs/traj_*_*.h5'):
    identifier = tuple(fn.split('/')[-1].split('.')[0].split('_')[1:])
    with h5py.File(fn, 'r') as fn:
        js_dv[identifier] = calc_JS(fn)

number_of_layers = len(set([int(k[0]) for k in js_dv.keys()]))

fig,axes = pt.subplots(number_of_layers, figsize=(10,20))
if number_of_layers==1:
    axes = [axes]

for i in range(number_of_layers):
    axes[i].set_prop_cycle(color=colors)
    vals = {k:v for k,v in js_dv.items() if k[0]==str(i)}
    vals = sorted(vals.items(), key=lambda item: item[1])
    axes[i].bar(["_".join(v[0]) for v in vals],[v[1] for v in vals],width=0.75)
    for tick in axes[i].get_xticklabels():
        tick.set_rotation(45)
    axes[i].set_ylabel('JS divergence of 0:50 to 50:100')


axes[-1].set_xlabel('Trajectory identifier')
pt.savefig('js_plot.pdf',bbox_inches='tight')
pt.close()


