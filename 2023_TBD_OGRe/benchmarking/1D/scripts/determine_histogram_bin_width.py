import numpy as np
import glob,os
import h5py


std = 1.
for dir in glob.glob('scans_*_*_*/'):
    for traj in glob.glob(os.path.join(dir,'database/*.h5')):
        with h5py.File(traj,'r') as f:
            tr = f['trajectory/cv_values'][:].reshape((1,-1,1))
            std = np.min([std,np.std(tr)])

print(std)