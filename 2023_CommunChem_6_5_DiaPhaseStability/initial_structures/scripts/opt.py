import sys
import os

import numpy as np
np.random.seed(1891)

from yaff import System, ForceField
from yaff.sampling.dof import StrainCellDOF, CartesianDOF
from yaff.sampling.opt import CGOptimizer

from molmod.units import angstrom

def parse():
    if not len(sys.argv) == 2 or not sys.argv[1].endswith('.chk'):
        raise IOError('Expected 1 chk file as argument')
    return sys.argv[1]

def get_ff_kwargs():
    # Return force field arguments
    kwargs = {'rcut': 15.0*angstrom,
              'alpha_scale': 2.86,
              'gcut_scale': 1.0,
              'smooth_ei': True,
              'tailcorrections': True}
    return kwargs

def get_lammps_kwargs(natom):
    # Return lammps arguments
    kwargs = {'fn_table': 'table.dat'}
    if natom > 2000:
        kwargs['kspace'] = 'pppm'
    return kwargs

if __name__ == '__main__':
    fn_sys = parse()
    fn_cluster = ['pars_cluster.txt'] # add polysix where necessary
    fn_uff = 'pars_uff.txt'
    fn_out = fn_sys.replace('.chk', '_opt.chk')
    ff_kwargs = get_ff_kwargs()

    # Initialize
    system = System.from_file(fn_sys)
    ff = ForceField.generate(system, fn_uff, **ff_kwargs)
    dof = StrainCellDOF(ff)
    opt = CGOptimizer(dof)
    opt.run()
    ff = ForceField.generate(system, fn_cluster, **ff_kwargs)
    dof = StrainCellDOF(ff)
    opt = CGOptimizer(dof)
    while True:
        opt.run(10000)
        if opt.dof.converged:
            opt.dof.ff.system.to_file(fn_out)
            break
        if not opt.counter % 10000 == 0:
            break
        sys.stdout.flush()
        opt.dof.ff.system.to_file('save.chk')

