import os
from optparse import OptionParser
import h5py

from molmod.units import angstrom, meter

from yaff.system import System
from yaff.pes.ff import ForceField
from yaff.sampling.dof import CartesianDOF, StrainCellDOF
from yaff.sampling.opt import CGOptimizer
from yaff.sampling.io import XYZWriter, HDF5Writer

def parse():
    usage = '%prog [options] system.chk'
    descr = 'Yaff optimization of a system'
    parser = OptionParser(usage = usage, description = descr)
    parser.add_option(
            '-c', '--cov', default = 'pars_yaff.txt',
            help = 'Covalent Force Field'
    )
    parser.add_option(
            '-e', '--ei', action = 'append', dest = 'fns_ff', default = [],
            help = 'Electrostatic Force Field'
    )
    parser.add_option(
            '-m', '--mm3', action = 'append', dest = 'fns_ff', default = [],
            help = 'MM3 Force Field'
    )
    parser.add_option(
            '-f', '--ff', action = 'append', dest = 'fns_ff',
            help = 'Other force field files'
    )
    parser.add_option(
            '-o', '--output', default = None,
            help = 'Output optimized CHK System file'
    )
    parser.add_option(
            '-x', '--xyz', action = 'store_true', default = False,
            help = 'If present, an xyz trajectory file is printed'
    )
    parser.add_option(
            '--xyz_steps', type = int, default = 1,
            help = 'Number of steps between each snapshot in the XYZ trajectory'
    )
    parser.add_option(
            '--xyz_name', default = None,
            help = 'Name of the XYZ trajectory file'
    )
    parser.add_option(
            '--hdf5', action = 'store_true', default = False,
            help = 'If present, a HDF5 file is printed'
    )
    parser.add_option(
            '--hdf5_steps', type = int, default = 1,
            help = 'Number of steps between each snapshot in the HDF5 trajectory'
    )
    parser.add_option(
            '--hdf5_name', default = None,
            help = 'Name of the HDF5 trajectory file'
    )
    parser.add_option(
            '-N', type = int, default = None,
            help = 'Number of optimization steps'
    )
    options, args = parser.parse_args()
    if not len(args) == 1 or not args[0].endswith('.chk'):
        raise IOError('Exactly one argument expected: the CHK System file')
    # File names
    fn_sys = args[0]
    path = os.path.dirname(fn_sys)
    fns_ff = [options.cov] + options.fns_ff
    if options.output == None:
        fn_out = fn_sys.replace('.chk', '_opt.chk')
    else:
        fn_out = options.output
    # hooks
    xyz_flag = options.xyz
    xyz_steps = options.xyz_steps
    if options.xyz_name == None:
        xyz_name = os.path.join(path, 'traj_opt.xyz')
    elif options.xyz_name.endswith('.xyz'):
        xyz_name = options.xyz_name
    else:
        if '.' in options.xyz_name:
            raise IOError('XYZ trajectory file {} has a wrong extension (needs .xyz)'.format(options.xyz_name))
        else:
            xyz_name = options.xyz_name + '.xyz'
    hdf5_flag = options.hdf5
    hdf5_steps = options.hdf5_steps
    if options.hdf5_name == None:
        hdf5_name = os.path.join(path, 'traj_opt.h5')
    elif options.hdf5_name.endswith('.h5'):
        hdf5_name = options.hdf5_name
    else:
        if '.' in options.hdf5_name:
            raise IOError('HDF5 trajectory file {} has a wrong extension (needs .h5)'.format(options.hdf5_name))
        else:
            hdf5_name = options.hdf5_name + '.h5'
    # Simulation options
    N = options.N
    return fn_sys, fns_ff, fn_out, [xyz_flag, xyz_steps, xyz_name], [hdf5_flag, hdf5_steps, hdf5_name], N

def optimize(N, system, ff, fn_out, xyz, hdf5):
    if system.cell.nvec == 0:
        dof = CartesianDOF(ff)
    else:
        dof = StrainCellDOF(ff)
    hooks = []
    if hdf5[0]:
        f = HDF5Writer(h5py.File(hdf5[2], mode = 'w'), step = hdf5[1])
        hooks.append(f)
    if xyz[0]:
        f = XYZWriter(xyz[2], step = xyz[1])
        hooks.append(f)
    opt = CGOptimizer(dof, hooks = hooks)
    opt.run(N)
    if opt.dof.converged:
        opt.dof.ff.system.to_file(fn_out)
    else:
        print('System is not converged after {} steps'.format(N))

def main():
    fn_sys, fns_ff, fn_out, xyz, hdf5, N = parse()
    system = System.from_file(fn_sys)
    if system.cell.nvec == 0:
        rcut = 1*meter
    else:
        rcut = 15*angstrom
    ff = ForceField.generate(system, fns_ff, rcut = rcut, smooth_ei = True, alpha_scale = 3.2, gcut_scale = 1.5)
    optimize(N, system, ff, fn_out, xyz, hdf5)
    

if __name__ == '__main__':
    main()


