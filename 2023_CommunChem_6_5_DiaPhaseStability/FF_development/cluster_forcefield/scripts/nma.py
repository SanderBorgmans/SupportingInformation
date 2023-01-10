from optparse import OptionParser

import numpy as np

from molmod.units import angstrom, centimeter, meter
from molmod.constants import lightspeed
from molmod.io.chk import load_chk
from molmod.unit_cells import UnitCell

from yaff.system import System
from yaff.pes.ff import ForceField
from yaff.sampling.harmonic import estimate_cart_hessian

from tamkin import Molecule, NMA, ConstrainExt
from tamkin.io.trajectory import dump_modes_gaussian

def parse():
    usage = "%prog [options] system.chk"
    descr = "Use TAMkin to do a NMA analysis on the given system"
    parser = OptionParser(usage = usage, description = descr)
    parser.add_option(
            '-c', '--cov', action = 'append', dest = 'fns_ff', default = [],
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
            '--ff', action = 'append', dest = 'fns_ff', default = None,
            help = 'Force Field file with all parameters, overwrites all other files'
    )
    parser.add_option(
            '-o', '--output', default = None,
            help = 'The Molmod CHK file in which all TAMkin information is stored'
    )
    parser.add_option(
            '-f', '--freq', default = None,
            help = 'The to which the NMA frequencies are dumped (Gaussian LOG file or Numpy TXT file)'
    )
    parser.add_option(
            '-p', '--pertnegfreq', action = 'store_true', default = False,
            help = 'Perturb the system along the largest negative frequency, if one exists'
    )
    options, args = parser.parse_args()
    if not len(args) == 1 or not args[0].endswith('.chk'):
        raise IOError('Exactly 1 argument expected: the CHK System file')
    
    # System file names
    fn_sys = args[0]
    if options.output == None:
        fn_out = fn_sys.replace('.chk', '_nma.chk')
    else:
        fn_out = options.output
    if options.freq == None:
        fn_freq = fn_sys.replace('.chk', '_freq.txt')
    else:
        fn_freq = options.freq
    
    # Force field file names
    fns_ff = options.fns_ff

    return fn_sys, fns_ff, fn_freq, fn_out, options.pertnegfreq

def get_hessian(fn_sys, fns_ff):
    system = System.from_file(fn_sys)
    if system.cell.nvec == 0:
        rcut = 1*meter
    else:
        rcut = 15*angstrom
    ff = ForceField.generate(system, fns_ff, rcut = rcut, smooth_ei = True, alpha_scale = 3.2, gcut_scale = 1.5)
    gpos = np.zeros(system.pos.shape, float)
    energy = ff.compute(gpos = gpos)
    hessian = estimate_cart_hessian(ff)
    return hessian, gpos

def get_molecule(fn_sys, fns_ff):
    system_chk = load_chk(fn_sys)
    # Define unit cell
    if 'rvecs' in system_chk.keys() and len(system_chk['rvecs']) > 0:
        cell = UnitCell(system_chk['rvecs'].T)
    else:
        cell = None
    # Define hessian
    if 'hessian' in system_chk.keys():
        N = len(system_chk['numbers'])
        hessian = system_chk['hessian'].reshape(3*N, 3*N)
        gpos = system_chk['gradient']
    else:
        hessian, gpos = get_hessian(fn_sys, fns_ff)
    mol = Molecule(system_chk['numbers'], system_chk['pos'], system_chk['masses'], 0.0, gpos, hessian, unit_cell = cell)
    return mol

def main():
    # Read arguments
    fn_sys, fns_ff, fn_freq, fn_out, pert = parse()
    
    # Initialize NMA analysis
    molecule = get_molecule(fn_sys, fns_ff)
    #nma = NMA(molecule, ConstrainExt(gradient_threshold=0.001))
    if pert:
        nma = NMA(molecule, ConstrainExt(), do_modes = True)
        pert = np.array([nma.modes[i][0] for i in range(len(nma.modes))])
        sys = System.from_file(fn_sys)
        sys.pos = molecule.coordinates + 1*angstrom*pert.reshape(-1, 3)
        sys.to_file('pert.chk')
    else:
        nma = NMA(molecule, ConstrainExt())

    # Print frequencies
    print('\nThe normal mode frequencies [1/cm] are:\n')
    print(nma.freqs/(lightspeed/centimeter))
    if fn_freq.endswith('.txt'):
        np.savetxt(fn_freq, nma.freqs/(lightspeed/centimeter), header = 'NMA Frequencies in cm-1')
    elif fn_freq.endswith('.log'):
        dump_modes_gaussian(fn_freq, nma)
    else:
        print('Extension of file {} not recognized, frequencies are not stored'.format(fn_freq))
    
    # Output Molmod CHK file
    nma.write_to_file(fn_out)

if __name__ == '__main__':
    main()

