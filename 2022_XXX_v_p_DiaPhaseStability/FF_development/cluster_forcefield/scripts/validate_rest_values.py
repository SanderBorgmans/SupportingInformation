#! /usr/bin/env python3

from __future__ import print_function

import os

from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt

from molmod.units import angstrom, deg
from molmod.ic import _bond_length_low, _bend_angle_low, _dihed_angle_low, _opdist_low

from yaff.system import System

def parse():
    usage = '%prog [options] fn_ai.chk fn_ff.chk'
    descr = 'Compare the rest values of two systems (system 1 = AI (reference), system 2 = FF, but can also be used for other system comparisons)'
    parser = OptionParser(usage = usage, description = descr)
    # Option for the Output file
    parser.add_option(
            '-o', '--output', default = None,
            help = 'File to write all the comparisons to'
    )

    # Options of the FF output files
    parser.add_option(
            '--ff_bond', default = None,
            help = 'File to write the FF bond lengths to'
    )
    parser.add_option(
            '--ff_bend', default = None,
            help = 'File to write the FF bend angles to'
    )
    parser.add_option(
            '--ff_dihed', default = None,
            help = 'File to write the FF dihedral angles to'
    )
    parser.add_option(
            '--ff_oop', default = None,
            help = 'File to write the FF oop distances to'
    )

    # Options of the AI output files
    parser.add_option(
            '--ai_bond', default = None,
            help = 'File to write the AI bond lengths to'
    )
    parser.add_option(
            '--ai_bend', default = None,
            help = 'File to write the AI bend angles to'
    )
    parser.add_option(
            '--ai_dihed', default = None,
            help = 'File to write the AI dihedral angles to'
    )
    parser.add_option(
            '--ai_oop', default = None,
            help = 'File to write the AI oop distances to'
    )

    # Options of the comparison figures
    parser.add_option(
            '--fig_bond', default = None,
            help = 'Figure name of the comparison of the bond lengths'
    )
    parser.add_option(
            '--fig_bend', default = None,
            help = 'Figure name of the comparison of the bend angles'
    )
    parser.add_option(
            '--fig_dihed', default = None,
            help = 'Figure name of the comparison of the dihedral angles'
    )
    parser.add_option(
            '--fig_oop', default = None,
            help = 'Figure name of the comparison of the oop distances'
    )

    options, args = parser.parse_args()
    if not len(args) == 2 or not args[0].endswith('.chk') or not args[1].endswith('.chk'):
        raise IOError('Exactly two arguments expected: the AI and FF Molmod CHK Files (in that order)')
    fn_ai, fn_ff = args

    if options.ff_bond == None:
        options.ff_bond = fn_ff.replace('.chk', '_bond.txt')
    if options.ff_bend == None:
        options.ff_bend = fn_ff.replace('.chk', '_bend.txt')
    if options.ff_dihed == None:
        options.ff_dihed = fn_ff.replace('.chk', '_dihed.txt')
    if options.ff_oop == None:
        options.ff_oop = fn_ff.replace('.chk', '_oop.txt')

    if options.ai_bond == None:
        options.ai_bond = fn_ai.replace('.chk', '_bond.txt')
    if options.ai_bend == None:
        options.ai_bend = fn_ai.replace('.chk', '_bend.txt')
    if options.ai_dihed == None:
        options.ai_dihed = fn_ai.replace('.chk', '_dihed.txt')
    if options.ai_oop == None:
        options.ai_oop = fn_ai.replace('.chk', '_oop.txt')
    
    path = os.path.dirname(fn_ff)
    if options.fig_bond == None:
        options.fig_bond = os.path.join(path, 'FF_AI_bond.pdf')
    if options.fig_bend == None:
        options.fig_bend = os.path.join(path, 'FF_AI_bend.pdf')
    if options.fig_dihed == None:
        options.fig_dihed = os.path.join(path, 'FF_AI_dihed.pdf')
    if options.fig_oop == None:
        options.fig_oop = os.path.join(path, 'FF_AI_oop.pdf')

    if options.output == None:
        fn_compare = os.path.join(path, 'FF_AI_rest_values.txt')
    else:
        fn_compare = options.output

    output_ff_names = [options.ff_bond, options.ff_bend, options.ff_dihed, options.ff_oop]
    output_ai_names = [options.ai_bond, options.ai_bend, options.ai_dihed, options.ai_oop]
    figs_compare = [options.fig_bond, options.fig_bend, options.fig_dihed, options.fig_oop]
    return fn_ff, fn_ai, output_ff_names, output_ai_names, figs_compare, fn_compare

def get_relative(system, i, j):
    dr_ij = system.pos[i] - system.pos[j]
    if system.cell.nvec > 0:
        system.cell.mic(dr_ij)
    return dr_ij

def get_bond_length(system, i, j):
    # Value
    bond = get_relative(system, i, j)
    value = _bond_length_low(bond, 0)[0]
    # FFatype
    itype = system.get_ffatype(i)
    jtype = system.get_ffatype(j)
    ffatype = '.'.join(sorted([itype, jtype]))
    return value, ffatype

def get_bend_angle(system, i, j, k):
    # Value
    bond_ij = get_relative(system, i, j)
    bond_kj = get_relative(system, k, j)
    value = _bend_angle_low(bond_ij, bond_kj, 0)[0]
    # FFatype
    itype = system.get_ffatype(i)
    jtype = system.get_ffatype(j)
    ktype = system.get_ffatype(k)
    if itype < ktype:
        ffatype = '{}.{}.{}'.format(itype, jtype, ktype)
    else:
        ffatype = '{}.{}.{}'.format(ktype, jtype, itype)
    return value, ffatype

def get_dihedral_angle(system, i, j, k, l):
    # Value
    bond_ij = get_relative(system, i, j)
    bond_kj = get_relative(system, k, j)
    bond_lk = get_relative(system, l, k)
    value = _dihed_angle_low(bond_ij, bond_kj, bond_lk, 0)[0]
    # FFatype
    itype = system.get_ffatype(i)
    jtype = system.get_ffatype(j)
    ktype = system.get_ffatype(k)
    ltype = system.get_ffatype(l)
    if itype < ltype:
        ffatype = '{}.{}.{}.{}'.format(itype, jtype, ktype, ltype)
    else:
        ffatype = '{}.{}.{}.{}'.format(ltype, ktype, jtype, itype)
    return value, ffatype

def get_oop_dist(system, i, j, k, l):
    # l is the central atom
    # Value
    bond_il = get_relative(system, i, l)
    bond_jl = get_relative(system, j, l)
    bond_kl = get_relative(system, k, l)
    value = _opdist_low(bond_il, bond_jl, bond_kl, 0)[0]
    # FFatype
    itype = system.get_ffatype(i)
    jtype = system.get_ffatype(j)
    ktype = system.get_ffatype(k)
    ltype = system.get_ffatype(l)
    ffatype = '{}.{}'.format('.'.join(sorted([itype, jtype, ktype])), ltype)
    return value, ffatype

def get_cell_params(system):
    if system.cell.nvec > 0:
        lengths, angles = system.cell.parameters
        volume = system.cell.volume
        return lengths, angles, volume
    else:
        return None, None, None

def get_rest_values(system, ffatypes = False):
    bonds = []
    bends = []
    diheds = []
    oops = []
    if ffatypes:
        bond_ffatypes = []
        bend_ffatypes = []
        dihed_ffatypes = []
        oop_ffatypes = []
    for i, j in system.iter_bonds():
        bond_length, bond_ffatype = get_bond_length(system, i, j)
        bonds.append(bond_length)
        if ffatypes:
            bond_ffatypes.append(bond_ffatype)
    for i, j, k in system.iter_angles():
        bend_angle, bend_ffatype = get_bend_angle(system, i, j, k)
        bends.append(bend_angle)
        if ffatypes:
            bend_ffatypes.append(bend_ffatype)
    for i, j, k, l in system.iter_dihedrals():
        dihed_angle, dihed_ffatype = get_dihedral_angle(system, i, j, k, l)
        diheds.append(dihed_angle)
        if ffatypes:
            dihed_ffatypes.append(dihed_ffatype)
    for i, j, k, l in system.iter_oops():
        oop_dist, oop_ffatype = get_oop_dist(system, i, j, k, l)
        oops.append(oop_dist)
        if ffatypes:
            oop_ffatypes.append(oop_ffatype)
    lengths, angles, volume = get_cell_params(system)
    if ffatypes:
        return bonds, bond_ffatypes, bends, bend_ffatypes, diheds, dihed_ffatypes, oops, oop_ffatypes, lengths, angles, volume
    else:
        return bonds, bends, diheds, oops, lengths, angles, volume

def correct_dihedrals(diheds_0, diheds_1):
    corrected_diheds_0 = []
    corrected_diheds_1 = []
    for dihed_0, dihed_1 in zip(diheds_0, diheds_1):
        if abs(dihed_0 - dihed_1) > 180*deg:
            if dihed_0 > 0 and dihed_1 < 0:
                dihed_1 += 360*deg
            elif dihed_0 < 0 and dihed_1 > 0:
                dihed_0 += 360*deg
        corrected_diheds_0.append(dihed_0)
        corrected_diheds_1.append(dihed_1)
    return corrected_diheds_0, corrected_diheds_1

def dump_to_files(all_ffatypes, all_ics, fn_names, fn_sys, units = [angstrom, deg, deg, angstrom]):
    # Dump all internal coordinates to a file
    for i in range(4):
        ffatypes = all_ffatypes[i]
        ics = all_ics[i]
        fn_name = fn_names[i]
        unit = units[i]
        ic_types = ['Bond distances', 'Bend angles', 'Dihedral angles', 'Out-of-plane distances']
        unit_names = ['A', 'deg', 'deg', 'A']
        data = []
        for ffatype, ic in zip(ffatypes, ics):
            data.append([ffatype, ic/unit])
        data = np.array(data)
        np.savetxt(fn_names[i], data, fmt = '%s', header = '{} [{}] for the system {}'.format(ic_types[i], unit_names[i], fn_sys))

def compare(ffatypes, ic_ff, ic_ai, figs_compare, fn):
    fn = open(fn, 'w')
    if not ic_ff[4] == None:
        # Print cell characteristics
        lengths_ai = ic_ai[4]
        angles_ai = ic_ai[5]
        volume_ai = ic_ai[6]
        lengths_ff = ic_ff[4]
        angles_ff = ic_ff[5]
        volume_ff = ic_ff[6]
        print('')
        print('{:83s}[{:^4s}]  |       AI       |       FF         '.format('CELL PARAMS', 'UNIT'), file=fn)
        print('------------------------------------------------------------------------------|----------------|----------------', file=fn)
        print('{:83s}[{:^4s}]  |                |                  '.format('  LENGTHS', 'A'), file=fn)
        print('{:89s}  |     {:6.3f}     |     {:6.3f}     '.format('    a', lengths_ai[0], lengths_ff[0]), file=fn)
        print('{:89s}  |     {:6.3f}     |     {:6.3f}     '.format('    b', lengths_ai[1,0], lengths_ff[1]), file=fn)
        print('{:89s}  |     {:6.3f}     |     {:6.3f}     '.format('    c', lengths_ai[2,0], lengths_ff[2]), file=fn)
        print('{:89s}  |                |                '.format(''), file=fn)
        print('{:83s}[{:^4s}]  |                |                '.format('ANGLES', 'deg'), file=fn)
        print('{:89s}  |   {:8.3f}     |   {:8.3f}      '.format('    alpha', angles_ai[0], angles_ff[0]), file=fn)
        print('{:89s}  |   {:8.3f}     |   {:8.3f}      '.format('    beta', angles_ai[1], angles_ff[1]), file=fn)
        print('{:89s}  |   {:8.3f}     |   {:8.3f}      '.format('    gamma', angles_ai[2], angles_ff[2]), file=fn)
        print('{:89s}  |                |                '.format(''), file=fn)
        print('{:83s}[{:^4s}]  |   {:8.3f}     |   {:8.3f}      '.format('VOLUME', 'A^3', volume_ai, volume_ff), file=fn)
        print('', file=fn)
    
    # Print rest values
    print('{:83s}[{:^4s}]  |      RMSD      |       MD       |      RMSE      '.format('REST VALUES', 'UNIT'), file = fn)
    print('-------------------------------------------------------------------------------------------|----------------|----------------|----------------', file=fn)

    rest_value_names = ['BONDS', 'BENDS', 'DIHEDRALS', 'OOPDISTS']
    unit_names = ['A', 'deg', 'deg', 'A']
    units = [angstrom, deg, deg, angstrom]
    for i in range(4):
        rest_value_name = rest_value_names[i]
        unit = units[i]
        unit_name = unit_names[i]
        rest_value_ff = np.array(ic_ff[i])/unit
        rest_value_ai = np.array(ic_ai[i])/unit
        if len(rest_value_ff) > 0:
            rmsd = np.sqrt(((rest_value_ff - rest_value_ai)**2).mean())
            md   = (rest_value_ff - rest_value_ai).mean()
            rmse = np.sqrt(((rest_value_ff - rest_value_ai - md)**2).mean())
            print('{:83s}[{:^4s}]  |    {: 5.3e}   |   {: 5.3e}   |    {: 5.3e}    '.format(rest_value_name, unit_name, rmsd, md, rmse), file = fn)
            ffatype_dict = {}
            for j in range(len(rest_value_ff)):
                ffatype = ffatype_dict.setdefault(ffatypes[i][j], [])
                ffatype.append([rest_value_ai[j], rest_value_ff[j]])
            for ffatype, data in ffatype_dict.items():
                data = np.array(data)
                ffatype_rmsd = np.sqrt(((data[:,1]-data[:,0])**2).mean())
                ffatype_md   = (data[:,1]-data[:,0]).mean()
                ffatype_vmd  = np.sqrt(((data[:,1]-data[:,0]-md)**2).mean())
                print('{:83s}({:^4d})  |    {: 5.3e}   |   {: 5.3e}   |    {: 5.3e}    '.format(ffatype, len(data), ffatype_rmsd, ffatype_md, ffatype_vmd), file=fn)
            print('{:89s}  |                 |                |'.format(''), file=fn)
        else:
            print('{:83s}[{:^4s}]  |         -      |        -       |          -     '.format(rest_value_name, unit_name), file=fn)
    
        # Plot rest values
        if rest_value_name in ['BONDS', 'BENDS', 'DIHEDRALS', 'OOPDISTS'] and len(rest_value_ff) > 0:
            rmsd = np.sqrt(((rest_value_ff - rest_value_ai)**2).mean())
            md   = (rest_value_ff - rest_value_ai).mean()
            rmse = np.sqrt(((rest_value_ff - rest_value_ai - md)**2).mean())
            props = dict(boxstyle='round', facecolor = 'none', alpha=0.5)
            text = '\n'.join(['RMSD = {: 5.3e}'.format(rmsd), 'MD = {: 5.3e}'.format(md), 'RMSE = {: 5.3e}'.format(rmse)])
            fig, ax = plt.subplots()
            ax.text(0.05, 0.95, text, transform = ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
            plt.title('AI vs FF rest {}'.format(rest_value_name.lower()))
            plt.scatter(rest_value_ai, rest_value_ff)
            lims = plt.xlim()
            plt.plot(lims, lims, 'k--', linewidth = 0.7)
            plt.xlim(lims)
            plt.ylim(lims)
            plt.xlabel('Ab Initio rest {} [{}]'.format(rest_value_name.lower(), unit_name))
            plt.ylabel('Force Field rest {} [{}]'.format(rest_value_name.lower(), unit_name))
            plt.savefig(figs_compare[i], bbox_inches = 'tight')
    fn.close()

def main():
    fn_ff, fn_ai, output_ff_names, output_ai_names, figs_compare, fn_compare = parse()
    
    # Create systems
    system_ff = System.from_file(fn_ff)
    system_ai = System.from_file(fn_ai)
    
    # Extract rest values
    bonds_ff, bond_ffatypes, bends_ff, bend_ffatypes, diheds_ff, dihed_ffatypes, oops_ff, oop_ffatypes, lengths_ff, angles_ff, volume_ff = get_rest_values(system_ff, ffatypes = True)
    bonds_ai, bends_ai, diheds_ai, oops_ai, lengths_ai, angles_ai, volume_ai = get_rest_values(system_ai)
    diheds_ff, diheds_ai = correct_dihedrals(diheds_ff, diheds_ai)
    ffatypes = [bond_ffatypes, bend_ffatypes, dihed_ffatypes, oop_ffatypes]
    ic_ff = [bonds_ff, bends_ff, diheds_ff, oops_ff, lengths_ff, angles_ff, volume_ff]
    ic_ai = [bonds_ai, bends_ai, diheds_ai, oops_ai, lengths_ai, angles_ai, volume_ai]

    # Write to files + compare
    dump_to_files(ffatypes, ic_ff, output_ff_names, fn_ff)
    dump_to_files(ffatypes, ic_ai, output_ai_names, fn_ai)
    compare(ffatypes, ic_ff, ic_ai, figs_compare, fn_compare)

if __name__ == '__main__':
    main()
