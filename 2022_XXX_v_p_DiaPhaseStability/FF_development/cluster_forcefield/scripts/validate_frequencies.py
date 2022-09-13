#! /usr/bin/env python3

from __future__ import print_function

import os
from optparse import OptionParser

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def parse():
    usage = "%prog [options] fn_ai.txt fn_ff.txt"
    descr = "Compare the frequencies given in fn_ai.txt and fn_nma.txt"
    parser = OptionParser(usage = usage, description = descr)
    parser.add_option(
            '-o', '--output', default = None,
            help = 'Name of the file where all frequencies are compared'
    )
    parser.add_option(
            '-f', '--fig', default = None,
            help = 'Name of the figure where the frequencies are plotted to'
    )
    options, args = parser.parse_args()
    if not len(args) == 2 and not args[0].endswith('.txt') and not args[1].endswith('.txt'):
        raise IOError('Exactly two arguments expected: the AI and the NMA frequencies Numpy TXT files')
    
    # Frequency files
    fn_ai, fn_ff = args
    path = os.path.dirname(fn_ai)
    
    # Output files
    if options.output == None:
        fn_out = os.path.join(path, 'AI_FF_freqs.txt')
    else:
        fn_out = options.output
    if options.fig == None:
        fn_fig = os.path.join(path, 'AI_FF_freqs.pdf')
    else:
        fn_fig = options.fig
    return fn_ai, fn_ff, fn_out, fn_fig

def compare(freq_ai, freq_ff, fn_out, fn_fig):
    # File
    f = open(fn_out, 'w')
    freqs_rmsd = np.sqrt(((freq_ff - freq_ai)**2).mean())
    freqs_md   = (freq_ff - freq_ai).mean()
    freqs_rmse  = np.sqrt(((freq_ff - freq_ai -freqs_md)**2).mean())
    freq_ranges = [
        ['   0 -  100',    0,  100],
        [' 100 -  500',  100,  500],
        [' 500 - 1000',  500, 1000],
        ['1000 - 3000', 1000, 3000],
        ['3000 -  inf', 3000, np.inf],
    ]
    print('FREQUENCIES     [1/cm] |       RMSD       |        MD        |       RMSE     ', file=f)
    print('  All                  |    {: 5.3e}    |    {: 5.3e}    |     {:5.3e}    '.format(freqs_rmsd, freqs_md, freqs_rmse), file=f)
    for label, lower, upper in freq_ranges:
        mask = (lower <= freq_ai) & (freq_ai <= upper)
        N = len(freq_ai[mask])
        if N>0:
            fs_rmsd = np.sqrt(((freq_ff[mask] - freq_ai[mask])**2).mean())
            fs_md   = (freq_ff[mask] - freq_ai[mask]).mean()
            fs_rmse  = np.sqrt(((freq_ff[mask] - freq_ai[mask] - fs_md)**2).mean())
            print('  {:s} ({:3d})    |    {: 5.3e}    |    {: 5.3e}    |    {: 5.3e}   '.format(label, N, fs_rmsd, fs_md, fs_rmse), file=f)
        else:
            print('  {:s}          |         -      |        -       |          -     '.format(label), file=f)
    f.close()

    # Figure
    fig, ax = plt.subplots()

    ax.scatter(freq_ai, freq_ff)
    ax.plot([-100,3500],[-100,3500], '--') # Create reference line x = y
    # Set limits
    ax.set_xlim(-100,3500)
    ax.set_ylim(-100,3500)
    # Set labels
    ax.set_xlabel("Ab initio frequencies (cm$^{-1}$)")
    ax.set_ylabel("NMA frequencies (cm$^{-1}$)")
    
    ## Create zoomed inset plot
    zoom = 2.8
    axins = zoomed_inset_axes(ax, zoom, bbox_to_anchor=(410,170))
    axins.scatter(freq_ai, freq_ff)
    axins.plot([-100,3500],[-100,3500], '--') # Create reference line x = y
    # Set labels
    axins.set_xlim(0, 500)
    axins.set_ylim(0, 500)
    
    # Draw where the inset comes from
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
    # Write RMSD, MD and RMSE
    text = '\n'.join(['RMSD = {: 5.3e}'.format(freqs_rmsd), 'MD = {: 5.3e}'.format(freqs_md), 'RMSE = {: 5.3e}'.format(freqs_rmse)])
    props = dict(boxstyle='round', facecolor = 'none', alpha=0.5)
    ax.text(0.05, 0.95, text, transform = ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.savefig(fn_fig, bbox_inches = 'tight')

def main():
    # Read arguments
    fn_ai, fn_ff, fn_out, fn_fig = parse()
    
    # Load frequencies
    freq_ai = np.loadtxt(fn_ai).ravel()
    freq_ff = np.loadtxt(fn_ff).ravel()

    # Compare frequencies
    compare(freq_ai, freq_ff, fn_out, fn_fig)

if __name__ == '__main__':
    main()
