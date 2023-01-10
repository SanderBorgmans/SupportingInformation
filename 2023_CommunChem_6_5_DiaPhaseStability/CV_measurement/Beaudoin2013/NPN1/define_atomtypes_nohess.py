#! /usr/bin/env python

from molmod.io import dump_chk, load_chk
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph
from molmod.molecular_graphs import HasAtomNumber, HasNumNeighbors, HasNeighborNumbers, HasNeighbors
from molmod.graphs import CritAnd, CritNot, CritOr
from molmod.periodic import periodic as pt
from molmod.units import angstrom

from yaff import System as YaffSystem
from yaff.log import log

import sys, os, numpy as np

from collections import Counter

fn_chk = sys.argv[1]
fn_fin = sys.argv[1]


C_C4_TNM         =   CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,6,6))
C_C2N_0_TNM      =   CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,7))
C_C3_TNM         =   CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,6))
C_HC2_HC2N_0_TNM =   CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), C_C2N_0_TNM))
C_HC2_HC3_TNM    =   CritAnd(HasAtomNumber(6), HasNeighbors(HasAtomNumber(1), HasAtomNumber(6), C_C3_TNM))
H_0_TNM          =   CritAnd(HasAtomNumber(1), HasNeighbors(C_HC2_HC2N_0_TNM))
H_1_TNM          =   CritAnd(HasAtomNumber(1), HasNeighbors(C_HC2_HC3_TNM))
N_0_TNM          =   HasAtomNumber(7)
O_0_TNM          =   HasAtomNumber(8)


afilters = [('C_C4_TNM', C_C4_TNM), ('C_C2N_0_TNM', C_C2N_0_TNM), ('C_C3_TNM', C_C3_TNM), ('C_HC2_HC2N_0_TNM', C_HC2_HC2N_0_TNM), ('C_HC2_HC3_TNM', C_HC2_HC3_TNM),
            ('H_0_TNM', H_0_TNM), ('H_1_TNM', H_1_TNM), ('N_0_TNM', N_0_TNM), ('O_0_TNM', O_0_TNM),
            ]



def get_atypes(yaffsystem):
    graph = MolecularGraph(yaffsystem.bonds, yaffsystem.numbers)
    ffatypes = [0]*len(yaffsystem.numbers)
    test = [False]*len(yaffsystem.numbers)  # check whether every atom has been assigned a type
    teller = -1
    for ffatype, filter in afilters:
        print(ffatype)
        teller = teller + 1
        for iatom, number in enumerate(yaffsystem.numbers):
            if filter(iatom, graph) and not test[iatom]:
                test[iatom] = True
                ffatypes[iatom] = ffatype
                print(ffatype, iatom)
    for iatom, number in enumerate(yaffsystem.numbers):
        if not test[iatom]:
            print('No atom type found for atom %i(%s)' %(iatom, pt[number].symbol))
            print('This atom has neighbors:')
            for neighbor in graph.neighbors[iatom]:
                print('  %i (%s)' %(neighbor, pt[yaffsystem.numbers[neighbor]].symbol))
    return ffatypes


def read_sample(fn,fn_fin):
    #log.set_level(log.silent)
    yaffsystem = YaffSystem.from_file(fn)

    #GaGabond = {(40,40):1.0*angstrom,(40,6):0.5*angstrom,(1,40):0.5*angstrom}

    yaffsystem.detect_bonds()


    # Set atomtypes
    ffatypes = get_atypes(yaffsystem)
    list_atypes = dict(Counter(ffatypes))
    print("\n\n Number of occurrences per atom type:")
    for key in sorted(list_atypes.keys()):
        try:
            print("{:10s}:{:5d}".format(key,list_atypes[key]))
        except ValueError:
            print("An error occured, attention!")

    # Finish
    temp = load_chk(fn)
    dump_chk(fn_fin, {
    'numbers': yaffsystem.numbers,
    'pos': yaffsystem.pos,
    'ffatypes': ffatypes,
    'bonds': yaffsystem.bonds,
    'rvecs': yaffsystem.cell.rvecs,
    'masses': yaffsystem.masses,
    #'hessian': temp['hessian'],
    #'gradient': temp['gradient'],
    })

def main():
    #READ SAMPLE
    read_sample(fn_chk,fn_fin)

if __name__=='__main__':
    main()
