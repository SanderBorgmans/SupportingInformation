import sys
sys.path.insert(1, '../../../scripts')

import numpy as np

from molmod.unit_cells import UnitCell
from molmod.units import deg

from Topology import Node, WyckoffSet, Topology, split_name, get_suffix, get_diff, get_suffix_difference

def add_edge(top, name, n1, n2):
    prefix, wyckoff_name, number, suffix = split_name(name)
    if type(n1) == str:
        n1 = top.get_node(n1)
    if type(n2) == str:
        n2 = top.get_node(n2)
    frac_pos = (n1.frac_pos + n2.frac_pos)/2
    edge = Node(prefix + wyckoff_name + str(number), 2, frac_pos, neighbors = [n1, n2])
    for neighbor in [n1, n2]:
        neighbor_prefix, neighbor_wyckoff_name, neighbor_number, neighbor_suffix = split_name(neighbor.name)
        unit_neighbor = top.get_unit_node(neighbor.name)
        suffix_difference = get_suffix_difference(edge, neighbor)
        if suffix_difference == '000':
            unit_neighbor.add_neighbor(edge)
        else:
            unit_neighbor.add_neighbor(edge.copy(suffix_difference))
    wyckoff_set = top.wyckoff_edges.setdefault(prefix + wyckoff_name, WyckoffSet(prefix + wyckoff_name))
    wyckoff_set.add_node(edge)

def get_dia_cN(n, do_odd = False):
    # n is the degree of interpenetration
    assert n % 1 == 0 and n > 0, "Expected a positive integer value for n, got {}".format(n)
    if n == 1:
        name = 'dia'
    elif n == 2:
        name = 'dia-c'
    else:
        name = 'dia-c{}'.format(n)
    if n % 2 == 1 or do_odd:
        # Odd degree of interpenetration
        # do_odd can be turned on to have a consistent structure
        nvertices = 4
        nedges = 8
        a = 1.6330 # sqrt(2)/2*2.3094
        b = 1.6330
        c = 2.3094/n
        unit_cell = UnitCell.from_parameters3([a, b, c], [90*deg, 90*deg, 90*deg])
        unit_cell_c1 = UnitCell.from_parameters3([1.6330, 1.6330, 2.3094], [90*deg, 90*deg, 90*deg])

        # Initialize topology
        top = Topology(name, unit_cell, dimension = '3D')

        # Nodes
        frac_pos_c1 = [
                np.array([0.0, 0.75, 0.375]),
                np.array([0.5, 0.75, 0.125]),
                np.array([0.5, 0.25, 0.875]),
                np.array([0.0, 0.25, 0.625])
        ]
        cart_pos_c1 = [unit_cell_c1.to_cartesian(frac_pos) for frac_pos in frac_pos_c1]
        wyckoff_set = WyckoffSet('A')
        for i, cart_pos in enumerate(cart_pos_c1):
            frac_pos = unit_cell.to_fractional(cart_pos)
            frac_pos_unit = frac_pos - np.floor(frac_pos)
            wyckoff_set.add_node(Node('A{}'.format(i+1), 4, frac_pos_unit))
        top.wyckoff_vertices = {'A': wyckoff_set}
        top.nvertices = nvertices

        # Edges
        edges_c1 = [
                ['_A1', 'A1000', 'A4000'],
                ['_A2', 'A10-0', 'A4000'],
                ['_A3', 'A2000', 'A300-'],
                ['_A4', 'A20-0', 'A300-'],
                ['_A5', 'A1000', 'A2000'],
                ['_A6', 'A3000', 'A4000'],
                ['_A7', 'A1+00', 'A2000'],
                ['_A8', 'A3000', 'A4+00']
        ]
        for edge_name, n1, n2 in edges_c1:
            cart_pos_n1 = cart_pos_c1[int(n1[1]) - 1] + unit_cell_c1.to_cartesian(get_diff(n1[-3:]))
            cart_pos_n2 = cart_pos_c1[int(n2[1]) - 1] + unit_cell_c1.to_cartesian(get_diff(n2[-3:]))
            assert abs(np.linalg.norm(cart_pos_n1 - cart_pos_n2) - 1.0) < 1e-3
            cart_pos = (cart_pos_n1 + cart_pos_n2)/2
            frac_pos = unit_cell.to_fractional(cart_pos)
            frac_pos_unit = frac_pos - np.floor(frac_pos)
            diff = frac_pos - frac_pos_unit
            cart_pos_n1 = cart_pos_n1 - unit_cell.to_cartesian(diff)
            cart_pos_n2 = cart_pos_n2 - unit_cell.to_cartesian(diff)
            assert abs(np.linalg.norm(cart_pos_n1 - cart_pos_n2) - 1.0) < 1e-3
            diff_n1 = cart_pos_n1 - unit_cell.to_cartesian(top.wyckoff_vertices['A'].nodes[int(n1[1]) - 1].frac_pos)
            diff_n2 = cart_pos_n2 - unit_cell.to_cartesian(top.wyckoff_vertices['A'].nodes[int(n2[1]) - 1].frac_pos)
            suffix_n1 = get_suffix(np.round(unit_cell.to_fractional(diff_n1), 3))
            suffix_n2 = get_suffix(np.round(unit_cell.to_fractional(diff_n2), 3))
            if not suffix_n1 == '000':
                n1 = 'A{}{}'.format(n1[1], suffix_n1)
            else:
                n1 = 'A{}'.format(n1[1])
            if not suffix_n2 == '000':
                n2 = 'A{}{}'.format(n2[1], suffix_n2)
            else:
                n2 = 'A{}'.format(n2[1])
            add_edge(top, edge_name, n1, n2)
        top.nedges = nedges
    else:
        raise NotImplementedError
    top.check()
    return top

if __name__ == '__main__':
    for i in range(1, 10):
        print(i)
        dia = get_dia_cN(i, do_odd = True)

