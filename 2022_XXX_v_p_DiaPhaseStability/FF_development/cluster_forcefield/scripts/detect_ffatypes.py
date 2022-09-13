import numpy as np
from collections import Counter

from optparse import OptionParser

from yaff import System, log
from molmod.periodic import periodic
from molmod.molecular_graphs import MolecularGraph
from quickff.tools import set_ffatypes as qff_set_ffatypes

def parse():
    usage = 'python %prog [options] sys'
    descr = 'Define the atomtypes for a given system'
    parser = OptionParser(usage = usage, description = descr)
    parser.add_option(
        '--output', default = None, help = 'Name of the output file'
    )
    parser.add_option(
        '-l', '--level', default = None,
        help = 'Level at which the ffatypes have to be detected '
               '[default = symmetric].'
    )
    parser.add_option(
        '-i', '--inner', action = 'append',
        help = 'Define the inner geometry of the cluster (comma separated)'
    )
    parser.add_option(
        '-o', '--outer', action = 'append',
        help = 'Define the outer terminations of the cluster '
               'These will receive the ffatype Xd_term, with '
               'X the atomic symbol and d the minimal distance '
               'to the inner geometry.'
    )
    parser.add_option(
        '-s', '--suffix', default = None,
        help='Suffix that is pasted after every ffatype. Enforced atom '
        'types (terminations) are not given the suffix, except if it ends '
        'with a \'!\' (which is no part of the suffix).'
    )
    parser.add_option(
        '-q', '--quiet', default = False,
        action = 'store_true', help = 'No yaff output is given'
    )
    options, args = parser.parse_args()
    assert len(args) == 1, 'Exactly one argument expected: input system (got {})'.format(args)
    sys_fn = args[0]
    if options.output is None:
        options.output = sys_fn.split('.')[0] + '_ffatypes.chk'
    if options.level is None:
        options.level = 'symmetric'
    inner = []
    for indices in options.inner:
        inner.extend([int(index) for index in indices.split(',')])
    terms = []
    for index in options.outer:
        assert len(index.split(',')) == 1, 'Only 1 index allowed for outer terms'
        terms.append(int(index))
    if options.quiet:
        log.set_level(0)
    if options.suffix == None:
        options.suffix = ''
    return sys_fn, options.output, options.level, inner, terms, options.suffix

def set_ffatypes(system, level, enforce={}):
    print('Setting ffatypes with level {} and {} atoms to enforce'.format(level, len(enforce)))
    if level in ['low', 'medium', 'high', 'higest']:
        qff_set_ffatypes(system=system, level=level, enforce=enforce)
    elif level == 'symmetric':
        # Initiate enforced atom types
        atypes = [None for i in range(system.natom)]
        for i, number in enumerate(system.numbers):
            if i in enforce.keys():
                atypes[i] = enforce[i]
            elif periodic[number].symbol in enforce.keys():
                atypes[i] = enforce[periodic[number].symbol]
        # Find other atomtypes with symmetric algorithm
        if system.bonds is None:
            system.detect_bonds()
        graph = MolecularGraph(system.bonds, system.numbers)
        # Initial atypes describing the atom and the neighbors and second neighbors
        atype_sets = {}
        work = np.zeros(system.natom)
        for i in range(system.natom):
            if work[i] == 1: continue
            atype = '{}'.format(periodic[system.numbers[i]].symbol)
            for j, neighs in enumerate([system.neighs1, system.neighs2]):
                atype += '_'
                counter = Counter([system.numbers[k] for k in neighs[i]])
                for number in sorted(counter.keys()):
                    if counter[number] == 1:
                        atype += '{}'.format(periodic[number].symbol)
                    else:
                        atype += '{}{}'.format(periodic[number].symbol, counter[number])
            sets = atype_sets.setdefault(atype, [])
            equiv = set([j for j in graph.equivalent_vertices[i] if atypes[j] == None])
            if not equiv in sets:
                sets.append(equiv)
                for j in graph.equivalent_vertices[i]:
                    work[j] = 1
        # Reduce the obtained atomtypes where possible: not needed C_C2H_C2H2 if simply C is enough
        for atype in list(atype_sets.keys()):
            for i in [1, 2]:
                atype_red = '_'.join(atype.split('_')[:i])
                if not np.any([not (atype == atype2) and atype2.startswith(atype_red + '_')
                               for atype2 in atype_sets.keys()]):
                    atype_sets[atype_red] = atype_sets.pop(atype)
                    break
        # Asign atypes taking into account multiple sets can exist for same label
        for atype, sets in atype_sets.items():
            if len(sets) > 1:
                for i, equiv in enumerate(sets):
                    for index in equiv:
                        atypes[index] = '{}_{}'.format(atype, i)
            else:
                for index in sets[0]:
                    atypes[index] = atype
        system.ffatypes = np.array(atypes)
        system._init_derived_ffatypes() 
    else:
        raise ValueError('Invalid level, recieved %s' % level)

def get_termination_atomtypes(system, inner, outer):
    # Only implemented for terminations that are along 1 bond
    # The inner atoms should be inside the inner geometry, the outer outside
    if system.bonds is None:
        system.detect_bonds()
    graph = MolecularGraph(system.bonds, system.numbers)
    enforce = {}
    outer_indices = graph.get_part(outer, inner)
    for i in outer_indices:
        min_dist = min([graph.distances[i, j] for j in range(system.natom) if j not in outer_indices])
        enforce[i] = '{}{}_term'.format(periodic[system.numbers[i]].symbol, min_dist)
    print('Found {} atom in the termination'.format(len(enforce)))
    return enforce

def add_ffatype_suffix(system, suffix, enforce):
    if suffix == '': return
    do_enforce = False
    if suffix[-1] == '!':
        do_enforce = True
        suffix = suffix[:-1]
    if not suffix[0] == '_':
        suffix = '_' + suffix
    atypes = []
    for i in range(system.natom):
        if do_enforce or not i in enforce.keys():
            atypes.append(system.get_ffatype(i) + suffix)
        else:
            atypes.append(system.get_ffatype(i))
    system.ffatypes = None
    system.ffatype_ids = None
    system.ffatypes = np.array(atypes)
    system._init_derived_ffatypes()

def to_file(system, fn):
    print('Printing ffatypes to {}'.format(fn))
    if fn.endswith('.chk'):
        system.to_file(fn)
    elif fn.endswith('.txt'):
        with open(fn, 'w') as f:
            for i in range(system.natom):
                f.write('{}\n'.format(system.get_ffatype(i)))
    else:
        print('Can not write to file {}, expected chk or txt extensions'.format(fn))

def check(system):
    if type(system.ffatypes) == None: return
    for ffatype in system.ffatypes:
        if len(ffatype) > 20:
            print('WARNING: ffatype {} had more than 20 characters. Try to reduce it.'.format(ffatype))

if __name__ == '__main__':
    sys_fn, fn_out, level, inner, terms, suffix = parse()
    print('Reading system from {}'.format(sys_fn))
    system = System.from_file(sys_fn)
    enforce = {}
    for outer in terms:
        enforce.update(get_termination_atomtypes(system, inner, outer))
    set_ffatypes(system, level, enforce)
    add_ffatype_suffix(system, suffix, enforce)
    check(system)
    to_file(system, fn_out)
    
