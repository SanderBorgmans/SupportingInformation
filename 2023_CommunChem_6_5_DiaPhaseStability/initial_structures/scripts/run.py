import sys
sys.path.insert(1, '../input/topologies/')

from dia_cN import get_dia_cN

# For SBU and GeometryConstructor modules: see follow-up paper
from SBU import SBU
from Construct import GeometryConstructor
SBU.sbu_path = '../input/sbus'


# Topologies
tops = {'dia-c{}'.format(i): get_dia_cN(i, do_odd = True) for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]}
# SBUs
tam = SBU.load('TAM', '_PrimaryAmine_Imine')
tna = SBU.load('TNA', '_Nitroso_Azodioxy')
tnm = SBU.load('TNM', '_Nitroso_Azodioxy')
bdc = SBU.load('BDC', '_Aldehyde_Imine')
bpdc = SBU.load('BPDC', '_Aldehyde_Imine')

# Load FF parameters
for sbu in [tam, tna, tnm, bdc, bpdc]:
    fn_cluster = ['pars_yaff.txt', 'pars_ei.txt', 'pars_mm3.txt']
    fn_uff = ['pars_cov_uff.txt', 'pars_ei.txt', 'pars_lj_uff.txt']
    sbu.load_parameters(fns = fn_cluster, uff_fns = fn_uff)

# Structures
building_blocks = {
        'COF-300': [tam, bdc],
        'COF-320': [tam, bpdc],
        'NPN-1': [tnm, None],
        'NPN-3': [tna, None]
        }
for top_name, top in tops.items():
    for cof, (node4, node2) in building_blocks.items():
        if not cof == 'COF-320' and int(top_name.split('-c')[1]) > 8: continue
        geom = GeometryConstructor(top, [node4, node2])
        geom.get_nucleation_structure()
        geom.output('../output/{}/{}_{}_{}'.format(cof, top_name, node4.name, str(node2.name)))

