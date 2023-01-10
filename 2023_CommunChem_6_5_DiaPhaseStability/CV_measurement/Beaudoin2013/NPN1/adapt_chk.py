import numpy as np
from yaff import *

# This script uses the PBCs to move the relevant atoms (for the CVs) to their appropriate position

system = System.from_file('init.chk')



move_plus_b = [93,85,29,37,21,13,45,5,0,69,61,50,77]
move_plus_a = [90,82,26,18,34,10,2,42,66,58,52,74]

move_min_a = [181,189,125,117,133,141,109,101,165,157,173,151]
move_min_b = [182,190,126,118,134,142,110,102,166,158,174,149]


a,b,c = system.cell.rvecs
a = a/2

system.pos[move_plus_a] += (a-b)
system.pos[move_plus_b] += (a+b)
system.pos[move_min_a] -= (a-b)
system.pos[move_min_b] -= (a+b)

system.to_file('tmp.chk')
