import numpy as np
from yaff import *
# This script uses the PBCs to move the relevant atoms (for the CVs) to their appropriate position

system = System.from_file('init.chk')


move_pbc_a_plus = [130,106]
move_pbc_a_min = [247,271]

move_plus_b = [136,128,48,112,40,56,104,96,32,24,64,120,139,50,42,98
,34,26,66,122,58,114,93,100,132,124,45,37,29,21,61,53,109,117,134,127,39,47,55,31,23,63,119,111,
10,2,8,18,7,15,13,5,16,1,72,74,80,88,79,87,77,85,82,90,71,69,95,103,107,
247,271,
279,131]
move_plus_a = [124,132,135,100]

move_min_a = [265,273,274,241,]
move_min_b = [278,130,106]


a,b,c = system.cell.rvecs
a = a/2


system.pos[move_pbc_a_plus] += 2*a
system.pos[move_pbc_a_min] -= 2*a

system.pos[move_plus_a] += (a-b)
system.pos[move_plus_b] += (a+b)
system.pos[move_min_a] -= (a-b)
system.pos[move_min_b] -= (a+b)

system.to_file('tmp.chk')
