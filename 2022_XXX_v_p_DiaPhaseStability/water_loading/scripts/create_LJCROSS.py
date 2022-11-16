import itertools
import numpy as np


header = """# Lennard-Jones potential
# ======================
#
# Mathematical form:
# E_LJCROSS = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
#
LJCROSS:UNIT SIGMA A
LJCROSS:UNIT EPSILON kcalmol

LJCROSS:SCALE 1 0.0
LJCROSS:SCALE 2 0.0
LJCROSS:SCALE 3 1.0
#------------------------------------------------------------------------------------------
# KEY        label0                label1                SIGMA            EPSILON          
#------------------------------------------------------------------------------------------
"""



parameters_H2O = {'TO':(3.164400000000000,   0.18520764558814762)} # sigma, epsilon
parameters_framework = {
  'C_C4_C8_TAM'   : (3.4308509636e+00, 1.0500000000e-01), 
  'C_C2N_H2C3_TAM': (3.4308509636e+00, 1.0500000000e-01), 
  'C_HC2_HC2N_TAM': (3.4308509636e+00, 1.0500000000e-01), 
  'C_HC2_HC3_TAM' : (3.4308509636e+00, 1.0500000000e-01), 
  'C_C3_H2C5_TAM' : (3.4308509636e+00, 1.0500000000e-01), 
  'H1_C_C2_TAM'   : (2.5711337006e+00, 4.4000000000e-02), 
  'H0_C_C2_TAM'   : (2.5711337006e+00, 4.4000000000e-02), 
  'N_C2_HC3_TAM'  : (3.2606893084e+00, 6.9000000000e-02), 
  'C_C3_H3C2N_BDC': (3.4308509636e+00, 1.0500000000e-01), 
  'C_HC2_HC3_BDC' : (3.4308509636e+00, 1.0500000000e-01), 
  'H_C_C2_BDC'    : (2.5711337006e+00, 4.4000000000e-02), 
  'C_HCN_C3_BDC'  : (3.4308509636e+00, 1.0500000000e-01), 
  'H_C_CN_BDC'    : (2.5711337006e+00, 4.4000000000e-02), 
}


fn = 'pars_ljcross.txt'
with open(fn,'w') as f:
    f.write(header)
    for p1,v1 in parameters_H2O.items():
        for p2,v2 in parameters_framework.items():
            sigma = 0.5*(v1[0]+v2[0]) #
            epsilon = np.sqrt(v1[1]*v2[1]) #
            f.write("LJCROSS:PARS {:21s} {:21s} {:16s} {:16s}\n".format(p1,p2,str(sigma),str(epsilon)))
    f.write('\n')
