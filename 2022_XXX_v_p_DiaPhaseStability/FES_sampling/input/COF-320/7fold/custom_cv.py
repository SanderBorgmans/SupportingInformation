#!/usr/bin/env python

import numpy as np
from molmod.units import *
from yaff import *


class CVDistanceProjection(CollectiveVariable):
    '''Compute the length of the vector connecting two centers of masses
       projected onto the xy plane where the cell vectors have been rotated such that
       * the x-axis corresponds to the a direction
       * the y-axis corresponds to the direction perpendicular to a and (a x b)
       * the z-axis corresponds to the direction perpendicular to the ab-plane

       if another plane projection is required, simply rotate the cell.

       cv=norm((r_{COM}^{B}-r_{COM}^{A})[:2])
       and r_{COM} is a vector with centers of mass of groups A and B:

       Note that periodic boundary conditions are NOT taken into account

            * the centers of mass are computed using absolute positions; this is
              most likely the desired behavior
            * the center of mass difference can in principle be periodic, but
              the periodicity is not the same as the periodicity of the system,
              because of the projection on a selected vector
    '''
    def __init__(self, system, groups):
        '''
           **Arguments:**

           system
                An instance of the ``System`` class

           groups
                List of 2 arrays, each array containing atomic indexes
                used to compute one of the centers of mass

        '''
        CollectiveVariable.__init__(self, 'CVDistanceProjection', system)
        # Safety checks
        assert len(groups)==2, "Exactly 2 groups need to be defined"
        assert system.cell.nvec==3, "Only 3D periodic systems are supported"
        # Masses need to be defined in order to compute centers of mass
        if self.system.masses is None:
            self.system.set_standard_masses()
        # Define weights w_i such that difference of centers of mass can be
        # computed as sum_i w_i r_i
        self.weights = np.zeros((system.natom))
        self.weights[groups[0]] = -self.system.masses[groups[0]]/np.sum(self.system.masses[groups[0]])
        self.weights[groups[1]] = self.system.masses[groups[1]]/np.sum(self.system.masses[groups[1]])

    def get_conversion(self):
        return log.length.conversion

    def compute(self, gpos=None, vtens=None):
        '''
           Consider a rotation of the entire system such that the ``a`` vector
           is aligned with the X-axis, the ``b`` vector is in the XY-plane, and
           the ``c`` vector chosen such that a right-handed basis is formed.
           The rotated cell is lower-diagonal in the Yaff notation.

           In this rotated system, it is fairly simple to compute the required
           projection length and derivatives, because the projection length is
           simply the length of the vector using the first two Cartesian
           components. Values obtained in the rotated system are then
           transformed back to the original system.
        '''
        # Compute rotation that makes cell lower diagonal
        _, R = cell_lower(self.system.cell.rvecs)
        # The projected vector of centers of mass difference (aka the
        # collective variable) in the rotated system
        cv_orig = np.sum(self.weights.reshape((-1,1))*self.system.pos, axis=0)
        # Transform back to the original system
        cv = np.dot(R, cv_orig)
        self.value = np.linalg.norm(cv[:2])
        if gpos is not None:
            gpos[:] = 0.0
            gpos[:,:2] = np.outer(self.weights,cv[:2])/self.value
            # Forces (vector) need to be rotated back to original system
            gpos[:] = np.einsum('ij,kj', gpos, R.T)
        if vtens is not None:
            vtens[:] = 0.0
            vtens[:] = np.outer(cv,cv)/self.value
            vtens[2,2] = 0
            # Virial (tensor) needs to be rotated back to original system
            vtens[:] = np.dot(R.T,np.dot(vtens[:],R))
        return self.value

class CVLinComb(CollectiveVariable):
    '''
       A linear combination of CVs:
        cv = w0*ic0 + w1*ic1 + ...
    '''
    def __init__(self, cvs, weights):
        '''
           **Arguments:**

           ics
                A list of CollectiveVariable instances.

           weights
                A list defining the weight of each CollectiveVariable that is
                used when computing the linear combination.


        '''
        assert len(weights)==len(cvs)
        self.cvs = cvs
        self.weights = weights

    def get_conversion(self):
        # Units depend on the particular linear combination of internal
        # coordinates
        return 1.0

    def compute(self, gpos=None, vtens=None):
        self.value = 0.0
        if gpos is not None: gpos[:] = 0.0
        if vtens is not None: vtens[:] = 0.0

        if gpos is not None:
            init_gpos = gpos.copy()

        if vtens is not None:
            init_vtens = vtens.copy()

        for n,cv in enumerate(self.cvs):
            tmp_gpos = init_gpos.copy() if gpos is not None else None
            tmp_vtens = init_vtens.copy() if vtens is not None else None
            self.value += self.weights[n]*cv.compute(gpos=tmp_gpos, vtens=tmp_vtens)

            if gpos is not None:
                gpos += tmp_gpos*self.weights[n]
            if vtens is not None:
                vtens += tmp_vtens*self.weights[n]

        return self.value


def get_cv(ff):
    cv_ohb = CVDistanceProjection(ff.system, [np.array([227]),np.array([368])])
    cv_oht = CVDistanceProjection(ff.system, [np.array([42]),np.array([183])])

    cv_ivl = CVDistanceProjection(ff.system, [np.array([228]),np.array([44])])
    cv_ivr = CVDistanceProjection(ff.system, [np.array([370]),np.array([182])])

    cv_ihb = CVDistanceProjection(ff.system, [np.array([230]),np.array([369])])
    cv_iht = CVDistanceProjection(ff.system, [np.array([41]),np.array([184])])
     
    cv_ovl = CVDistanceProjection(ff.system, [np.array([229]),np.array([43])])
    cv_ovr = CVDistanceProjection(ff.system, [np.array([371]),np.array([185])])

    cv1 = CVLinComb([cv_ohb,cv_oht,cv_ivl,cv_ivr],[0.25]*4)
    cv2 = CVLinComb([cv_ihb,cv_iht,cv_ovl,cv_ovr],[0.25]*4)

    return [cv1,cv2]

def adapt_structure(simulator,ff):
    identity = 'init_structures/{}-{}.chk'.format(simulator.grid,simulator.nr)
    ff.system.to_file(identity)
