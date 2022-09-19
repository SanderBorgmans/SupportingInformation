from yaff import *
import numpy as np
import sys
import itertools

log.set_level(0)


def get_ICs(system, ics):
    num = 0
    ICl = InternalCoordinateList(DeltaList(system))
    for ic in ics:
        ICl.add_ic(ic)
        num+=1
    ICl.dlist.forward()
    ICl.forward()

    values = np.zeros(num)
    for n in range(num):
        values[n] = ICl.ictab[n][-2]

    return values

def get_bond_dist(system,ref,probes):
    bonds = []
    if not isinstance(ref,list):
        ref = [ref]

    for bond in itertools.product(ref,probes):
        bonds.append(Bond(*bond))

    return np.array(get_ICs(system, bonds))


def get_idx(do_return=True, fname='init.chk'):
    sys = System.from_file(fname)
    rvecs = sys.cell.rvecs
    nl_id = np.where(sys.ffatypes == 'N_L')[0]
    nl_idx = np.where(sys.ffatype_ids == nl_id)[0]

    c4_id = np.where(sys.ffatypes == 'C4')[0]
    c4_idx = np.where(sys.ffatype_ids == c4_id)[0]

    sets = {}
    labels = {}

    xm = np.average(sys.pos[c4_idx,0],axis=0)
    ym = np.average(sys.pos[c4_idx,1],axis=0)
    #print(xm,ym)


    for n in range(4):
        #substitute c4_idx by LB,RB,RT,LT (leftbottom, rightbottom, righttop, lefttop)
        label = ''
        if sys.pos[c4_idx[n]][0] < xm:
            label+= 'L'
        else:
            label+= 'R'

        if sys.pos[c4_idx[n]][1] < ym:
            label+= 'B'
        else:
            label+= 'T'

        # Get four closest nl_idx (use InternalCoordinate to account for MIC)
        dist = get_bond_dist(sys,c4_idx[n],nl_idx)
        closest_nl_idx = nl_idx[np.argsort(dist)[:4]]
        sets[label] = tuple(closest_nl_idx)
        labels[label] = c4_idx[n]
    print(labels)
    print(sets)

    # Sort such that the sequence for each vertex is left, bottom, right, top
    # For each c4, move the n_l that are a periodic image

    for k,v in sets.items():
        c4_i = labels[k]
        vecs = []
        for i in v:
            vec = sys.pos[i] - sys.pos[c4_i]
            if np.linalg.norm(vec - np.sign(vec[0])*rvecs[0]) < np.linalg.norm(vec):
                vec -= np.sign(vec[0])*rvecs[0]
                #sys.pos[i] -= np.sign(vec[0])*rvecs[0]
            if np.linalg.norm(vec - np.sign(vec[1])*rvecs[1]) < np.linalg.norm(vec):
                vec -= np.sign(vec[1])*rvecs[1]
                #sys.pos[i] -= np.sign(vec[1])*rvecs[1]
            vecs.append(vec)

        left_idx = np.argmin([vec[0] for vec in vecs])
        right_idx = np.argmax([vec[0] for vec in vecs])

        top_idx = np.argmax([vec[1] for vec in vecs])
        bottom_idx = np.argmin([vec[1] for vec in vecs])

        sets[k] = (v[left_idx],v[bottom_idx],v[right_idx],v[top_idx])


    ohb = (sets['LB'][1], sets['RB'][1])
    oht = (sets['LT'][3], sets['RT'][3])

    ivl = (sets['LB'][2], sets['LT'][2])
    ivr = (sets['RB'][0], sets['RT'][0])

    ihb = (sets['LB'][3], sets['RB'][3])
    iht = (sets['LT'][1], sets['RT'][1])

    ovl = (sets['LB'][0], sets['LT'][0])
    ovr = (sets['RB'][2], sets['RT'][2])



    print('cv_ohb = CVDistanceProjection(ff.system, [np.array([{}]),np.array([{}])]) '.format(ohb[0],ohb[1]))
    print('cv_oht = CVDistanceProjection(ff.system, [np.array([{}]),np.array([{}])]) '.format(oht[0],oht[1]))

    print('cv_ivl = CVDistanceProjection(ff.system, [np.array([{}]),np.array([{}])]) '.format(ivl[0],ivl[1]))
    print('cv_ivr = CVDistanceProjection(ff.system, [np.array([{}]),np.array([{}])]) '.format(ivr[0],ivr[1]))

    print('cv_ihb = CVDistanceProjection(ff.system, [np.array([{}]),np.array([{}])]) '.format(ihb[0],ihb[1]))
    print('cv_iht = CVDistanceProjection(ff.system, [np.array([{}]),np.array([{}])]) '.format(iht[0],iht[1]))

    print('cv_ovl = CVDistanceProjection(ff.system, [np.array([{}]),np.array([{}])]) '.format(ovl[0],ovl[1]))
    print('cv_ovr = CVDistanceProjection(ff.system, [np.array([{}]),np.array([{}])]) '.format(ovr[0],ovr[1]))

    idx = [ohb,oht,ivl,ivr,ihb,iht,ovl,ovr]

    if do_return:
        return idx


if __name__=='__main__':
    if len(sys.argv)>1:
        get_idx(do_return=False,fname=sys.argv[1])
