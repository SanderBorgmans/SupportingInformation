from yaff import *
import numpy as np

def get_idx(do_return=True):
    sys = System.from_file('init.chk')
    nl_id = np.where(sys.ffatypes == 'N_0_TNA')[0]
    nl_idx = np.where(sys.ffatype_ids == nl_id)[0]

    c4_id = np.where(sys.ffatypes == 'C_H2C2_TNA')[0]
    c4_idx = np.where(sys.ffatype_ids == c4_id)[0]

    if len(c4_idx)>4:
        c4_idx = c4_idx[::len(c4_idx)//4]

    sets = {}
    labels = {}
    for n in range(4):
        #substitute c4_idx by LB,RB,RT,LT (leftbottom, rightbottom, righttop, lefttop)
        label = ''
        if sys.pos[c4_idx[n]][0] < 5*angstrom:
            label+= 'L'
        else:
            label+= 'R'

        if sys.pos[c4_idx[n]][1] < 10*angstrom:
            label+= 'B'
        else:
            label+= 'T'

        print(label)

        sets[label] = tuple(nl_idx[n*4:(n+1)*4])
        labels[label] = c4_idx[n]

    print(labels)
    print(sets)

    # Sort such that the sequence for each vertex is left, bottom, right, top

    for k,v in sets.items():
        c4_i = labels[k]
        vecs = [sys.pos[i] - sys.pos[c4_i] for i in v]
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
    get_idx(do_return=False)
