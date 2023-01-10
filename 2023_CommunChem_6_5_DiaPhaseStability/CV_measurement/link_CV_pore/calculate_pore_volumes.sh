#!/bin/bash
#

run_up=999
step_size=100
total_run_time=2999
material=$1

for traj in ${material}/trajs/*.h5; do
    mkdir ${traj%.h5}/
    for n in `seq ${run_up} ${step_size} ${total_run_time}`; do
        c3f.py convert ${traj}:${n} ${traj%.h5}/${n}.cif
        zeo -ha -vol 1.2 1.2 50000 ${traj%.h5}/${n}.cif
    done
done
