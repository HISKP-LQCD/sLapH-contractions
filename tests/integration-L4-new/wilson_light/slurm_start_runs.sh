#!/bin/bash

conf_start=1000
conf_step=4
conf_end=1000

rv_start=0
rv_end=9

for i in $( seq ${conf_start} ${conf_step} ${conf_end} ); do
  j=`printf %04d $i`
  for rv in $( seq ${rv_start} ${rv_end} ); do
    echo "starting config $i rv $rv"
    rnd2=$( printf %02d ${rv} )
    wdir=cnfg${j}/rnd_vec_${rnd2}

    pushd ${wdir}
    jscr=job.slurm.${j}_${rnd2}.cmd
    sbatch ${jscr}
    popd

  done
done
