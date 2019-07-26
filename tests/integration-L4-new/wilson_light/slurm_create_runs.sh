#!/bin/bash

conf_start=1000
conf_step=4
conf_end=1000

rv_start=0
rv_end=9

seeds=( 12345 98765 34567 13579 86420 67891 54321 19876 97531 56798  )

flavour="wilson_light"

RUNDIR="/hiskp4/bartek/peram_generation/test/test4x4x4x4/${flavour}\/cnfg"
EVDIR="/hiskp4/eigensystems/test/test4x4x4x4/nev_16"
GCONFBASE="/hiskp4/gauges/test/test4x4x4x4/conf"
EXEC="/qbigwork/bartek/build/bleeding_edge/kepler/peram_gen_multigpu.tmLQCD.etmc.quda_develop/main/main"
JOBNAME="test4x4x4x4_${flavour}"
QUDA_RSC_PATH="/qbigwork/bartek/quda_resources/kepler_405d5bf1ac9cdbccbc11ac957e07d822065ac36e"

for i in $( seq ${conf_start} ${conf_step} ${conf_end} ); do
  echo "creating config $i"
  j=`printf %04d $i`

  mkdir -p cnfg${j}/
  mkdir -p cnfg${j}/outputs
  
  for rv in $( seq ${rv_start} ${rv_end} ); do
    seed=${seeds[${rv}]}
    rnd2=$( printf %02d ${rv} )
    wdir=cnfg${j}/rnd_vec_${rnd2}
    mkdir -p ${wdir}

    ifile=${wdir}/invert.input
  
    cp templates/invert.input ${ifile}
    sed -i "s@NSTORE@${i}@g" ${ifile}
    sed -i "s@GCONFBASE@${GCONFBASE}@g" ${ifile}
    
    laphin=LapH_${j}_${rnd2}.in
    jscr=${wdir}/job.slurm.${j}_${rnd2}.cmd
    outfile="../outputs/run_${j}_${rnd2}.out"
  
    cp templates/job.slurm.cmd ${jscr}
    sed -i "s@RUNDIR@${RUNDIR}${j}/rnd_vec_${rnd2}@g" ${jscr}
    sed -i "s@JOBNAME@${JOBNAME}_${j}_${rv}@g" ${jscr}
    sed -i "s@INFILE@${laphin}@g" ${jscr}
    sed -i "s@OUTFILE@${outfile}@g" ${jscr}
    sed -i "s@EXEC@${EXEC}@g" ${jscr}
    sed -i "s@QUDA_RSC_PATH@${QUDA_RSC_PATH}@g" ${jscr}
  
    cp templates/LapH.in ${wdir}/${laphin}
    sed -i "s@NSTORE@${i}@g" ${wdir}/${laphin} 
    sed -i "s@NB_RND@${rv}@g" ${wdir}/${laphin}
    sed -i "s@SEED@${seed}@g" ${wdir}/${laphin}
    sed -i "s@EVDIR@${EVDIR}@g" ${wdir}/${laphin}

  done
done
