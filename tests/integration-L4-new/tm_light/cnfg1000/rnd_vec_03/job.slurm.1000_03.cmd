#!/bin/bash -x
#SBATCH --job-name=test4x4x4x4_tm_light_1000_3
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=bartosz_kostrzewa@fastmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem_bind=verbose
#SBATCH --time=24:00:00
#SBATCH --mem=64200
#SBATCH --gres=gpu:kepler:1
#SBATCH --partition=kepler
#SBATCH --reservation=testing

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/qbigwork/bartek/libs/bleeding_edge/kepler/quda_develop/lib:/opt/cuda/lib64

rundir=/hiskp4/bartek/peram_generation/test/test4x4x4x4/tm_light/cnfg1000/rnd_vec_03
exe=/qbigwork/bartek/build/bleeding_edge/kepler/peram_gen_multigpu.tmLQCD.etmc.quda_develop/main/main
outfile=../outputs/run_1000_03.out
infile=LapH_1000_03.in
export QUDA_RESOURCE_PATH=/qbigwork/bartek/quda_resources/kepler_405d5bf1ac9cdbccbc11ac957e07d822065ac36e

if [ ! -d ${QUDA_RESOURCE_PATH} ]; then
  mkdir -p ${QUDA_RESOURCE_PATH}
fi

cd ${rundir}
date > ${outfile}
QUDA_RESOURCE_PATH=${QUDA_RESOURCE_PATH} OMP_NUM_THREADS=2 \
  QUDA_ENABLE_GDR=1 QUDA_ENABLE_P2P=1 QUDA_ENABLE_TUNING=1 \
  QUDA_ENABLE_DEVICE_MEMORY_POOL=0 \
  srun ${exe} -LapHsin ${infile} | tee -a ${outfile}

date >> ${outfile}

