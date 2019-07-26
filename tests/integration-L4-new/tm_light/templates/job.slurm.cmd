#!/bin/bash -x
#SBATCH --job-name=JOBNAME
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

rundir=RUNDIR
exe=EXEC
outfile=OUTFILE
infile=INFILE
export QUDA_RESOURCE_PATH=QUDA_RSC_PATH

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

