#!/bin/bash 
#$ -S /bin/bash
#$ -l h_rt=120:00:00
#$ -j y
#$ -m n
#$ -cwd
#$ -pe openmpi-fill 28

. ~/.bashrc
module load uge compiler/intel-17.0.4 openmpi/1.10.6/intel-17.0.4.lp mpich/3.2/intel-17.0.4.lp

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

#do the preprocessing typing:
#$ preprocessVMC.sh

#### mVMC groud state calculations:
date
mpirun -np 28 vmc_new.out namelist.def
#mpirun -np 28 vmc_new.out namelist.def ./output/zqp_opt.dat  # to add some iterations
date

#### mVMC dynamical Green calculations:
makeExcitation.py
mpirun -np 28 vmc_new.out namelist_G.def ./output/zqp_opt.dat
date


