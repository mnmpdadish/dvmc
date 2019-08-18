#!/bin/sh
#QSUB -queue F36cpu
#QSUB -node 36
#QSUB -mpi 144
#QSUB -omp 6
#QSUB -place pack
#QSUB -over false
#PBS -N run_vmc240.sh
#PBS -l walltime=23:30:00
cd ${PBS_O_WORKDIR}

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=6

#### mVMC groud state calculations:
date
#mpijob vmc_new.out namelist.def 
mpijob vmc_new.out namelist.def ./output/zqp_opt.dat
date

#### mVMC dynamical Green calculations:
mpijob vmc_new.out namelist_G.def ./output/zqp_opt.dat
date

