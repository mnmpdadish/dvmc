#!/bin/sh
#QSUB -queue F144cpu
#QSUB -node 144
#QSUB -mpi 144
#QSUB -omp 24
#QSUB -place pack
#QSUB -over false
#PBS -N run_vmc.sh
#PBS -l walltime=23:55:00
cd ${PBS_O_WORKDIR}

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=24

#init_groundstate.py opt.in
#mkdir output
#mv zqp_opt.dat output/.

#### mVMC groud state calculations:
date
mpijob dvmc.out namelist.def ./output/zqp_opt.dat
date
mpijob dvmc.out namelist_G.def ./output/zqp_opt.dat
date

# careful with the postprocessing.sh. This cannot be
# applied automatically, this example requires more attention.
# use function filter_bin.sh to show which file does not respect
# the sum rule. Keep only the file that respect the sum rule
# (0.96 < int A(k,w) < 1.01), i.e. keep only the file that does
# not produce a warning by the function filter_bin.py.
# This process need to be refined at some point.

