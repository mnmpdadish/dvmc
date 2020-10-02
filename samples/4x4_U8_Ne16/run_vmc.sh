#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node 18
#QSUB -mpi 18
#QSUB -omp 24
#QSUB -place pack
#QSUB -over false
#PBS -N run_vmc.sh
#PBS -l walltime=00:29:30
cd ${PBS_O_WORKDIR}

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=24

#This example is ready to run on issp sekirei, with module:
module rm intel/18.0.5.274
module add gnuplot/5.2.3 intel/16.0.4.258 intel-mpi/5.1.3.258 intel-mkl/16.0.4.258

#before running this submission, type:
#$ preprocess_dVMC.sh
#to generate all the input files.

#### mVMC groud state calculations:
date
mpijob vmc_new.out namelist.def
date

#### mVMC dynamical Green calculations:
mpijob vmc_new.out namelist_G.def ./output/zqp_opt.dat
date

#to generate the graph, type:
#$ postprocess_dVMC.sh

#to generate the data and gnuplot files:
#$ gnuplot plot_Akw.gp
#$ gnuplot plot_allAkw.gp
#$ gnuplot plot_dos.gp

