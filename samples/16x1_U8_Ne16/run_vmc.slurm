#!/bin/bash
#SBATCH --job-name=run_vmc
#SBATCH --time=00:29:00
#SBATCH --account=def-charleb1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem-per-cpu=2000MB

cd $SLURM_SUBMIT_DIR

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#This example is ready to run on cedar in Canada, with module:
# module add cmake/3.16.3 python/2.7.14 scipy-stack/2019a icc/.2017.5.239 gcccore/.6.4.0 ifort/.2017.5.239 imkl/2017.4.239 intel/2017.5 boost/1.65.1 impi/2017.4.239

#before running this submission, type:
#$ preprocess_dVMC.sh
#to generate all the input files.

#### mVMC groud state calculations:
date
srun --cpu_bind=cores dvmc.out namelist.def
date

#### mVMC dynamical Green calculations:
srun --cpu_bind=cores dvmc.out namelist_G.def ./output/zqp_opt.dat
date

#to generate the graph, type:
#$ preprocess_dVMC.sh

#to generate the data and gnuplot files:
#$ gnuplot plot_Akw.gp
#$ gnuplot plot_allAkw.gp
#$ gnuplot plot_dos.gp

