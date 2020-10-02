#!/bin/bash
#SBATCH --job-name=run_vmc
#SBATCH --time=00:29:00
#SBATCH --account=def-charleb1
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=2000MB

cd $SLURM_SUBMIT_DIR

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=32

#This example is ready to run on cedar in Canada, with module:
# module reset 

#before running this submission, type:
#$ preprocess_dVMC.sh
#to generate all the input files.

#### mVMC groud state calculations:
date
mpirun -np 32 dvmc.out namelist.def
date

#### mVMC dynamical Green calculations:
mpirun -np 32 dvmc.out namelist_G.def ./output/zqp_opt.dat
date

#to generate the graph, type:
#$ postprocess_dVMC.sh

#to generate the data and gnuplot files:
#$ gnuplot plot_Akw.gp
#$ gnuplot plot_allAkw.gp
#$ gnuplot plot_dos.gp


