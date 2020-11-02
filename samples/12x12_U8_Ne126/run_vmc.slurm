#!/bin/bash
#SBATCH --job-name=run_vmc
#SBATCH --time=23:59:00
#SBATCH --account=def-charleb1
#SBATCH -n 108
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

# careful with the postprocessing.sh. This cannot be 
# applied automatically, this example requires more attention.
# use function filter_bin.sh to show which file does not respect
# the sum rule. Keep only the file that respect the sum rule 
# (0.96 < int A(k,w) < 1.01), i.e. keep only the file that does 
# not produce a warning by the function filter_bin.py.
# This process need to be refined at some point. 
