#!/bin/sh
#QSUB -queue L18acc
#QSUB -node 18
#QSUB -mpi 216
#QSUB -omp 2
#QSUB -place pack
#QSUB -over false
#PBS -N run_vmc.sh
#PBS -l walltime=72:00:00
cd ${PBS_O_WORKDIR}

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=2

####  mVMC preparations :
vmcdry_new.out StdFace.def 

date
sed -i "s/NExUpdatePath  0/NExUpdatePath  1/"     modpara.def
sed -i "s/NDataQtySmp    1/NDataQtySmp    10/"    modpara.def

cp namelist.def namelist_G.def
cp modpara.def modpara_G.def
sed -i "s/CDataFileHead  zvo/CDataFileHead  zbo/"           modpara.def
sed -i "s/NVMCCalMode    0/NVMCCalMode    3/"               modpara_G.def
sed -i "s/ModPara  modpara.def/ModPara  modpara_G.def/"     namelist_G.def
echo "      Excitation  excitation.def" >> namelist_G.def

#### mVMC groud state calculations:
date
#mpirun -np 28 vmc_new.out namelist.def
#mpirun -np 28 vmc_new.out namelist.def ./output/zqp_opt.dat  # to add some iterations
mpijob vmc_new.out namelist.def 
date

#### mVMC dynamical Green calculations:
makeExcitation.py
mpijob vmc_new.out namelist_G.def ./output/zqp_opt.dat
#mpirun -np 28 vmc_new.out namelist_G.def ./output/zqp_opt.dat
date

