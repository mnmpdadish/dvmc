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
mpirun -np 28 vmc_new.out namelist.def
#mpirun -np 28 vmc_new.out namelist.def ./output/zqp_opt.dat  # to add some iterations
date

#### mVMC dynamical Green calculations:
makeExcitation.py
mpirun -np 28 vmc_new.out namelist_G.def ./output/zqp_opt.dat
date


