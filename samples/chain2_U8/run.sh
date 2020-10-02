#!/bin/bash 

python3 analytical_mVMC.py
vmc_new.out namelist.def output/zqp_opt.dat
python3 compare.py

