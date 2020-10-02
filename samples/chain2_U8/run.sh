#!/bin/bash 

python3 analytical_mVMC.py
dvmc.out namelist.def output/zqp_opt.dat
python3 compare.py

