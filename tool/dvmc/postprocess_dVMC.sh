#!/bin/sh

####  mVMC postprocessing :

#merge different Monte Carlo iteration file output/zvo_nCHAm_nAHCm_0*.bin
#in order to produce the python files:
#nACm.npy, nAHCm.npy, nCAm.npy, nCHAm.npy
#in output directory
mergeOutputBin.py output/zvo_nCHAm_nAHCm_0*.bin

#calculate the full spectrum and generate the files:
#Akw_all.dat, Akw_e_all.dat, Akw_h_all.dat, dos.dat
#in output directory
dvmc_spectrum.py spectrumpara.def output

#copy and adapt the gnuplot template files:
#plot_Akw.gp  plot_allAkw.gp  plot_dos.gp
#in the current directory:
generate_template_gnuplot.py spectrumpara.def

#calulate chemical potential mu and 
#adapt the files:
#plot_Akw.gp  plot_allAkw.gp  plot_dos.gp
#in the current directory.
#And generate the LHB and UHB partial dos files:
#dos_lhb.dat  dos_p.dat  dos_uhb.dat
#in the output directory
dos_analysis.py 

#generate the files 
#Akw.dat, Akw_e.dat, Akw_h.dat
#from the full Akw files:
#Akw_all.dat, Akw_e_all.dat, Akw_h_all.dat, dos.dat
#by selecting the k-points defined in spectrumpara.def
select_kPath.py spectrumpara.def

#after all this, we can plot simply using the commands:
#
#$ gnuplot plot_Akw.gp
#$ gnuplot plot_allAkw.gp
#$ gnuplot plot_dos.gp
#
#you can also visualize the Fermi surface in real time
#by using the program (you need to have gnuplot installed):
#
#$ dvmc_gnuplot
#
#instructions are shown upon typing this command. Only make 
#sense in 2D



