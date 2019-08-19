# mVMC: dynamical Green

Author: Maxime Charlebois

This code is based on the original mVMC open source package 
[source](https://github.com/issp-center-dev/mVMC) 
and [arXiv:1711.11418](https://arxiv.org/abs/1711.11418).
Please refer to the documentation in the open source package
as this will not be covered here.

It implements a new feature where you can calculate the frequency 
dependent Green function (option NVMCCalMode = 3), hence the
spectral weight. For now there is a number of condition on this 
new calculation mode. It is restricted to (for now):
- 1D or 2D model (easy to generalize to 3D, to be done)
- real calculation
- OrbitalAntiParallel (not OrbitalGeneral or OrbitalParallel)


# Usage

Detailed usage is not covered here. Instead many examples can be found
in the "./samples/spectrum/" subdirectory. In summary, the code requires
a converged ground state obtained by the "NVMCCalMode 0" option. From 
this the spectral weight can be calculated with option "NVMCCalMode 3".

This mode requires a new input file "excitation.def". This can be generated
in an automated way from the python code:
"./tool/dynamicalGreenProcessing/makeExcitation.py"
together with the input file "spectrumpara.def". In this last file, the
values "dr1_x, dr1_y, dr2_x, dr2_y", defines the difference in position
of the charge excitation relatif to site "i" or "j" of the operator 
"c^dagger_i" or  "c_j" for the Green function <n|c^dagger_i cj|m>. 
This is technical (may requires more information soon, or reference 
to the arxiv paper). For example:

"
dr1_x       -1:1
dr1_y       -1:1
dr2_x       -1:2
dr2_y       -1:2
"

will generate 244 different charge charge excitations in the file
"excitation.def". Once the file is generated, the code in this 
mode can be run through this command:

$ vmc_new.out namelist_G.def ./output/zqp_opt.dat
or
$ mpirun -np $nproc vmc_new.out namelist_G.def ./output/zqp_opt.dat

for example. The files will be generated in the "./output/" subdirectory
with names "zvo_nCHAm_nAHCm_###.bin". These can then be treated to produce
a graph of the density of states with the functions in:
"./tool/dynamicalGreenProcessing/*"

In order

$ mergeOutput.py output/zvo_nCHAm_nAHCm_0*   # to average the different iterations in 4 numpy data array.
$ print_spectrum.py                          # to output a matrix of A(k,w) in the files Akw.dat, Akw_e.dat, Akw_h.dat

This last function read the parameter from many files including the frequency
axis parameters in "spectrumpara.def".

After that, the spectral density A(k,w) can be plotted through 
the simple gnuplot command:
> plot 'output/Akw.dat' matrix notitle w image

A better (prettier) version of the same gnuplot graph can be generated

$ generate_template_gnuplot.py
$ gnuplot plot_singleAkw.gp
$ gnuplot plot_allAkw.gp

although, some generic parameters might need to be adjusted to obtain a 
good result. For example, the an offset must be applied on the frequency
axis so that the Fermi level (omega = 0) is placed manually between the
electron minimum and hole maximum density of state.  

Gnuplot is the software of choice here because it can embed an image
of the A(k,w), resulting in smaller file, faster to process and easier
to scroll by in a PDF file.


# Samples description:

"./samples/spectrum/chain2_U8"  -- trivial example for comparison with analytical results
more to come


# Requirements

Have this new mVMC code compiled first (follow original mVMC 
documentation). Put it in a directory available in the PATH environment 
variable (usually ~/bin/.). Name this executable "vmc_new.out". 
Also copy the files contained in 
./tool/dynamicalGreenProcessing/*
to a path contained in the PATH environment variable.

The libraries required to compile the code are:
- mpi
- lapack or MKL
- blas or MKL






############################################################################
# the rest of the README is the one contained in the original mVMC package #
############################################################################


# mVMC

A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method

### What is mVMC ?

mVMC (many-variable Variational Monte Carlo method)
is a software for performing the highly-accurate 
variational Monte Carlo calculations
with the simple and flexible user interface.
mVMC also supports the large-scale parallelization.
For the conventional models in strongly correlated electron systems such as the Hubbard model, the Heisenberg model, and the Kondo-lattice model,
users can perform the calculation by preparing the one input files whose length is shorter than ten lines.  
By using the same input file, users can perform the exact diagonalization through [HPhi](https://github.com/QLMS/HPhi/releases).
Thus, it is easy to examine the accuracy of the variational calculation for small system sizes
and to perform the calculations 
for large system sizes that can not be treated 
by the exact diagonalization.
A broad spectrum of users including experimental scientists is cordially welcome.


### Methods

many-variable variational Monte Carlo method


### Target models

Hubbard model, Heisenberg model, Kondo lattice model, multi-orbital Hubbard model

### Available physical quantities

specific heat, susceptibility, ground state energy, structure factors


## Requirement

- C compiler (intel, Fujitsu, GNU, etc. ) 
- ScaLAPACK library (intel MKL, Fujitsu, ATLAS, etc.) 
- MPI library

## Install

You can install mVMC and also get a manual for mVMC from a [release note](https://github.com/issp-center-dev/mVMC/releases).


## Licence

GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)). 

The mVMC package is developed based on the [mVMC-mini](https://github.com/fiber-miniapp/mVMC-mini) program. The license of mVMC-mini is "The BSD 3-Clause License".

We would appreciate if you cite the following article in your research with mVMC:  
mVMC - Open-source software for many-variable variational Monte Carlo method, Takahiro Misawa, Satoshi Morita, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Yuichi Motoyama, Kota Ido, Takahiro Ohgoe, Masatoshi Imada, Takeo Kato, [arXiv:1711.11418](https://arxiv.org/abs/1711.11418).

## Tutorials

Lecture notes and sample scripts used in Hands-on
are available at [mVMC-tutorial](https://github.com/issp-center-dev/mVMC-tutorial)

## Authors

Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Yuichi Motoyama, Mitsuaki Kawamura, Kazuyoshi Yoshimi, Takeo Kato, Masatoshi Imada.
