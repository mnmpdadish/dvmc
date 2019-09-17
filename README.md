# mVMC: dynamical Green

Author: Maxime Charlebois

This code is based on the original mVMC open source package 
[source](https://github.com/issp-center-dev/mVMC) 
and [arXiv:1711.11418](https://arxiv.org/abs/1711.11418).
Please refer to the documentation in the open source package
as this will not be covered here.

It implements a new feature where you can calculate the frequency 
dependent Green function (option NVMCCalMode = 3). 
For now there is a number of condition on this 
new calculation mode. It is restricted to (for now):
- 1D or 2D model (easy to generalize to 3D, to be done)
- real calculation
- OrbitalAntiParallel (not OrbitalGeneral or OrbitalParallel)
- Hubbard model

# Usage

Detailed usage is not covered here. Instead many examples can be found
in the "./samples/spectrum/" subdirectory. The easiest way to understand
how to use it is to run these examples. Go there to read the README
and run the few examples. However, the code must be installed first
(see next section)

# Installation:

To install everything, this new version of mVMC must be compiled. 
This can be done using cmake. In the currect (root) directory, do:

$ mkdir build && cd build
$ cmake -DCONFIG=intel -DCMAKE_INSTALL_PREFIX:PATH=~ ..
$ make 
$ make install

This should install everything with correct linking in the bin directory 
in your home directory. We use a custom binary directory because on a 
remote cluster, it might be impossible or prohibited to install using 
sudo. For this to work, you need to create the ~/bin/ and make this 
directory execuable by adding this line to the ~/.bashrc 
(if this file does not exist, create it):

export PATH=$PATH:$HOME/bin/

after that, type the command line:

$ source ~/.bashrc 

Note, the -DCONFIG=intel option requires that you have icc, ifort, etc.
If you want to use gcc, gfortran, etc. instead, use the line instead:
$ cmake -DCONFIG=gcc -DCMAKE_INSTALL_PREFIX:PATH=~ ..

# Requirements

The libraries required to compile the code are:
- mpi
- lapack or MKL
- blas or MKL

It is recommended to consult the original mVMC documentation to learn
more about dependencies and parameters names and idea. Note however that
not all the case in mVMC are covered by the present code, as stated above.


# Details:

Two more README are encourage to read after this one:

"./samples/spectrum/README" - to learn about the code usage (with functionnal examples).
"./tool/dvmc/README"        - to learn about the tools that can be used to analyse the data.

This complete the documentation.


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
