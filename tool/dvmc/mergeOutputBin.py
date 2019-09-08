#!/usr/bin/python
# zlib license:

# Copyright (c) 2019 Maxime Charlebois

# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.

# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:

# 1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
# 2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.


import sys, os
import numpy as np
import re
from ctypes import cdll, c_int, c_float, pointer, c_char_p

full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)


L=W=1
Excitation =  open('excitation.def').read().replace(' ','')
n_exc = int(re.compile('NExcitation([0-9]*)').findall(Excitation)[0])
L = int(re.compile('L([0-9]*)').findall(Excitation)[0])
W = int(re.compile('W([0-9]*)').findall(Excitation)[0])
print n_exc, L, W

def X(ri):
  return ri%W

def Y(ri):
  return (ri/W)%L

def main():
  dirOutput = './output/'
  
  fileIn = []
  n_file = len(sys.argv[:])-1
  n_file2 = 0
  print sys.argv[:]
  print n_file
  
  for nn in range(0,n_file):
    n_file2 +=1
    fileIn.append(sys.argv[1+nn])
  
  print fileIn
  if (n_file <1):
    print 'error: no input files.\nexample\n$ mergeOutputBin.py output/zvo_nCHAm_nAHCm_00*\n'
    exit()
    
  Nsite = L*W
  #use compiled version if present:
  if( not os.path.isfile(pythonPathCode+'/libdvmc_speedup.so')):
    print 'Error: the shared object library "libdvmc_speedup.so" was not found.\nTo obtain it, go in the original directory of this python code (probably: '+ pythonPathCode+'/) and do make.'
    exit(0)
  c_ArrayFloatN = c_float * (n_exc*n_exc*Nsite)

  phys_CA_averaged = c_ArrayFloatN()
  phys_AC_averaged = c_ArrayFloatN()
  phys_CHA_averaged = c_ArrayFloatN()
  phys_AHC_averaged = c_ArrayFloatN()
  c_NExcitation = c_int(n_exc)
  c_W = c_int(W)
  c_L = c_int(L)
  c_ncomp = c_int(n_exc*n_exc*Nsite*4)

  argList = (c_char_p * len(fileIn))()
  argList[:] = fileIn
  
  #Loading and using our home made module:
  lib1 = cdll.LoadLibrary(pythonPathCode+'/libdvmc_speedup.so')
  lib1.mergeOutputBin(len(fileIn), argList, c_NExcitation, c_L, c_W, 
              phys_CA_averaged,phys_AC_averaged,phys_CHA_averaged,phys_AHC_averaged)
  
  print 'after read'
  def convert_c2numpy(phys_averaged):
    numpy_averaged = np.zeros(n_exc*n_exc*Nsite)
    numpy_averaged[:] = phys_averaged[:]
    numpy_reshaped = numpy_averaged.reshape((Nsite,n_exc,n_exc))
    #print test[0,0,:]
    numpyOut = np.zeros([n_exc,n_exc,W,L], dtype='f')
    for rj in range(Nsite):
      rjx = X(rj)
      rjy = Y(rj)
      numpyOut[:,:,rjx,rjy] = numpy_reshaped[rj,:,:]
    return numpyOut

  nCAm_up = convert_c2numpy(phys_CA_averaged)
  nACm_up = convert_c2numpy(phys_AC_averaged)
  nCHAm_up = convert_c2numpy(phys_CHA_averaged)
  nAHCm_up = convert_c2numpy(phys_AHC_averaged)

  print 'after reshape'
    
  print dirOutput+'nCAm'
  print dirOutput+'nACm'
  print dirOutput+'nCHAm'
  print dirOutput+'nAHCm'
  np.save(dirOutput+'nCAm',nCAm_up)
  np.save(dirOutput+'nACm',nACm_up)
  np.save(dirOutput+'nCHAm',nCHAm_up)
  np.save(dirOutput+'nAHCm',nAHCm_up)

if __name__ == "__main__":
   main()
