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

import numpy as np
from numpy import linalg as la
import sys, os, re

def main():
  
  S_CA_k = FFT('output/nCAm.npy')
  S_AC_k = FFT('output/nACm.npy')
  H_CA_k = FFT('output/nCHAm.npy')
  H_AC_k = FFT('output/nAHCm.npy')
  np.save('output/nCAm_k.npy',S_CA_k)
  np.save('output/nACm_k.npy',S_AC_k)
  np.save('output/nCHAm_k.npy',H_CA_k)
  np.save('output/nAHCm_k.npy',H_AC_k)
  
  exit()


def FFT(dataFileName):
  print 'fft of '+ dataFileName + '.'
  sys.stdout.flush()
  data_up = np.load(dataFileName)
  
  
  #print data_up
  data_k  = np.fft.fft2(data_up) # fft2 only on the last 2 indices
  
  print data_k
  print 
  for ii in range(3):
   for jj in range(3):
    for kk in range(3):
     print '% 3.5f '% data_up[ii,jj,0,kk]
  return data_k


if __name__ == "__main__":
   main()

