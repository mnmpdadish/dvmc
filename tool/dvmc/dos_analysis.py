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
from ctypes import cdll, c_int, c_double

#import matplotlib
#import matplotlib.pyplot as plt


StdFace = open('StdFace.def').read().replace(' ','')
U=0.
L=W=1
Ne=0
if(StdFace.find('L=')>=0): L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
if(StdFace.find('W=')>=0): W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])
if(StdFace.find('U=')>=0): U = float(re.compile('U=([0-9.]*)').findall(StdFace)[0])
if(StdFace.find('nelec=')>=0): Ne = float(re.compile('nelec=([0-9.]*)').findall(StdFace)[0])
N = W*L

w_min_p=0.0
w_max_p=0.0
outputDir='output/'

if (len(sys.argv)==1):
  print 'no energy window.'
elif (len(sys.argv)==2):
  print 'no energy window.'
  outputDir=sys.argv[1]+'/'
elif (len(sys.argv)==3):
  w_min_p=float(sys.argv[1])
  w_max_p=float(sys.argv[2])
else: 
  #print("example (w_min = 1.2, w_max = 1.4):\n$ dos_partial.py 1.2 1.4")
  #print("example (w_window_min=1.2, w_window_max=1.3):\n$ dos_partial.py 1.2 1.3")
  print("example:\n$ dos_partial.py\nor\n$ dos_partial.py output3\nor\n$ dos_partial.py 1.2 1.3\nIn the last example, the two values are frequencies that define an energy windows (to calculate and plot partial dos)")
  sys.exit()

def main():

  file_dos   = open(outputDir+'dos.dat').read()
  file_dos_p = open(outputDir+'dos_p.dat','w')
  file_dos_uhb = open(outputDir+'dos_uhb.dat','w')
  file_dos_lhb = open(outputDir+'dos_lhb.dat','w')
  
  spectrumpara = open('spectrumpara.def').read()
  
  NN   = 0
  norm = 0.0  
  sum8 = np.zeros(8,np.double)
  
  window = 0.0
  LHB = 0.0
  UHB = 0.0
  mu = 0.0
  
  Nw = 0
  for line in  file_dos.split('\n'):
    if len(line)>0:
      term = line.split()
      freq = float(term[0])
      norm += float(term[1])
      Nw+=1

  cumulative_dos = np.zeros(Nw)
  frequencies = np.zeros(Nw)
  cumulated_dos = 0.0

  iw = 0
  for line in  file_dos.split('\n'):
    if len(line)>0:
      term = line.split()
      frequencies[iw] = float(term[0])
      cumulative_dos[iw] = cumulated_dos
      cumulated_dos += float(term[1])/norm
      iw+=1
  
  print frequencies
  print (frequencies[-1]-frequencies[0])/len(frequencies)
      
  mu = np.interp(0.5*(float(Ne)/float(N)),cumulative_dos,frequencies)
  
  N_mu = np.interp(mu,frequencies,range(len(frequencies)))
  
  freq_min =-1000.0
  freq_min =-1000.0
  for freq in frequencies:
    if len(line)>0:
      term = line.split()
      freq = float(term[0])
      
      if(freq<(w_min_p+mu)) :
        freq_min = freq      
      if(freq<=(w_max_p+mu)) :
        freq_max = freq
         
        
      #else:
  
  file_dos_p.write('% 7.6f   % 7.6f  \n'%((w_min_p+mu), 0.0))
  file_dos_uhb.write('% 7.6f   % 7.6f  \n'%(mu, 0.0))
  file_dos_lhb.write('% 7.6f   % 7.6f  \n'%(-1000.0, 0.0))
  
  for line in  file_dos.split('\n'):
    if len(line)>0:
      term = line.split()
      freq = float(term[0])
      
      if(freq>=(w_min_p+mu) and freq<=(w_max_p+mu) ):
        #print freq
        #file_dos_p.write(line+'\n')
        file_dos_p.write('% 7.6f   % 7.6f \n'%(freq, float(term[1])/norm))
  
        window += float(term[1])/norm
        NN+=1
      
      if(freq<mu):
        LHB += float(term[1])/norm
        #file_dos_lhb.write(line+'\n')
        file_dos_lhb.write('% 7.6f   % 7.6f \n'%(freq, float(term[1])/norm))
      else:
        UHB += float(term[1])/norm
        #file_dos_uhb.write(line+'\n')
        file_dos_uhb.write('% 7.6f   % 7.6f \n'%(freq, float(term[1])/norm))
    
  file_dos_p.write('% 7.6f   % 7.6f  \n'%((w_max_p+mu), 0.0))
  file_dos_lhb.write('% 7.6f   % 7.6f  \n'%(mu, 0.0))
  file_dos_uhb.write('% 7.6f   % 7.6f  \n'%(1000.0,0.0))
        
  #print '\n partial dos sum = % 4.5f' %(
  print 'LHB ~= % 4.5f' %(LHB)
  print 'UHB ~= % 4.5f' %(UHB)
  print 'window_states = % 4.5f' %(window)
  print 
  print 'chemical potential is at frame% 3.2f in dvmc_gnuplot' % N_mu
  print 'chemical potential (precise interpolation) = % 4.5f' %(mu)
  
  
  file_dos_p.close()
  file_dos_lhb.close()
  file_dos_uhb.close()
  
  
  def replaceTemplateFileValues(fileName):
    if os.path.isfile(fileName):
      print('Replacing line "mu = 0.0" by "mu = % 4.3f" in file %s' %(mu,fileName))  
      gnuplotString = open(fileName).read()
      gnuplotString = gnuplotString.replace('mu = 0.0','mu = % 4.5f'%(mu))
      
      fileOut = open('./'+fileName,'w')
      fileOut.write(gnuplotString)
      fileOut.close()
    return;
  
  replaceTemplateFileValues('plot_allAkw.gp')
  replaceTemplateFileValues('plot_Akw.gp')
  replaceTemplateFileValues('plot_dos.gp')
  
  #print ''
  exit()

if __name__ == "__main__":
   main()

