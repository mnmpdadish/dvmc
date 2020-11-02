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

full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)



StdFace = open('StdFace.def').read().replace(' ','')
U=0.
L=W=1
if(StdFace.find('L=')>=0): L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
if(StdFace.find('W=')>=0): W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])
if(StdFace.find('U=')>=0): U = float(re.compile('U=([0-9.]*)').findall(StdFace)[0])

def Xi(ri):
  return ri%W

def Yi(ri):
  return (ri/W)%L

outputDir='output/'
spectrumparaFileName='spectrumpara.def'
verbose_read = 1

sum_rule_max_ok = 1.01
sum_rule_min_ok = 0.96

if (len(sys.argv)>=2):
  spectrumparaFileName=sys.argv[1]
if (len(sys.argv)>=3):
  outputDir=sys.argv[2]+'/'
if (len(sys.argv)>=4):
  verbose_read=int(sys.argv[3])
if (len(sys.argv)>=5):
  sum_rule_min_ok=float(sys.argv[4])
if (len(sys.argv)>=6):
  sum_rule_max_ok=float(sys.argv[5])
if (len(sys.argv)>=7):
  print("example:\n$ vmc_spectrum.py \nor:\n$ vmc_spectrum.py spectrumpara.def\nor:\n$ vmc_spectrum.py spectrumpara.def output/")
  sys.exit()




def dvmc_spectrum(verbose=1):

  sum_rule_max = sum_rule_max_ok
  sum_rule_min = sum_rule_min_ok
#def main():
  zqp_opt_dat = open(outputDir+'zqp_opt.dat').read()
  Omega = float((zqp_opt_dat.split())[0])

  if(os.path.isfile(pythonPathCode+'/libdvmc_speedup.so')):
    lib1 = cdll.LoadLibrary(pythonPathCode+'/libdvmc_speedup.so')
  else:
    print 'Error: the shared object library "libdvmc_speedup.so" was not found.'
    print 'To obtain it, go in the original directory of this python code'
    print 'and compile the code with "make" (choose to link with lapack or mkl).'
    exit(-1)

  w_min_data =-15.0
  w_max_data = 15.0
  eta = 0.2
  Nw = 2000
  exc_choice = [0,1]
  calculateAll = True
  kPath = range(W*L)
  spectrumpara = open(spectrumparaFileName).read()

  for line in  spectrumpara.split('\n'):
    if len(line)>0:
     if line[0]!='#':
      term = line.split()
      if term[0]=='w_min_data' : w_min_data = float(term[1])
      if term[0]=='w_max_data' : w_max_data = float(term[1])
      if term[0]=='eta'   :   eta = float(term[1])
      if term[0]=='Nw'    :    Nw = int(term[1])
#      if (term[0][:]=='kPath' or term[0][:-1]=='kPath'): 
#        if term[1] == 'all':
#          kPath = range(W*L)
#          calculateAll = True
#        elif term[1][0:6]=='range(':
#          kPath = range(int(term[1][6:-1]))
#        else: 
#          kPath = ReadRange(term[1])
      if(term[0][:]=='exc_choice'): 
        if term[1][0:3] == 'all':
          excitation_def = open('excitation.def').read()
          line = excitation_def.split('\n')[1]
          terms = line.split()
          assert(terms[0]=='NExcitation')
          NExcitation = int(terms[1])
          exc_choice = range(NExcitation)
        elif term[1][0:6]=='range(':
          exc_choice = range(int(term[1][6:-1]))
        else: 
          exc_choice = ReadRange(term[1])
        if exc_choice[0] !=0:
          print 'error: first excitation must ALWAYS be 0.'
          exit()
  
  if(verbose):
    print 'excitations chosen:'
    print exc_choice
    
    print '\nOmega=', Omega
    print 'U=', U
    print 'W=', W
    print 'L=', L
  
  n_exc_choice = len(exc_choice)
  Nsite=W*L
  dw = (w_max_data-w_min_data)/(Nw-1)
  w_ = np.array(range(Nw))*dw + w_min_data
  
  lib1 = cdll.LoadLibrary(pythonPathCode+'/libdvmc_speedup.so')
  
  g_ac = np.zeros((Nw),np.cdouble)
  g_ca = np.zeros((Nw),np.cdouble)
  g_ac2 = np.zeros((Nw),np.cdouble)
  g_ca2 = np.zeros((Nw),np.cdouble)
  
  c_darray = np.ctypeslib.ndpointer(np.double)
  c_carray = np.ctypeslib.ndpointer(np.cdouble)
  #lib1.greenFrom_H_and_S.argtypes = [c_int,c_int, c_darray, c_int, c_carray, c_carray, c_double, c_double, c_double, c_carray]
  lib1.greenFrom_e_U_Uinv_S.argtypes = [c_int,c_int, c_darray, c_int, c_carray, c_carray, c_carray, c_double, c_double, c_double, c_carray]

  
  
  S_CA = FFT_selection(outputDir+'nCAm.npy', exc_choice,kPath,verbose)
  S_AC = FFT_selection(outputDir+'nACm.npy', exc_choice,kPath,verbose)
  H_CA = FFT_selection(outputDir+'nCHAm.npy',exc_choice,kPath,verbose)
  H_AC = FFT_selection(outputDir+'nAHCm.npy',exc_choice,kPath,verbose)
  
  np.set_printoptions(precision=9)  
  spectrum_hole = np.zeros([len(kPath),Nw], dtype='float')
  spectrum_elec = np.zeros([len(kPath),Nw], dtype='float')
  
  total_sum = 0.0
  partial_sum = 0.0
  k_label = u'%3s' % 'k#'
  if(verbose):
    print u'\n k#/Nk  --   kx/pi  ky/pi:  sumRule: int dw A(%s,w)  == 1.00000 ' % k_label
    print u' ------------------------------------------------------------'

  for kk in range(len(kPath)):#range(0,2*Nsite):
    k_label = u'%3s' % ('k%d' % kk)
    if(verbose): 
      print u' %2d/%2d ' % (kk+1,len(kPath)),
      sys.stdout.flush()

    #lib1.greenFrom_H_and_S( 1, Nw, w_, n_exc_choice, H_AC[kk], S_AC[kk], Omega, U/2., eta, g_ac)
    #lib1.greenFrom_H_and_S(-1, Nw, w_, n_exc_choice, H_CA[kk], S_CA[kk], Omega, U/2., eta, g_ca)
    
    e_ac,u_ac = la.eig(np.dot(H_AC[kk],la.inv(S_AC[kk])))
    e_ca,u_ca = la.eig(np.dot(H_CA[kk],la.inv(S_CA[kk])))
    
    if(verbose): print '--',
    sys.stdout.flush()
    
    u_ac_m1 = la.inv(u_ac)
    u_ca_m1 = la.inv(u_ca)
    
    us_ac = np.dot(u_ac_m1,S_AC[kk])
    us_ca = np.dot(u_ca_m1,S_CA[kk])
    
    sumRule = 0.0

    # lib1.greenFrom_e_U_Uinv_S compute Charlebois and Imada 2019, arxiv v1 (equation A4) in two parts.
    lib1.greenFrom_e_U_Uinv_S( 1, Nw, w_, n_exc_choice, e_ac, u_ac, us_ac, Omega, U/2., eta, g_ac)
    lib1.greenFrom_e_U_Uinv_S(-1, Nw, w_, n_exc_choice, e_ca, u_ca, us_ca, Omega, U/2., eta, g_ca)
    
    g_tot = g_ac + g_ca # Charlebois and Imada 2019, arxiv v1 (equation A3)
    spectrum_hole[kk,:] = -g_ca[:].imag/(np.pi)
    spectrum_elec[kk,:] = -g_ac[:].imag/(np.pi)
    sumRule = -dw*(g_tot.imag).sum()/(np.pi)

    total_sum += sumRule/float(len(kPath))

    flag = False
    if(sum_rule_min_ok > sumRule):
      flag=True
    if(sum_rule_max_ok < sumRule):
      flag=True
    suffix = ' '  #+ str(sum_rule_min) + ' ' + str(sum_rule_max)
    if(flag): suffix = '   <-- check sum rule' #+ str(sum_rule_min) + ' ' + str(sum_rule_max)
    
    if(verbose): print (u' %2d/%2d  %2d/%2d :           int dw A(%s,w)  = % 5.5f ' \
                     % (Xi(kPath[kk]),W/2,Yi(kPath[kk]),L/2,k_label, sumRule)) +suffix

    if(sum_rule_min > sumRule):
      sum_rule_min = sumRule
    if(sum_rule_max < sumRule):
      sum_rule_max = sumRule
    
  #outputting dos:
  suffix = ''
  totalAkw = spectrum_hole + spectrum_elec
  if calculateAll:
    suffix = '_all'
  
    if(verbose):
      #print u'\n\n                          int dw dk A(k,w)  = % 5.5f ' % (total_sum)
      ##print u'\n\n                                    \u222Bint dw dk A(k,w)  = % 5.5f ' % (total_sum)
      print u'\n\ndos printing'
      print u'\n  correction factor:                int dw dk A(k,w) = % 5.5f ' % (total_sum)
    
    dos = np.zeros([3,Nw], dtype='float')
    dos[0,:] = totalAkw.sum(axis=0)
    dos[1,:] = spectrum_elec.sum(axis=0)
    dos[2,:] = spectrum_hole.sum(axis=0)
    file_dos   = open(outputDir+'dos.dat','w')
    #file_dos_p = open(outputDir+'dos_p'+suffix+'.dat','w')
    for ii in range(Nw):
     file_dos.write('% 7.6f   '  %w_[ii])
     for kk in range(3):
      file_dos.write('% 7.6f '  %(dos[kk,ii]))#/(total_sum)))
     file_dos.write('\n')

  #outputting A(k,w):
  file_green_e = open(outputDir+'Akw_e'+suffix+'.dat','w')
  file_green_h = open(outputDir+'Akw_h'+suffix+'.dat','w')
  file_green   = open(outputDir+'Akw'+suffix+'.dat','w')
  for ii in range(Nw):
   for kk in range(len(kPath)):
    file_green_e.write('% 7.6f '%spectrum_elec[kk,ii])
    file_green_h.write('% 7.6f '%spectrum_hole[kk,ii])
    file_green.write('% 7.6f '  %totalAkw[kk,ii])
   file_green_e.write('\n')
   file_green_h.write('\n')
   file_green.write('\n')


  if((sum_rule_min < sum_rule_min_ok) or (sum_rule_max > sum_rule_max_ok)):
  
    print '\nThe sum rule is not within the range: 0.96 < int dw A(k,w) < 1.00.'
    print 'check range of integration (w_min_data and w_max_data)'
    print 'or increase the number of sampling'
    print 'or exclude this ".bin" file of the analysis.\n'
  else:
  
    print 'OK\n'
  #print ''
  #exit()
  #return sum_rule_min, sum_rule_max
  
  
#######################################################################
######################    SUB ROUTINES    #############################
#######################################################################

def dotdot(a,b,c):
    return np.dot(np.dot(a,b),c)

def FFT_selection(dataFileName,exc_choice,kPath, verbose = 1):
  if(verbose): print 'treatment of '+ dataFileName + '.',
  sys.stdout.flush()
  n_exc_choice = len(exc_choice)
  data_up = np.load(dataFileName)
  
  data_k  = np.fft.fft2(data_up) # fft2 only on the last 2 indices
  dataListOfMatrices = []

  for kk in range(len(kPath)):
    if((kk==len(kPath)-1) and(verbose)): 
        print '.',
        sys.stdout.flush()
    kk1 = kPath[kk]
    kx1 = Xi(kk1)
    ky1 = Yi(kk1)

    tmp1 = data_k[exc_choice,:,kx1,ky1] #slicing is faster than for loops
    tmp2 = tmp1[:,exc_choice]
    dataListOfMatrices.append(tmp2) 
  if(verbose): print ''
  return dataListOfMatrices

def ReadRange(inputStr):
  error=0
  rangeOut = []
  range1 = inputStr.split(',')
  #print range1
  for info in range1:
    element = info.split(':')
    if len(element) ==1:
      rangeOut += [int(element[0])]
    elif len(element) ==2:
      rangeOut += range(int(element[0]),int(element[1])+1)
    else:
      print 'Your definition cannot be interpreted properly.'
      print 'Please input numbers or pairs of numbers'
      print 'sepated by comma (no spaces).'
      print 'Pairs must contains only one ":"' 
      print 'example:'
      print '0:3,5,6,8:11'
      print 
      print 'will be interpreted as the list'
      print '[0,1,2,3,5,6,8,9,10,11]'
      exit()
  return rangeOut


if __name__ == "__main__":
   dvmc_spectrum(verbose_read)



