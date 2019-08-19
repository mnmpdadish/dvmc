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

def main():
  zqp_opt_dat = open('output/zqp_opt.dat').read()
  Omega = float((zqp_opt_dat.split())[0])

  if(os.path.isfile(pythonPathCode+'/c_speedup.so')):
    lib1 = cdll.LoadLibrary(pythonPathCode+'/c_speedup.so')
  else:
    print 'Error: the shared object library "c_speedup.so" was not found.'
    print 'To obtain it, go in the original directory of this python code'
    print 'and do make.'
    exit(-1)

  
  #default parameters
  kPath1 = [0] # default is k=(0,0)
  kPath2 = [0]
  w_min =-8.0
  w_max = 8.0
  eta = 0.2
  Nw = 1000
  exc_choice = [0,1]

  spectrumpara = open('spectrumpara.def').read()

  for line in  spectrumpara.split('\n'):
    if len(line)>0:
     if line[0]!='#':
      term = line.split()
      if term[0]=='w_min' : w_min = float(term[1])
      if term[0]=='w_max' : w_max = float(term[1])
      if term[0]=='eta'   :   eta = float(term[1])
      if term[0]=='Nw'    :    Nw = int(term[1])
      if (term[0][:]=='kPath' or term[0][:-1]=='kPath'): 
        if term[1] == 'all':
          kPath = range(W*L)
        elif term[1][0:6]=='range(':
          #print 'salut'
          kPath = range(int(term[1][6:-1]))
        else: 
          kPath = ReadRange(term[1])
        
        if term[0][-1] == '1': kPath1 = kPath
        elif term[0][-1] == '2': kPath2 = kPath
        else :
          kPath1 = kPath
          kPath2 = kPath
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
  
  print kPath1,kPath2
  print 'excitations chosen:'
  print exc_choice
  
  print '\nOmega=', Omega
  print 'U=', U
  print 'W=', W
  print 'L=', L
  
  n_exc_choice = len(exc_choice)
  Nsite=W*L
  dw = (w_max-w_min)/(Nw-1)
  w_ = np.array(range(Nw))*dw + w_min
  
  #Sif(os.path.isfile(pythonPathCode+'/c_speedup.so')):
  lib1 = cdll.LoadLibrary(pythonPathCode+'/c_speedup.so')
  #c_Nexc = c_int(n_exc_choice)
  #c_eta = c_float(eta)
  #c_ArrayFloatN = c_float * (n_exc_choice)
  #c_green_out = c_ArrayFloatN()
  
  g_ac = np.zeros((Nw),np.cdouble)
  g_ca = np.zeros((Nw),np.cdouble)
  
  #test1 = np.array([[1,2.],[3.,4+4.j]],np.cdouble)
  #lib.foo.argtypes = [np.ctypeslib.ndpointer(np.complex64, ndim=1, flags='C'), c_int]
  c_darray = np.ctypeslib.ndpointer(np.double)
  c_carray = np.ctypeslib.ndpointer(np.cdouble)
  ##lib1.smallTest.argtypes = [c_int, c_darray, c_int, c_carray]
  #lib1.smallTest(c_int(Nw), w_.ctypes, c_int(2), test1.ctypes)
  #lib1.smallTest(Nw, w_, 2, test1)
  #print test1
  #exit()
  lib1.greenFrom_H_and_S.argtypes = [c_int,c_int, c_darray, c_int, c_carray, c_carray, c_double, c_double, c_double, c_carray]

#    lib = ctypes.cdll[plat_info['lib']]
#    fct = lib.lowPassFilterGen
#    fct.restype = ctypes.c_double
#    fct.argtypes = (ndpointer(dtype=ctypes.c_double, shape=(len(x),)), 
#                    ctypes.c_uint64, 
#                    ctypes.c_double, 
#                    ctypes.c_double, 
#                    ctypes.c_double)
  
  
  S_CA = FFT_selection('output/nCAm.npy', exc_choice,kPath1)
  S_AC = FFT_selection('output/nACm.npy', exc_choice,kPath1)
  H_CA = FFT_selection('output/nCHAm.npy',exc_choice,kPath1)
  H_AC = FFT_selection('output/nAHCm.npy',exc_choice,kPath1)
  np.save('tmpS_CA',S_CA)
  np.save('tmpH_CA',H_CA)
  np.save('tmpS_AC',S_AC)
  np.save('tmpH_AC',H_AC)
  
  #S_CA=np.load('tmpS_CA.npy')
  #H_CA=np.load('tmpH_CA.npy')
  #S_AC=np.load('tmpS_AC.npy')
  #H_AC=np.load('tmpH_AC.npy')
  
  #print S_CA
  #print S_CA2
  #print S_CA == S_CA2
  #exit()
  
  #if(kPath2!=kPath1):
  
  
  np.set_printoptions(precision=9)  
  spectrum_hole = np.zeros([len(kPath1),Nw], dtype='float')
  spectrum_elec = np.zeros([len(kPath1),Nw], dtype='float')
  
  k_label = u'%3s' % 'k#'
  print u'\n k#/Nk  --   kx/\u03C0   ky/\u03C0 :  sumRule: \u222Bdw A(%s,w)  \u225F  1.00000 ' % k_label
  print u' ------------------------------------------------------------'
  for kk in range(len(kPath1)):#range(0,2*Nsite):
    k_label = u'%3s' % ('k%d' % kk)
    print u' %2d/%2d  --' % (kk+1,len(kPath1)),
    sys.stdout.flush()
    
    #H = np.ascontiguousarray(H_AC[kk], dtype=np.cdouble)
    #S = np.ascontiguousarray(S_AC[kk], dtype=np.cdouble)
    lib1.greenFrom_H_and_S(1, Nw, w_, n_exc_choice, H_AC[kk], S_AC[kk], Omega, U/2., eta, g_ac)
    
    #H = np.ascontiguousarray(H_CA[kk], dtype=np.cdouble)
    #S = np.ascontiguousarray(S_CA[kk], dtype=np.cdouble)
    lib1.greenFrom_H_and_S(-1, Nw, w_, n_exc_choice, H_CA[kk], S_CA[kk], Omega, U/2., eta, g_ca)
    
    
    sumRule = 0.0

    g_tot = g_ac + g_ca
    spectrum_hole[kk,:] = -g_ca[:].imag/(np.pi)
    spectrum_elec[kk,:] = -g_ac[:].imag/(np.pi)
    sumRule = -dw*(g_tot.imag).sum()/(np.pi)
      
      

      
    print u' %2d/%2d  %2d/%2d :           \u222Bdw A(%s,w)  = % 5.5f ' \
             % (Xi(kPath1[kk]),W/2,Yi(kPath1[kk]),L/2,k_label, sumRule)
    
  file_green_e = open('output/Akw_e.dat','w')
  file_green_h = open('output/Akw_h.dat','w')
  file_green   = open('output/Akw.dat','w')
  for ii in range(Nw):
   for kk in range(len(kPath1)):
    file_green_e.write('% 7.6f '%spectrum_elec[kk,ii])
    file_green_h.write('% 7.6f '%spectrum_hole[kk,ii])
    file_green.write('% 7.6f '%(spectrum_hole[kk,ii]+spectrum_elec[kk,ii]))
   file_green_e.write('\n')
   file_green_h.write('\n')
   file_green.write('\n')

  #print ''
  exit()
  
  #X,Y = np.meshgrid(w_,np.array(range(0,1+len(kPath1))))
  #Z = np.zeros([len(kPath1),Nw], dtype='float')
  
  #for ii in range(len(kPath1)):
  #  for jj in range(Nw):
  #    Z[ii,jj] = spectrum_hole[ii,jj] + spectrum_elec[ii,jj]
  
  #fName = 'figure.pdf'
  
  #fig= plt.figure(figsize=(5,6))
  #ax= fig.add_axes([0.1,0.1,0.8,0.8])
  
  #Z2 = np.log(np.abs(Z))
  #plt.pcolormesh(np.transpose(Y),np.transpose(X),np.transpose(Z2),cmap=cm.afmhot_r)
  ##plt.pcolormesh(np.transpose(Z),cmap=cm.afmhot_r)
  #plt.xticks(np.array([1./33.,32./33.])*len(kPath1), (r'$0$', r'$\pi$'))
  #plt.yticks([-8,-4,0,4,8])
  #plt.clim(-4.0,Z2.max()*1.0);

  #y0 = np.array(range(0,1+len(kPath1)))
  #x0 = 0*y0
  #ax.plot(y0, x0,'-', color = 'black',lw=1.4)
  #ax.plot(y0, x0,'--', color = 'white',lw=1.5)
  #plt.ylim([w_min,w_max])
  #plt.savefig(fName)
  
  #plt.show()
  # paste -d" " <( awk '{printf "% e \n", $7 }' output/zvo_physAC_001.dat) <(awk '{printf "% e \n", $5}' output/zvo_phys_nAHCm_001.dat) > tmp1; head tmp1
  
#######################################################################
######################    SUB ROUTINES    #############################
#######################################################################

def dotdot(a,b,c):
    return np.dot(np.dot(a,b),c)

def FFT_selection(dataFileName,exc_choice,kPath):
  print 'treatement of '+ dataFileName + '.',
  sys.stdout.flush()
  n_exc_choice = len(exc_choice)
  data_up = np.load(dataFileName)
  
  data_k  = np.fft.fft2(data_up) # fft2 only on the last 2 indices
  dataListOfMatrices = []

  #for kk in range(0,len(kPath)):
  #  dataListOfMatrices.append(np.zeros([n_exc_choice,n_exc_choice], dtype='cfloat'))

  for kk in range(len(kPath)):
    if((kk)%(len(kPath)/5)==0): 
      print '.',
      sys.stdout.flush()
    kk1 = kPath[kk]
    #kk2 = kPath2[kk]
    #print kk, kk1
    kx1 = Xi(kk1)
    ky1 = Yi(kk1)
    #kx2 = Xi(kk2)
    #ky2 = Yi(kk2)

    #for nn in range(0,n_exc_choice):
    #  for mm in range(0,n_exc_choice):
    tmp1 = data_k[exc_choice,:,kx1,ky1] #slicing is faster than for loops
    tmp2 = tmp1[:,exc_choice]
    #dataListOfMatrices[kk] = 0.5*(data_k[exc_choice,exc_choice,kx1,ky1]  +  data_k[exc_choice,exc_choice,kx2,ky2])  #averaging 2 paths
    dataListOfMatrices.append(tmp2) 
  print ''
  return dataListOfMatrices

def ReadRange(inputStr):
  rangeOut = []
  range1 = inputStr.split(',')
  print range1
  for info in range1:
    element = info.split(':')
    if len(element) ==1:
      rangeOut += [int(element[0])]
    elif len(element) ==2:
      rangeOut += range(int(element[0]),int(element[1])+1)
    else:
      print 'Your definition cannot be interpreted properly.'
      print 'Please input numbers or pairs of numbers'
      print 'sepated by comma. Pairs must contains only one ":"' 
      print 'example:'
      print '0:3,5,6,8:11'
      print 
      print 'will be interpreted as'
      print '[0,1,2,3,5,6,8,9,10,11]'
      exit()
  return rangeOut


if __name__ == "__main__":
   main()
   
   
# mergeOutput.py output/zvo_nCHAm_nAHCm_0*
# print_spectrum.py
### generate_template_gnuplot.py ###
# gnuplot plot_singleAkw.gp
# gnuplot plot_allAkw.gp

