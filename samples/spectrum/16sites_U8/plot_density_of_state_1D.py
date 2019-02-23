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
#   misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.

import numpy as np
from scipy import linalg as LA
from copy  import deepcopy
import sys, os, re
from scipy.linalg import eigh,eig,eigvals,eigvalsh

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D


def main():

  excitations =[[0, -1],
                [0, 0],
                [0, 1],
                [0, 2],
                [1, 0],
                [1, 1],
                [1, 2],
                [2, 0],
                [2, 1],
                [2, 2],
                [3, 0],
                [3, 1],
                [3, 2],
                [4, 0],
                [4, 1],
                [4, 2]]
  
  #excitations_choice = [3,4,6,7,9,10,12,13]
  #excitations_choice = [0,4,7,10]
  excitations_choice = [0,1,2,3]
  
  directory = "./" 
  list_to_print = range(9)

  w_min =-7.0
  w_max = 7.0
  eta = 0.2
  Nw = 1000
  
  name_nACm  = directory+"output/zvo_phys_nACm_001.dat"
  name_nCAm  = directory+"output/zvo_phys_nCAm_001.dat"
  name_nAHCm = directory+"output/zvo_phys_nAHCm_001.dat"
  name_nCHAm = directory+"output/zvo_phys_nCHAm_001.dat"
  
  zqp_opt_dat = ReadFile(directory+"output/zqp_opt.dat")[0]
  Omega = float((zqp_opt_dat.split())[0])
  
  
  StdFace = ReadFile(directory+"StdFace.def")
  
  U=0.0
  W=1
  L=1
  
  def Xi(ri):
    return ri%W

  def Yi(ri):
    return (ri/W)%L

  def Ri(xi, yi):
    return yi*W + xi;
  
  for line in StdFace:
    if line[0]=='U': U = float(line.split('=')[1])
    if line[0]=='W': W = int(line.split('=')[1])
    if line[0]=='L': L = int(line.split('=')[1])
  
  mu = U/2
        
  
  print 'Omega=', Omega, '   U=', U, '   mu=', mu
  print 'W=', W, 'L=', L
  
  #exit()
  #print excitations
  #print excitations[excitations_choice[0]][0]
  #print excitations[excitations_choice[1]][0]
  
  n_exc = len (excitations_choice)
  
  Nsite=W*L
  sublatt = 1

  
  
  
  def ReadData(fileIn):
    matrix_list=[]
    for ii in range(0,n_exc*n_exc):
      matrix_list.append(np.zeros([W,L], dtype='f'))
    print len(matrix_list)
    for line in fileIn:
      vals = line.split()
      if vals[0] != '#':
        #print vals
        ri=int(vals[0])
        s= int(vals[1])
        rj=int(vals[2])
        drx = (Xi(rj)-Xi(ri))%W
        dry = (Yi(rj)-Yi(ri))%L
        fact = 1./Nsite

        #if ri!=0: break
        if(s==0):
          #print drx, dry, vals[4]
          for nn in range(0,n_exc):
            for mm in range(0,n_exc):
              
              #tmp = fact*float(vals[4+excitations_choice[nn]+excitations_choice[mm]*n_exc])
              #print tmp
              #matrix_list[nn+mm*n_exc][drx,dry]=0.0  
              matrix_list[nn+mm*n_exc][drx,dry]  += \
                             fact*float(vals[4+excitations_choice[nn]+excitations_choice[mm]*n_exc])
#              matrix_list[mm+nn*n_exc][drx,dry]  += \
#                            0.5*fact*float(vals[4+excitations_choice[nn]+excitations_choice[mm]*n_exc])
                            
    return deepcopy(matrix_list)

  nACm_up =ReadData(ReadFile(name_nACm))
  nCAm_up =ReadData(ReadFile(name_nCAm))
  nAHCm_up=ReadData(ReadFile(name_nAHCm))
  nCHAm_up=ReadData(ReadFile(name_nCHAm))
  print '\n\ntest\n'
  for ii in range(len(nACm_up)):
    print nACm_up[ii]
    print
    print nCAm_up[ii]
    print
    print nAHCm_up[ii]
    print
    print nCHAm_up[ii]
    print
    print
  #exit()

  nACm_up_k =[None]*n_exc*n_exc
  nCAm_up_k =[None]*n_exc*n_exc
  nAHCm_up_k=[None]*n_exc*n_exc
  nCHAm_up_k=[None]*n_exc*n_exc
  
  np.set_printoptions(precision=9)
  
  for ii in range(0,n_exc*n_exc):
    #printmatrix(MACM_up[ii])
    #print
    
    nACm_up_k[ii] = np.fft.fft2(nACm_up[ii])
    nCAm_up_k[ii] = np.fft.fft2(nCAm_up[ii])
    nAHCm_up_k[ii]= np.fft.fft2(nAHCm_up[ii])
    nCHAm_up_k[ii]= np.fft.fft2(nCHAm_up[ii])
    #nCHAm_up_k[ii]= np.fft.fft2(nCHAm_up[ii]-0.8)
    #nCHAm_up_k[ii]= np.fft.fft(nCHAm_up[ii][0])

  print '\n\ntest_k\n'
  for ii in range(len(nACm_up_k)):
    print nACm_up_k[ii]
    print
    print nCAm_up_k[ii]
    print
    print nAHCm_up_k[ii]
    print
    print nCHAm_up_k[ii]
    print
    print
  #exit()

  S_AC = []
  S_CA = []
  H_AC = []
  H_CA = []
  
  eigenval_AC = []
  eigenval_CA = []
  
  for kk in range(0,2*Nsite):
    print kk
    S_AC.append(np.zeros([n_exc,n_exc], dtype='cfloat'))
    S_CA.append(np.zeros([n_exc,n_exc], dtype='cfloat'))
    H_AC.append(np.zeros([n_exc,n_exc], dtype='cfloat'))
    H_CA.append(np.zeros([n_exc,n_exc], dtype='cfloat'))
  #print n_exc
  #exit()

  spectrum = np.zeros([len(list_to_print),Nw], dtype='float')
  
  for kk in range(0,2*Nsite):
    kx = Xi(kk)
    ky = Yi(kk)
    print "\n\n##################"
    print kk, kx, ky
    
    for nn in range(0,n_exc):
      for mm in range(0,n_exc): 
        S_AC[kk][nn,mm] =  nACm_up_k[nn+mm*n_exc][kx,ky] 
        H_AC[kk][nn,mm] = nAHCm_up_k[nn+mm*n_exc][kx,ky] 
        S_CA[kk][nn,mm] =  nCAm_up_k[nn+mm*n_exc][kx,ky] 
        H_CA[kk][nn,mm] = nCHAm_up_k[nn+mm*n_exc][kx,ky] 
    
    print H_CA[kk]
    print H_AC[kk]
    
    print S_CA[kk]
    print S_AC[kk]
    
    #print2matrices(H_CA[kk], H_AC[kk])
    #print2matrices(S_CA[kk], S_AC[kk])
      
    hermitizeMatrix(S_AC[kk])
    hermitizeMatrix(H_AC[kk])
    hermitizeMatrix(S_CA[kk])
    hermitizeMatrix(H_CA[kk])

  #exit()    
  
  for kk in range(0,2*Nsite):
    
    dw = (w_max-w_min)/(Nw-1)
    w_ = np.array(range(Nw))*dw+w_min
    
    if kk in list_to_print:
    
      small_value = 0.0000#1
      
      for ii in range(S_CA[kk].shape[0]):      
        S_CA[kk][ii,ii] += small_value
        S_AC[kk][ii,ii] += small_value
      
      print
      print
      print2matrices(S_AC[kk],S_CA[kk])
      print
      print2matrices(H_AC[kk],H_CA[kk])
      print
      
      e_ca1= eigvals(H_CA[kk], S_CA[kk])
      e_ac1= eigvals(H_AC[kk], S_AC[kk])
      for nn in range(n_exc):
        print pretty_c(e_ac1[nn]-Omega),
      for nn in range(n_exc):
        print pretty_c(e_ca1[nn]-Omega),
      print 
      print
      print kk, list_to_print.index(kk)
      indices = [i for i, x in enumerate(list_to_print) if x == kk]
      
      for ii in range(Nw):
        w = w_[ii]
        omeg = Omega
        
        G_AC = np.dot(S_AC[kk],np.dot(LA.inv((w + 1.0j*eta + omeg + mu)*S_AC[kk] - H_AC[kk]), S_AC[kk]))
        G_CA = np.dot(S_CA[kk],np.dot(LA.inv((w + 1.0j*eta - omeg + mu)*S_CA[kk] + H_CA[kk]), S_CA[kk]))        
        
        #G_AC = S_AC0[kk]*S_AC0[kk]/((w + 1.0j*eta + Omega + mu)*S_AC0[kk] - H_AC0[kk])
        #G_CA = S_CA0[kk]*S_CA0[kk]/((w + 1.0j*eta - Omega + mu)*S_CA0[kk] + H_CA0[kk])
        
        #G_AC = np.dot(S_AC[kk],LA.inv((w + 1.0j*eta + Omega + mu)*S_AC[kk] - H_AC[kk]))
        #G_CA = np.dot(S_CA[kk],LA.inv((w + 1.0j*eta - Omega + mu)*S_CA[kk] + H_CA[kk]))
        #G_AC = np.dot(LA.inv((w + 1.0j*eta + Omega)*S_AC[kk] - H_AC[kk]), S_AC[kk])
        #G_CA = np.dot(LA.inv((w + 1.0j*eta - Omega)*S_CA[kk] + H_CA[kk]), S_CA[kk])
        #G_AC = np.dot(LA.inv((w + 1.0j*eta + Omega + mu)*S_AC[kk] - H_AC[kk]), S_AC[kk])
        #G_CA = np.dot(LA.inv((w + 1.0j*eta - Omega + mu)*S_CA[kk] + H_CA[kk]), S_CA[kk])
         
        #S_AC_inv = LA.inv(S_AC[kk])
        #S_CA_inv = LA.inv(S_CA[kk])
        
        #G_AC = np.dot(LA.inv((w + 1.0j*eta + Omega + mu)*S_AC_inv - np.dot(S_AC_inv,np.dot(H_AC[kk], S_AC_inv) )), S_AC_inv)
        #G_CA = np.dot(LA.inv((w + 1.0j*eta - Omega + mu)*S_CA_inv + np.dot(S_CA_inv,np.dot(H_CA[kk], S_CA_inv) )), S_CA_inv)
        
        tmp = -G_CA[0,0] -G_AC[0,0]    
        
        for ind in indices:
          spectrum[ind,ii] = tmp.imag/(np.pi)
        
  print
  X,Y = np.meshgrid(w_,np.array(range(0,1+len(list_to_print))))
  
  Z = np.zeros([len(list_to_print),Nw], dtype='float')
  Z_sum = np.zeros(len(list_to_print), dtype='float')
  for ii in range(len(list_to_print)):
    for jj in range(Nw):
      Z_sum[ii] += spectrum[ii,jj]*dw
    print Z_sum[ii]
  
  for ii in range(len(list_to_print)):
    for jj in range(Nw):
      Z[ii,jj] = spectrum[ii,jj] #/Z_sum[ii]   #/max(spectrum[ii,:])
    
  fig, ax = plt.subplots()
  
  #print X
  #print Y
  #print Z
  
  Z2 = np.log(np.abs(Z))
  plt.pcolormesh(np.transpose(Y),np.transpose(X),np.transpose(Z2),cmap=cm.afmhot_r)
  #plt.pcolormesh(np.transpose(Z),cmap=cm.afmhot_r)
  plt.xticks([0,len(list_to_print)], (r'$0$', '$\pi$'))
  plt.yticks([-5,0,5])
  #print(Z2.min(),Z2.max())
  plt.clim(-3,Z2.max());
  #plt.clim(Z2.min(),Z2.max());
  #plt.clim(0,1.0);
  #plt.clim(0.0,0.5*Z2.max());
  #cbar = plt.colorbar()
  #cutoff = round(np.exp(Z2.max()))
  #cbar.set_ticks([])  # vertically oriented colorbar
  #cbar.ax.set_yticklabels([])  # vertically oriented colorbar
  plt.savefig("figure"+".pdf")
  plt.show()
  #ax.plot_surface(X, Y, spectrum, rstride=1, cstride=1, cmap=cm.coolwarm)
  #plt.show()



# paste -d" " <( awk '{printf "% e \n", $7 }' output/zvo_physAC_001.dat) <(awk '{printf "% e \n", $5}' output/zvo_phys_nAHCm_001.dat) > tmp1; head tmp1
  
###########################################################################################################
###########################################################################################################

















def ReadFile(fileName):
  fileExist = os.path.isfile(fileName)
  if not fileExist:
    print '\nerror: file '+fileName+' does not exist.'   
    print 'terminated.'
    exit()
   
  f = open(fileName, 'r')
  file1 = f.readlines()
  f.close
  
  file2 = []
  for lines in file1:
    file2.append(lines.strip())
  return file2

def hermitizeMatrix(matrix):
  for ii in range(matrix.shape[0]):
    for jj in range(ii,matrix.shape[1]):
      tmp = 0.5*(matrix[ii,jj] + np.conj(matrix[jj,ii]))
      matrix[ii,jj]=tmp
      matrix[jj,ii]=np.conj(tmp)

def makePositive(matrix):
  for ii in range(matrix.shape[0]):
    for jj in range(matrix.shape[1]):
      matrix[ii,jj]=abs(matrix[ii,jj])

def pretty(float1):
  if abs(float1)<1e-8:
    return '  . '#%(float1)   
  return '% 4.3f'%(float1)   
  
def choose_e_or_f(value):
  if abs(value)<1e-2:
    return '% 8.1e'%(value)
  else:
    return '% 8.4f'%(value)
    
def pretty_c(complex1):
  string = ' '
  if abs(complex1.real)<1e-8:
    string += '  .     '
  else:
    string += choose_e_or_f(complex1.real)
  string += ' '
  if abs(complex1.imag)<1e-8:
    string += '    '
  else:
    string += choose_e_or_f(complex1.imag)+'i'
  return string

def print2matrices(matrix1,matrix2):
  Nlines = matrix1.shape[0]
  Ncol1 =  matrix1.shape[1]
  Ncol2 =  matrix2.shape[1]
  assert(Nlines == matrix2.shape[0])
  for ii in range(Nlines):
    for jj in range(Ncol1):
      print pretty_c(matrix1[ii,jj]),
    print '    ',
    for jj in range(Ncol2):
      print pretty_c(matrix2[ii,jj]),
    print

def printmatrix(matrix1):
  Nlines = matrix1.shape[0]
  Ncol1 =  matrix1.shape[1]
  for ii in range(Nlines):
    for jj in range(Ncol1):
      print pretty_c(matrix1[ii,jj]),
    print


if __name__ == "__main__":
   main()
