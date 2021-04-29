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
import os

#some functions and classes:

def C(integerState, flavor):
      # 0: empty state
      # -1: void state
      # greater than zero: any state
   tmp = (1<<flavor)
   isNull = integerState | tmp == integerState
   if isNull:
      return -1  # -1 is the void
   else:
      return integerState + tmp

def A(integerState, flavor):
   tmp = (1<<flavor)
   isNull = integerState & tmp == 0
   if isNull:
      return -1  # -1 is the void
   else:
      return integerState - tmp

def CA(integerState, flavor_i,flavor_j):
   tmp0 = self.A(integerState,flavor_j)
   if tmp0 != -1:
      return Create_fast(tmp0,flavor_i)
   else: return -1  # -1 is the void



class qState:
   
   def __init__(self,Nsites,coeff): #nspin=2 always
      self.N = (1<<(Nsites*2))
      self.Nsites = Nsites
      self.ket = np.zeros((self.N), dtype=complex) #[0.00]*self.N
      self.ket[0] = coeff  #empty state with coefficient 
   
   def squareNorm(self):
      sum1=0.0
      for ii in range(self.N):
         sum1+=self.ket[ii]*np.conjugate(self.ket[ii])
      return sum1

   def scalarProd(self,qState2):
      assert(qState2.N==self.N)
      assert(qState2.Nsites==self.Nsites)
      sum1=0.0
      for ii in range(self.N):
         sum1+= np.conj(qState2.ket[ii]) * self.ket[ii]
      return sum1
   
   def equal(self,qState2,coeff):
      self.N = qState2.N
      self.Nsites = qState2.Nsites
      for ii in range(self.N):
         self.ket[ii] = coeff*qState2.ket[ii]
      #print(self.string2())
   
   def add(self,qState2,coeff):
      assert(qState2.N==self.N)
      assert(qState2.Nsites==self.Nsites)
      for ii in range(self.N):
         self.ket[ii] += coeff*qState2.ket[ii]
      #print(self.string())
   
   def create(self,flavor,coeff):
      #first: erase doubly occupied for flavor
      for ii in range(self.N):
         if C(ii,flavor) == -1:
            self.ket[ii] = 0.0
      for ii in range(self.N):
         newState = C(ii,flavor)
         numberOfBitsToPermute = bin(ii >> (flavor+1)).count('1')
         tmpCoeff = 1.0 
         if (numberOfBitsToPermute % 2 == 1): 
            tmpCoeff = -1.0   # anticommutation if commute with odd number of operator
         if (newState != -1) and (newState < self.N):
            self.ket[newState] = tmpCoeff*coeff*self.ket[ii]
            self.ket[ii]=0.0
      #print(self.string2())
            
   def annihi(self,flavor,coeff):
      #first: erase the voids
      for ii in range(self.N):
         if A(ii,flavor) == -1:
            self.ket[ii] = 0.0
      for ii in range(self.N):
         newState = A(ii,flavor)
         numberOfBitsToPermute = bin(ii >> (flavor+1)).count('1')
         tmpCoeff = 1.0 
         if (numberOfBitsToPermute % 2 == 1): 
            tmpCoeff = -1.0   # anticommutation if commute with odd number of operator
         if (newState != -1) and (newState < self.N):
            self.ket[newState] = tmpCoeff*coeff*self.ket[ii]
            self.ket[ii]=0.0
      #print(self.string2())
      

   def string(self):
      s1 = ""
      for ii in range(self.N):
         s1 += bin(ii)[2:].zfill(self.Nsites*2) + " "
      s1 += '\n'
      for ii in range(self.N):
         if(abs(self.ket[ii].real)>0.000001):
            s1 += "% 2.1f "%(self.ket[ii].real)
         else:
            s1 += "  .  "
      s1+="\n"
      for ii in range(self.N):
         if(abs(self.ket[ii].imag)>0.000001):
            s1 += "% 2.1f "%(self.ket[ii].imag)
         else:
            s1 += "  .  "
      return s1
   
   def string2(self):
      s1 = ""
      for ii in range(self.N):
         if(abs(self.ket[ii]) >0.0001):
            s1 = bin(ii)[2:].zfill(self.Nsites*2)+", % 2.1f % 2.1fi "%(self.ket[ii].real,self.ket[ii].imag)
      return s1
   
   def __str__(self):
      return self.string()
   
   def __add__(self,qState2):
      tmp_phi = qState(self.Nsites,0.0)  #recipient
      tmp_phi.equal(self,1.0)
      tmp_phi.add(qState2,1.0)
      return tmp_phi
      
   def __rmul__(self,some_float):
      assert((type(some_float) is float) or (type(some_float) is int)  or (type(some_float) is complex) )
      
      tmp_phi = qState(self.Nsites,0.0)  #recipient
      tmp_phi.equal(self,1.0)
      
      for ii in range(self.N):
         tmp_phi.ket[ii] = some_float*tmp_phi.ket[ii]
      
      return tmp_phi


      
class c_:

  def __init__(self,flavor): 
    self.flavor = flavor 
    self.coeff = 1.0 
    
  def __mul__(self,something_else):
    if(type(something_else) is qState):
      tmp_state = qState(something_else.Nsites,1.0)
      #tmp_state = qState(self.n_sites,0.0)
      tmp_state.equal(something_else,1.0)
      tmp_state.create(self.flavor,self.coeff)
      return tmp_state
    if(type(something_else) is a_):
      tmp_state = qState(something_else.Nsites,1.0)
      #tmp_state = qState(self.n_sites,0.0)
      tmp_state.equal(something_else,1.0)
      tmp_state.create(self.flavor,self.coeff)
      return tmp_state
    if((type(something_else) is float) or (type(something_else) is int)  or (type(something_else) is complex) ):
      tmp_c = c_(self.flavor)
      tmp_c.coeff = something_else*self.coeff
      return tmp_c
    else:
      print('?')
      exit(1)
    
    
class a_:

  def __init__(self,flavor): #nspin=2 always
    self.flavor = flavor 
    self.coeff = 1.0 
    
  def __mul__(self,something_else):
    if(type(something_else) is qState):
      tmp_state = qState(something_else.Nsites,1.0)
      #tmp_state = qState(self.n_sites,0.0)
      tmp_state.equal(something_else,1.0)
      tmp_state.annihi(self.flavor,self.coeff)
      return tmp_state
    if((type(something_else) is float) or (type(something_else) is int)  or (type(something_else) is complex) ):
      tmp_c = a_(self.flavor)
      tmp_c.coeff = something_else*self.coeff
      return tmp_c
    
      
class n_:

  def __init__(self,flavor): 
    self.flavor = flavor 
    self.coeff = 1.0 
    
  def __mul__(self,something_else):
    return self.coeff * (c_(self.flavor) * (a_(self.flavor) * something_else))
         
      
      
########################
##### main program #####
########################

#nSites = 2
#nSites = 3
nSites = 2
Read = 1
Write = 0
Symmetrize = 0

t = -1.0
U = 8.0

#U=0
#f_[1][1] = 1.632588035
#f_[1][0] = 3.927111666
#f_[0][1] = 3.927111666
#f_[0][0] = 4.0

H_T = np.zeros((nSites,nSites), dtype=float)
H_U = np.zeros((nSites), dtype=float)

for ii in range(nSites):
  for jj in range(nSites):
    if((ii-jj)%nSites==1 or (ii-jj)%nSites==3):
      if(nSites==2):
        H_T[ii][jj] += 2*t  
      else:
        H_T[ii][jj] += t  

for ii in range(nSites):
  H_U[ii] += U
  #H_T[ii][ii] += 2*t


for ii in range(nSites):
  for jj in range(nSites):
    print(ii,jj,H_T[ii][jj])

#exit()

f_ = np.zeros((nSites,nSites), dtype=complex)
#f_[1][1] = 1.0
#f_[1][0] = 1.0
#f_[0][1] = 1.0
#f_[0][0] = 1.0



if(Read):

  def ReadFile(fileName):
    fileExist = os.path.isfile(fileName)
    if not fileExist:
      print('\nerror: file '+fileName+' does not exist.')   
      print('terminated.')
      exit()
     
    f = open(fileName, 'r')
    file1 = f.readlines()
    f.close
    
    file2 = []
    for lines in file1:
      file2.append(lines.strip())
    return file2

  print('\nreading file ./output/zqp_opt.dat')
  data = ReadFile("./output/zqp_opt.dat")
  data_list = data[0].split()
  print(len(data_list))
  
  for ii in range(nSites):
    for jj in range(nSites):
      data_ii_re = data_list[len(data_list)-3-3*(ii+jj*nSites)]
      data_ii_im = data_list[len(data_list)-2-3*(ii+jj*nSites)]
      f_[ii][jj] += float(data_ii_re) + 1.j * float(data_ii_im)
      #f_[ii][jj] += 0.5*float(data_ii)
      #f_[jj][ii] += 0.5*float(data_ii)
      print(ii,jj,f_[ii][jj])
  
  print('')
  for ii in range(nSites):
    for jj in range(nSites):
      print(f_[ii][jj])
  #exit()


if(Write):
  # ### zqp_opt.dat ###
  print('\nprinting file ./zqp_opt.dat')
  ofile = open("./zqp_opt.dat",'w')
  ofile.write("{0:.10f} {1:.5f} {2:.5f}".format(0.0,0.0,0.0))
  for idx in range(4): #hmm
    ofile.write("{0:.5f} {1:.5f} {2:.5f} ".format(0.0,0.0,0.0))
  for ii in range(nSites):
    for jj in range(nSites):
      print(ii, jj, f_[ii][jj].real, f_[ii][jj].imag)
      ofile.write("{0:.10f} {1:.10f} {2:.5f} ".format(f_[ii][jj].real,f_[ii][jj].imag,0.0))
  ofile.close()

ftmp_ = np.zeros((nSites,nSites), dtype=complex)
if(Symmetrize):
  for ii in range(nSites):
    for jj in range(nSites):
      diff = (jj-ii)%nSites
      ftmp_[0][diff] += (1./nSites)*f_[ii][jj]
  for ii in range(nSites):
    for jj in range(nSites):
      diff = (jj-ii)%nSites
      ftmp_[ii][jj] = ftmp_[0][diff]
  
  f_ = ftmp_
  print('')
  for ii in range(nSites):
    for jj in range(nSites):
      print(f_[ii][jj])

#exit()

phi_0 = qState(nSites,1.0)
phi_pair = qState(nSites,0.0) #recipient 

for ii in range(0,nSites):
  for jj in range(0,nSites):
    #print("ii,jj= ",ii,jj)
    phi_pair += f_[ii][jj] * (c_(nSites+ii) * (c_(jj) * phi_0))
print(phi_pair)

#exit(0)
#print(n_(3)*phi_pair)

def flavor(i,nSites):
  site = i%nSites
  spin = i//nSites
  return '%1d %1d '%(site,spin)   

def pretty(float1):
  if abs(float1)<1e-8:
    return '  . '#%(float1)   
  return '% 4.3e '%(float1)   



for aa in range(0,nSites):
 for ss in range(0,2):
  for bb in range(0,nSites):
      s = ' %d %d %d %d '%(aa,ss,bb,ss)
      
      psi_2 = (c_(aa+ss*nSites)*(a_(bb+ss*nSites)*(phi_pair)))
      val_2 = psi_2.scalarProd(phi_pair)/phi_pair.squareNorm()      
      s+= ' % 4.4e  '%(val_2.real)

      psi_2 = (a_(aa+ss*nSites)*(c_(bb+ss*nSites)*(phi_pair)))
      val_2 = psi_2.scalarProd(phi_pair)/phi_pair.squareNorm()      
      s+= ' % 4.4e  '%(val_2.real)
      print(s)

print('')


# cat output/zvo_phys_nACm_001.dat; printf "\n\n"; cat output/zvo_phys_nCAm_001.dat; printf "\n\n"; cat output/zvo_phys_nAHCm_001.dat; printf "\n\n"; cat output/zvo_phys_nCHAm_001.dat; printf "\n"
def get_excitations(ii,tt):
  #psi1 = 1.0*(phi_pair)
  #psi2 = 1.0*(n_(ii+tt*nSites)*(phi_pair)) 
  psi0 = 1.0*(phi_pair)
  psi1 = 1.0*(n_((ii+0)%nSites+(1-tt)*nSites)*(phi_pair)) 
  
  #psi4 = 1.0*(n_((ii+1)%nSites+(  tt)*nSites)*(phi_pair)) 
  #psi5 = 1.0*(n_((ii+1)%nSites+(1-tt)*nSites)*(phi_pair)) 
  #psi6 = 1.0*(n_((ii+1)%nSites+(1-tt)*nSites)*(n_((ii+1)%nSites+(tt)*nSites)*(phi_pair))) 

  #psi14 = 1.0*(n_((ii+1)%nSites+(  tt)*nSites)*(n_((ii+0)%nSites+(1-tt)*nSites)*(phi_pair))) 
  #psi15 = 1.0*(n_((ii+1)%nSites+(1-tt)*nSites)*(n_((ii+0)%nSites+(1-tt)*nSites)*(phi_pair)))
  
  #psi16 = 1.0*(n_((ii+1)%nSites+(1-tt)*nSites)*(phi_pair))
  psi3 = 1.0*(n_((ii+1)%nSites+(1-tt)*nSites)*(n_((ii+1)%nSites+(1-tt)*nSites)*(phi_pair)))
  psi4 = 1.0*(n_((ii+1)%nSites+(1-tt)*nSites)*(n_((ii+1)%nSites+(tt)*nSites)*(phi_pair)))
  
  
  return [psi0,psi1,psi3, psi4]
  
nExc = len(get_excitations(0,0))

fileo = open('output/green_comp_ave.dat', 'w')
fileo.write('#<cHa>\n')
#print('\n< cHa >')



Ntot = nSites*nExc
S_AC = np.zeros([Ntot,Ntot])
S_CA = np.zeros([Ntot,Ntot])
H_AC = np.zeros([Ntot,Ntot])
H_CA = np.zeros([Ntot,Ntot])

for aa in range(0,nSites):
 #bb=0
 ss=0
 for bb in range(0,nSites):
   psi_n = get_excitations(aa,ss)
   psi_m = get_excitations(bb,ss)
   for nn in range(0,len(psi_n)):
    for mm in range(0,len(psi_m)):
    
      s = '%d %d %d %d '%(aa,ss,bb,ss)
      s+= ' %d %d '%(nn,mm)
      
      psi_2 = (a_(aa+ss*nSites)*(c_(bb+ss*nSites)*(psi_m[mm])))
      val_2 = psi_2.scalarProd(psi_n[nn])/phi_pair.squareNorm()      
      s+= '% 4.2e '%(val_2.real)

      S_AC[aa+nn*nSites,bb+mm*nSites] = val_2.real      

      psi_2 = (c_(aa+ss*nSites)*(a_(bb+ss*nSites)*(psi_m[mm])))
      val_2 = psi_2.scalarProd(psi_n[nn])/phi_pair.squareNorm()      
      s+= '% 4.2e '%(val_2.real)

      if (1):
       print('\n\n\n',aa, bb, '  ', nn, mm)
       
       print(2.*np.sqrt(1./phi_pair.squareNorm())*psi_n[nn])
       #print(2.*np.sqrt(1./phi_pair.squareNorm())*psi_m[mm])
       print('\n< S_CA | x >')
       print(2.*np.sqrt(1./phi_pair.squareNorm())*psi_2)
       print(val_2)      
      S_CA[aa+nn*nSites,bb+mm*nSites] = val_2.real      
      
      psi_2 = 0.0*phi_pair
      for ii in range(0,nSites):
       for uu in range(0,2):
        for jj in range(0,nSites):
         psi_2 += H_T[ii,jj]*(a_(aa+ss*nSites)*(c_(ii+uu*nSites)*(a_(jj+uu*nSites)*(c_(bb+ss*nSites)*(psi_m[mm])))))
         #psi_2 += H_T[ii,jj]*(c_(ii+uu*nSites)*(a_(aa+ss*nSites)*(c_(bb+ss*nSites)*(a_(jj+uu*nSites)*(psi_m[mm])))))
         #psi_2 += H_T[ii,jj]*(c_(bb+ss*nSites)*(a_(aa+ss*nSites)*(c_(ii+uu*nSites)*(a_(jj+uu*nSites)*(psi_m[mm])))))
      
      for ii in range(0,nSites):
        psi_2 += H_U[ii]*(a_(aa+ss*nSites)*(n_(ii)*(n_(ii+nSites)*(c_(bb+ss*nSites)*(psi_m[mm]))))) 
      val_2 = psi_2.scalarProd(psi_n[nn])/phi_pair.squareNorm()
      
      s+= '% 4.2e '%(val_2.real)

      H_AC[aa+nn*nSites,bb+mm*nSites] = val_2.real      




      psi_2 = 0.0*phi_pair
      for ii in range(0,nSites):
       for uu in range(0,2):
        for jj in range(0,nSites):
         #psi_2 += H_T[ii,jj]*(c_(aa+ss*nSites)*(a_(bb+ss*nSites)*(c_(ii+uu*nSites)*(a_(jj+uu*nSites)*(psi_m[mm])))))
         psi_2 += H_T[ii,jj]*(c_(aa+ss*nSites)*(c_(ii+uu*nSites)*(a_(jj+uu*nSites)*(a_(bb+ss*nSites)*(psi_m[mm])))))
         
      for ii in range(0,nSites):
        psi_2 += H_U[ii]*(c_(aa+ss*nSites)*(n_(ii)*(n_(ii+nSites)*(a_(bb+ss*nSites)*(psi_m[mm]))))) 
      val_2 = psi_2.scalarProd(psi_n[nn])/phi_pair.squareNorm()
      
      s+= '% 4.2e '%(val_2.real)
      
      if (1):
       
       print('\n< H_CA | x >')
       #print(2.*np.sqrt(1./phi_pair.squareNorm())*psi_m[mm])
       #print(2.*np.sqrt(1./phi_pair.squareNorm())*psi_n[nn])
       print(2.*np.sqrt(1./phi_pair.squareNorm())*psi_2)
       print(val_2)
      H_CA[aa+nn*nSites,bb+mm*nSites] = val_2.real      

      
      print(s)
      fileo.write(s+'\n')

      #print('  %d  %d '%(nn,mm), end='')
      #print( (' % 4.6e  % 4.6e '%(val_2.real,val_2.imag)), end='\n')

fileo.close()

np.set_printoptions(precision=2)
print('\nS_AC')
print(S_AC)
print('\nS_CA')
print(S_CA)
print('\nH_AC')
print(H_AC)
print('\nH_CA')
print(H_CA)
print('')

