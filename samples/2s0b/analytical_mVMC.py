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


      
      
################################################
##### main program #############################
################################################

# This only check for small 1D clusters

nSites = 2
Read = 1

# defining an Hamiltonian:
H_T = np.array([[   0.,  -1.],
                [  -1.,   0.]])

H_U = np.array([  8.,  8.])

# assigning the f_ij:
f_ = np.zeros((nSites,nSites), dtype=complex)
if(Read):
  print('\nreading file ./output/zqp_opt.dat')
  data = ReadFile("./output/zqp_opt.dat")
  data_list = data[0].split()
  print(len(data_list))

  for ii in range(nSites):
    for jj in range(nSites):
      data_ii_re = (data[0].split())[len(data_list)-3-3*(ii+jj*nSites)]
      data_ii_im = (data[0].split())[len(data_list)-2-3*(ii+jj*nSites)]
      f_[ii][jj] += float(data_ii_re) + 1.j * float(data_ii_im)
      #print(ii,jj,f_[ii][jj])
      
else:
  f_[0][0] = 1.0
  f_[0][1] = 2.0
  f_[1][0] = 3.0
  f_[1][1] = 4.0 # all arbitrary for now

print('\ni j  f_ij')
for ii in range(nSites):
  for jj in range(nSites):
    print(ii,jj,'% 4.3f' % f_[ii][jj].real) # limited to real f_ij for now




# phi_0 is the vacuum:
phi_0 = qState(nSites,1.0)
print('\n| phi_0 > (this is the vacuum):')
print(phi_0)

# |phi> = phi_pair = 0, it is a recipient for now:
phi_pair = qState(nSites, 0.0)

# | phi >  = sum_{ij} f_{ij} c_{i down} c_{j up} | phi_0 >  :
for ii in range(0,nSites):
 for jj in range(0,nSites):
  phi_pair += f_[ii][jj] * (c_(nSites+ii) * (c_(jj) * phi_0))
print('\n| phi >  = sum_{ij} f_{ij} c_{i down} c_{j up} | phi_0 >')
print(phi_pair)

# | psi >  = H |phi> :
psi = qState(nSites,0.0)

for ii in range(0,nSites):
 for ss in range(0,2): # ss for spin
  for jj in range(0,nSites):
   psi += H_T[ii,jj]*(c_(ii+ss*nSites)*(a_(jj+ss*nSites)*(phi_pair)))

for ii in range(0,nSites):
  #psi += H_U[ii]*(phi_pair) 
  psi += H_U[ii]*(n_(ii)*(n_(ii+nSites)*(phi_pair)))

print('\n| psi >  = H |phi>')
print(psi)

# ground state energy: E = <phi| H |phi> / <phi | phi>
E = psi.scalarProd(phi_pair)/phi_pair.squareNorm()

print('\nE = <phi| H |phi> / <phi | phi> = % 4.8f' % E.real) # E is supposed to be purely real


# -1.656854249492419440e+00
# -1.656854249492395459e+00


#def flavor(i,nSites):
#  site = i%nSites
#  spin = i//nSites
#  return '%1d %1d '%(site,spin)   

#def pretty(float1):
#  if abs(float1)<1e-8:
#    return '  . '#%(float1)   
#  return '% 4.3e '%(float1)   

