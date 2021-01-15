#!/usr/bin/python3
from numpy import *

def main():

   flag_symmetrize_index = False # f_ij  !=  f_ji
   #flag_symmetrize_index = True #v_ij   =  v_ji
   
   #symmetries = ['0 1 2 3 4 5 6 7 8',
   #              '3 4 2 1 5 8 7 6 0',
   #              '1 7 2 3 0 8 4 6 5']
   
   symmetries = ['1 0 3 2',
                 '2 3 0 1']
                 
   verbose = 1   
   set_printoptions(precision=5)
   find_f_ij(symmetries, flag_symmetrize_index, verbose)
   
   exit()
   
###########################################################################################################
#####################################  END OF THE PROGRAM - ###############################################
###########################################################################################################

def charHexa(digit):
   characters = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!?' #any cluster should not have above 36 sites, even in 2050 xD
   if digit < 64 and digit >= 0:
      char = characters[digit]
   else:
      char = '.'
   return char

def find_f_ij(symmetries, flag_symmetrize_index,verbose=2):
   gL = len(symmetries[0].split())
   
   symmetryRuleList = []
   for sym in symmetries:
     print(sym)
     sym_int = [int(i) for i in sym.split()] 
     
     symmetryRuleList.append(sym_int)
   
   def swapCouple(couple):
      return (couple[1],couple[0])

   f_ij_indexBasic = {}
   f_ij_index = {}
   for ii in range(0,gL):
      for jj in range(0,gL):
         f_ij_indexBasic[ii,jj] = ij(jj,ii) 
         couple = (ii,jj)
         
         i = ii
         j = jj
         if flag_symmetrize_index:
           if ii>jj:
             i = jj
             j = ii
             
         f_ij_index[ii,jj] = ij(i,j)
         
   
   if verbose>0:
      print('\ngeneral \nf_ij =')
      print(print_f_ij_indexToString(f_ij_indexBasic,gL))
      
      print('\ni->j, j->i symmetry: \nf_ij =')
      print(print_f_ij_indexToString(f_ij_index,gL))
   
   changed = True
   while changed:
      [changed,symmetryRuleList2] = symmetrizeElementOneByOne(symmetryRuleList,f_ij_index,gL,verbose)

   
   if verbose>0:
      print('\nsymmetrized\nf_ij =')
      print(print_f_ij_indexToString(f_ij_index,gL))

   listOfIndepGreenFound = {}
   N=0
   for ii in range(0,gL):
      for jj in range(0,gL):
         if f_ij_index[ii,jj].string() not in listOfIndepGreenFound:
            N+=1
            #print ii,jj
            listOfIndepGreenFound[f_ij_index[ii,jj].string()] = (ii,jj) 
   
   #print(N)
   

   return 

def symmetrizeElementOneByOne(symmetryRuleList,f_ij_index,gL,verbose):
   zero = ij(0,0)
   symmetryRuleList2 = []


   for sym in symmetryRuleList:
      symmetryRuleList2.append(sym)
      symmetryRuleList2
               
      for i in range(0,gL):
         for j in range(0,gL):
            permutation_i= i
            permutation_j= j
            nn = 0
            isFirst = True
            
            while ((permutation_i != i) or (permutation_j != j)) or isFirst: #while didn't loop on orbit
               nn+=1
               permutation_i= sym[permutation_i]
               permutation_j= sym[permutation_j]
               factor = 1.0               
               isConj = False

               
               if not(f_ij_index[permutation_j,permutation_i].isEqual(f_ij_index[j,i])): 
                 ij0 = f_ij_index[j,i].string()
                 ij1 = f_ij_index[permutation_j,permutation_i].string()
                 
                 ####this portion of code ensure that always the lowest name (in alphabetical order) is kept
                 if ij0 < ij1:
                    f_ij_index[permutation_j,permutation_i].equal(f_ij_index[j,i])
                 else:
                    f_ij_index[j,i].equal(f_ij_index[permutation_j,permutation_i])
                 ####
                 
                 return [True, symmetryRuleList2]
               isFirst = False
               if nn>100: 
                  print('error, too many iteration')
                  exit()
   return [False, symmetryRuleList2]        


class ij:
   def __init__(self,i,j):
      self.i = i
      self.j = j
   def equal(self,another):
      self.i = another.i
      self.j = another.j
   def isEqual(self,another):
      return (self.i == another.i) and (self.j == another.j)
   def string(self):
      return charHexa(self.i) + charHexa(self.j)


def print_f_ij_indexToString(f_ij_index,gL):
   maxLen = 1
   s1 = ''
   for jj in range(0,gL):
      for ii in range(0,gL):
         maxLen = max(maxLen,len(f_ij_index[ii,jj].string()))
   
   for jj in range(0,gL):
      for ii in range(0,gL):
         var=1
         if (f_ij_index[ii,jj].string())[0]!='-':
            s1 += ' '
            var=0 
         s1 += f_ij_index[ii,jj].string()
         for kk in range(0,maxLen-len(f_ij_index[ii,jj].string())+var):
            s1 += ' '
      if jj is not gL-1: s1 += '\n'
   return s1



###########################################################################################################

if __name__ == "__main__":
   main()
