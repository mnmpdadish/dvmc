#!/usr/bin/python3
from numpy import *

def main():

   L=4
   W=4
   Lsub=2
   Wsub=2

   verbose = 1   
   set_printoptions(precision=5)
   
   for ii in range(16):
    for jj in range(16):
     (fold_x,fold_y) = fold_vector(ii,2,2,4,4)
     print '%2d %2d  - %2d %2d         - %2d %2d' %(ii,jj,fold(ii,(fold_x,fold_y),2,2,4,4), fold(jj,(fold_x,fold_y),2,2,4,4),  fold_x, fold_y)
   
   exit()
   
###########################################################################################################
#####################################  END OF THE PROGRAM - ###############################################
###########################################################################################################

def xy(i, W, L):
   x = i % W
   y = i // W
   return (x, y)

def fold_vector(index, Wsub, Lsub, W, L):
   (x,y) = xy(index, W, L)
   
   fold_x = (x//Wsub)*Wsub
   fold_y = (y//Lsub)*Lsub
   
   return (fold_x,fold_y)

def fold(index, (fold_x, fold_y), Wsub, Lsub, W, L):
   (x,y) = xy(index, W, L)
   
   x2 = (x - fold_x) % W
   y2 = (y - fold_y) % L
   new_index = x2 + y2*W
   return new_index


if __name__ == "__main__":
   main()
