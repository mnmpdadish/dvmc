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
     (x,y)   = xy(ii, (W, L))
     (fold_x,fold_y) = fold_vector((x,y), (Wsub, Lsub))
     print '%2d %2d  - %2d %2d ' %(ii,jj,fold(ii,(fold_x,fold_y),(W,L)), fold(jj,(fold_x,fold_y),(W,L)))
   
   exit()
   
###########################################################################################################
#####################################  END OF THE PROGRAM - ###############################################
###########################################################################################################

def xy(i, (W, L)):
   x = i % W
   y = i // W
   return (x, y)

def fold_vector((x,y), (Wsub, Lsub)):
   
   fold_x = (x//Wsub)*Wsub
   fold_y = (y//Lsub)*Lsub
   
   return (fold_x, fold_y)

def fold(index, (fold_x, fold_y), (W, L)):
   (x,y) = xy(index, (W, L))
   
   new_x = (x - fold_x) % W
   new_y = (y - fold_y) % L
   new_index = new_x + new_y*W
   return new_index


if __name__ == "__main__":
   main()
