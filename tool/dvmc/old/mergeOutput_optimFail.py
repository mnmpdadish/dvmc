#!/usr/bin/python
#import numpy as np
#from scipy import linalg as LA
#from copy  import deepcopy
#import sys, os
#from scipy.linalg import eigh,eig
import sys, os
import numpy as np
import re

StdFace =  open('StdFace.def').read().replace(' ','')
L=W=1
if(StdFace.find('L=')>=0): L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
if(StdFace.find('W=')>=0): W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])

def X(ri):
  return ri%W

def Y(ri):
  return (ri/W)%L


def openSafe(fileName,mode):
  fileExist = os.path.isfile(fileName)
  if not fileExist:
    print '\nerror: file '+fileName+' does not exist.'   
    print 'terminated.'
    exit()
  return open(fileName, mode)


def ReadValue(valueName,stream):
  lines = stream.readline().split()
  #print lines[0]
  assert(lines[0]==valueName)
  return int(lines[1])
  

def ReadExcitation(fileName):

  f = openSafe(fileName, 'r')
  #file1 = f.readlines()
  
  NN=0
  f.readline()
  NExcitation = ReadValue('NExcitation',f)
  f.close()
  return NExcitation

def stringHasDot(test):
  return test.find('.')!=-1

def main():
  dirOutput = './output/'
  
  n_exc = ReadExcitation("./excitation.def")
  Nsite = W*L#exit()

  fileIn = []
  n_file = len(sys.argv[:])-1
  n_file2 = 0
  print sys.argv[:]
  print n_file
  
  for nn in range(0,n_file):
    n_file2 +=1
    fileIn.append(sys.argv[1+nn])
  
  print fileIn
  fileName1 =sys.argv[1]
  #fileOut = open(fileName1[:fileName1.find('_0')]+'_ave.dat', 'w' )
  
  files = []
  #for fileName in fileIn:
  #  fileTmp = openSafe(fileName, 'r')
  #  files.append(fileTmp)
  
  factor = 1./n_file
  
  end = 0
  NN=0
  
  nCAm_up  = np.zeros([n_exc,n_exc,W,L], dtype='f')
  nACm_up  = np.zeros([n_exc,n_exc,W,L], dtype='f')
  nCHAm_up = np.zeros([n_exc,n_exc,W,L], dtype='f')
  nAHCm_up = np.zeros([n_exc,n_exc,W,L], dtype='f')
  
  #nCAm_k  = np.zeros([n_exc,n_exc,W,L,W,L], dtype=complex)
  #nACm_k  = np.zeros([n_exc,n_exc,W,L,W,L], dtype=complex)
  #nCHAm_k = np.zeros([n_exc,n_exc,W,L,W,L], dtype=complex)
  #nAHCm_k = np.zeros([n_exc,n_exc,W,L,W,L], dtype=complex)
  
  #for ri in range(Nsite):
  for fileName in fileIn:
    print fileName
    with open(fileName) as f:
      for line in f:
      
        #print line,
        fileLines = []
        isComment = False
        isNothing = True
        
        if len(line)>0:
          isNothing = False
          if line[0]=='#':
            isComment = True
        
        lineTmp = line.split()
        
        if (not isComment):
          if (not isNothing):
            #print lineTmp[0]
            rj = int(lineTmp[0])
            nn = int(lineTmp[1])
            mm = int(lineTmp[2])
            
            NN+=1
            rjx = X(rj)
            rjy = Y(rj)
 
            nCAm_up[nn,mm,rjx,rjy] += factor*float(lineTmp[3])
            nACm_up[nn,mm,rjx,rjy] += factor*float(lineTmp[4])
            nCHAm_up[nn,mm,rjx,rjy] += factor*float(lineTmp[5])
            nAHCm_up[nn,mm,rjx,rjy] += factor*float(lineTmp[6])
          
          if (isNothing):
            break

  print dirOutput+'nCAm'
  print dirOutput+'nACm'
  print dirOutput+'nCHAm'
  print dirOutput+'nAHCm'
  np.save(dirOutput+'nCAm',nCAm_up)
  np.save(dirOutput+'nACm',nACm_up)
  np.save(dirOutput+'nCHAm',nCHAm_up)
  np.save(dirOutput+'nAHCm',nAHCm_up)
  #f.close

if __name__ == "__main__":
   main()
