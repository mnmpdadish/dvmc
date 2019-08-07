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
  for fileName in fileIn:
    fileTmp = openSafe(fileName, 'r')
    files.append(fileTmp)
  
  file_ref = files[0]
  
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
  for rj in range(Nsite):
    print '%d / %d' % (rj+1,Nsite)
    for n in range(n_exc):
      for m in range(n_exc):
      
        fileLines = []
        isComment = False
        isNothing = True
        for f in files:
          
          lineTmp = f.readline()
          #print len(lineTmp), lineTmp
          if len(lineTmp)>0:
            isNothing = False
            if lineTmp[0]=='#':
              isComment = True
          fileLines.append(lineTmp.split())
        
        if (not isComment):
          if (not isNothing):
            #ri = int(fileLines[0][0])
            rj = int(fileLines[0][0])
            nn = int(fileLines[0][1])
            mm = int(fileLines[0][2])
            
            phys_CA = phys_AC = phys_CHA = phys_AHC = 0.0
            for jj in range(0,n_file):
              phys_CA  += float(fileLines[jj][3])/n_file
              phys_AC  += float(fileLines[jj][4])/n_file
              phys_CHA += float(fileLines[jj][5])/n_file
              phys_AHC += float(fileLines[jj][6])/n_file
            
            NN+=1
            #rix = X(ri)
            #riy = Y(ri)
            rjx = X(rj)
            rjy = Y(rj)
 
            nCAm_up[nn,mm,rjx,rjy] += phys_CA
            nACm_up[nn,mm,rjx,rjy] += phys_AC
            nCHAm_up[nn,mm,rjx,rjy] += phys_CHA
            nAHCm_up[nn,mm,rjx,rjy] += phys_AHC
          
          if (isNothing):
            break
          #if NN%10000==0:
          #  print NN#, '/',len(matrix_list)

  #s = [W,L,W,L]
  #print '\nfft_CA';  nCAm_k = np.fft.fftn(nCAm_k,s)
  #print '\nfft_AC';  nACm_k = np.fft.fftn(nACm_k,s)
  #print '\nfft_CHA'; nCHAm_k = np.fft.fftn(nCHAm_k,s)
  #print '\nfft_AHC'; nAHCm_k = np.fft.fftn(nAHCm_k,s)
  #print nCAm_up

  print dirOutput+'nCAm'
  print dirOutput+'nACm'
  print dirOutput+'nCHAm'
  print dirOutput+'nAHCm'
  np.save(dirOutput+'nCAm',nCAm_up)
  np.save(dirOutput+'nACm',nACm_up)
  np.save(dirOutput+'nCHAm',nCHAm_up)
  np.save(dirOutput+'nAHCm',nAHCm_up)
  #f.close

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
  
if __name__ == "__main__":
   main()
