#!/usr/bin/python
#import numpy as np
#from scipy import linalg as LA
#from copy  import deepcopy
#import sys, os
#from scipy.linalg import eigh,eig
import sys, os
import numpy as np
import re
from ctypes import cdll, c_int, c_double, pointer, c_char_p

full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)



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
  if (n_file <1):
    print 'error: no input files.\nexample\n$ mergeOutput.py output/zvo_nCHAm_nAHCm_00*\n'
    exit()
    

  #use compiled version if present:
  if( os.path.isfile(pythonPathCode+'/mergeOutput.so')):
    c_ArrayDoubleN = c_double * (n_exc*n_exc*Nsite)

    phys_CA_averaged = c_ArrayDoubleN()
    phys_AC_averaged = c_ArrayDoubleN()
    phys_CHA_averaged = c_ArrayDoubleN()
    phys_AHC_averaged = c_ArrayDoubleN()
    c_NExcitation = c_int(n_exc)
    c_W = c_int(W)
    c_L = c_int(L)
    c_ncomp = c_int(n_exc*n_exc*Nsite*4)

    argList = (c_char_p * len(fileIn))()
    argList[:] = fileIn
    
    #Loading and using our home made module:
    lib1 = cdll.LoadLibrary(pythonPathCode+'/mergeOutput.so')
    lib1.mergeOutput(c_ncomp, len(fileIn), argList, c_NExcitation, c_L, c_W, 
                phys_CA_averaged,phys_AC_averaged,phys_CHA_averaged,phys_AHC_averaged)
    
    def convert_c2numpy(phys_averaged):
      numpy_averaged = np.zeros(n_exc*n_exc*Nsite)
      numpy_averaged[:] = phys_averaged[:]
      numpy_reshaped = numpy_averaged.reshape((Nsite,n_exc,n_exc))
      #print test[0,0,:]
      numpyOut = np.zeros([n_exc,n_exc,W,L], dtype='f')
      for rj in range(Nsite):
        rjx = X(rj)
        rjy = Y(rj)
        numpyOut[:,:,rjx,rjy] = numpy_reshaped[rj,:,:]
      return numpyOut

    nCAm_up = convert_c2numpy(phys_CA_averaged)
    nACm_up = convert_c2numpy(phys_AC_averaged)
    nCHAm_up = convert_c2numpy(phys_CHA_averaged)
    nAHCm_up = convert_c2numpy(phys_AHC_averaged)
    
  else:
    
    files = []
    for fileName in fileIn:
      fileTmp = openSafe(fileName, 'r')
      files.append(fileTmp)
    
    factor = 1.0/n_file
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
              rrj = int(fileLines[0][0])
              nn = int(fileLines[0][1])
              mm = int(fileLines[0][2])
              
              phys_CA = phys_AC = phys_CHA = phys_AHC = 0.0
              for jj in range(0,n_file):
                phys_CA  += float(fileLines[jj][3])
                phys_AC  += float(fileLines[jj][4])
                phys_CHA += float(fileLines[jj][5])
                phys_AHC += float(fileLines[jj][6])
              
              NN+=1
              rjx = X(rrj)
              rjy = Y(rrj)
   
              nCAm_up[nn,mm,rjx,rjy] = factor*phys_CA
              nACm_up[nn,mm,rjx,rjy] = factor*phys_AC
              nCHAm_up[nn,mm,rjx,rjy] = factor*phys_CHA
              nAHCm_up[nn,mm,rjx,rjy] = factor*phys_AHC
            
            if (isNothing):
              break
            #if NN%10000==0:
            #  print NN#, '/',len(matrix_list)

    #print '\nfft_CA';  nCAm_k  = np.fft.fft2(nCAm_up)
    #print '\nfft_AC';  nACm_k  = np.fft.fft2(nACm_up)
    #print '\nfft_CHA'; nCHAm_k = np.fft.fft2(nCHAm_up)
    #print '\nfft_AHC'; nAHCm_k = np.fft.fft2(nAHCm_up)
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

if __name__ == "__main__":
   main()
