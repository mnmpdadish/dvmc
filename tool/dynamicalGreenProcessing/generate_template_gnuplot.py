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

full_path = os.path.realpath(__file__)
pythonPathCode, file1 = os.path.split(full_path)

StdFace = open('StdFace.def').read().replace(' ','')
U=0.
L=W=1
if(StdFace.find('L=')>=0): L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
if(StdFace.find('W=')>=0): W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])
if(StdFace.find('U=')>=0): U = float(re.compile('U=([0-9.]*)').findall(StdFace)[0])

def main():

  spectrumpara = open('spectrumpara.def').read()

  for line in  spectrumpara.split('\n'):
    if len(line)>0:
     if line[0]!='#':
      term = line.split()
      if term[0]=='w_min' : w_min = float(term[1])
      if term[0]=='w_max' : w_max = float(term[1])
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

  
  
  def replaceTemplateFileValues(templateFileName):
  
    gnuplotString = open(pythonPathCode+'/'+templateFileName).read()
    #w_mid = 0.5*(w_min+w_max)
    gnuplotString = gnuplotString.replace('w_min_data = -13.0','w_min_data =% 4.5f'%(w_min))
    gnuplotString = gnuplotString.replace('w_max_data =  13.0','w_max_data =% 4.5f'%(w_max))
    gnuplotString = gnuplotString.replace('w_min = -8.0','w_min =% 4.5f'%(w_min))
    gnuplotString = gnuplotString.replace('w_max =  8.0','w_max =% 4.5f'%(w_max))
    gnuplotString = gnuplotString.replace('kRange = 32','kRange =% 4.5f'%(len(kPath1)-1))
    gnuplotString = gnuplotString.replace('Nw = 1500','Nw =% 4.5f'%(Nw))
    
    fileOut = open('./'+templateFileName[9:],'w')
    fileOut.write(gnuplotString)
    fileOut.close()
    return;
  
  replaceTemplateFileValues('template_plot_allAkw.gp')
  replaceTemplateFileValues('template_plot_singleAkw.gp')
  
    
  exit()
  
#######################################################################
######################    SUB ROUTINES    #############################
#######################################################################



def ReadRange(inputStr):
  rangeOut = []
  range1 = inputStr.split(',')
  #print range1
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
