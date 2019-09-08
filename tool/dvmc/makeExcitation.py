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

import re

R1x = [0,1,2,3,4,5]
#R1x = [0,1,2]
R1y = [0]
R2x = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]
#R2x = [-1,0,1]
R2y = [0]


StdFace =  open('StdFace.def').read().replace(' ','')
L=W=1
if(StdFace.find('L=')>=0): L = int(re.compile('L=([0-9]*)').findall(StdFace)[0])
if(StdFace.find('W=')>=0): W = int(re.compile('W=([0-9]*)').findall(StdFace)[0])


def WriteExcitation():

  spectrumpara = open('spectrumpara.def').read()
  for line in  spectrumpara.split('\n'):
    if len(line)>0:
     if line[0]!='#':
      term = line.split()
      if(term[0][:]=='dr1_x'):  dr1_x = ReadRange(term[1])
      if(term[0][:]=='dr1_y'):  dr1_y = ReadRange(term[1])
      if(term[0][:]=='dr2_x'):  dr2_x = ReadRange(term[1])
      if(term[0][:]=='dr2_y'):  dr2_y = ReadRange(term[1])
  
  s = ''
  s+= '%5d %6d %2d %6d %2d \n' % (0,0,0,0,0) 
  s+= '%5d %6d %2d %6d %2d \n' % (1,0,0,0,0) 
  NN=2
  s4= ''
  s3= ''
  s2= ''
  pairList_s3 = []
  pairList_s4 = []
  for r1_x in dr1_x:
   for r1_y in dr1_y:
    for r2_x in dr2_x:
     for r2_y in dr2_y:
      if ((r2_x%L!=0) or (r2_y%W!=0)):  #condition to prevent unallowed excitation
       if ((r1_x%L,r1_y%W,r2_x%L,r2_y%W) not in pairList_s4):   #condition to prevent redundant excitation
        NN+=1
        s4+= '%5d %6d %2d %6d %2d \n' % (4,r1_x%L,r1_y%W,r2_x%L,r2_y%W)
        pairList_s4.append((r1_x%L,r1_y%W,r2_x%L,r2_y%W))
        if ((r2_x%L,r2_y%W,r1_x%L,r1_y%W) not in pairList_s3):   #condition to prevent redundant excitation
         NN+=1#2
         s3+= '%5d %6d %2d %6d %2d \n' % (3,r1_x%L,r1_y%W,r2_x%L,r2_y%W)
         #s2+= '%5d %6d %2d %6d %2d \n' % (2,r1_x%L,r1_y%W,r2_x%L,r2_y%W)   # careful, if you restore this line, also restore NN+=2 line
         pairList_s3.append((r1_x%L,r1_y%W,r2_x%L,r2_y%W))
         
  f = open('excitation.def','w')
  f.write(
    "===============================\n"+
    "NExcitation       "+str(NN)+"\n"+
    "L                 "+str(L)+"\n"+
    "W                 "+str(W)+"\n"+
    "==================Excitations==\n"
    )
  
  f.write(s+s4+s3) 
  f.close()

def ReadRange(inputStr):
  error=0
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
      print 'sepated by comma (no spaces).'
      print 'Pairs must contains only one ":"' 
      print 'example:'
      print '0:3,5,6,8:11'
      print 
      print 'will be interpreted as the list'
      print '[0,1,2,3,5,6,8,9,10,11]'
      exit()
  return rangeOut


WriteExcitation()
