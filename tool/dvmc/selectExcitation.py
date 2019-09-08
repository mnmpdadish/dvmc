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

import re,sys


if len(sys.argv)!=3:
  print("example (min=1,max=4):\n$ selectExcitation.py 1 4")
  sys.exit()

min1=int(sys.argv[1])
max1=int(sys.argv[2])

print min1, max1


def SelectExcitation():

  spectrumpara = open('spectrumpara.def').read()
  for line in  spectrumpara.split('\n'):
    if len(line)>0:
     if line[0]!='#':
      term = line.split()
      if(term[0][:]=='dr1_x'):  dr1_x = ReadRange(term[1])
      if(term[0][:]=='dr1_y'):  dr1_y = ReadRange(term[1])
      if(term[0][:]=='dr2_x'):  dr2_x = ReadRange(term[1])
      if(term[0][:]=='dr2_y'):  dr2_y = ReadRange(term[1])

  #s = ''
  #s+= '%5d %6d %2d %6d %2d \n' % (0,0,0,0,0) 
  #s+= '%5d %6d %2d %6d %2d \n' % (1,0,0,0,0) 
  NN=2
  NNs4=0
  NNs3=0
  #s4= ''
  #s3= ''
  #s2= ''
  exc_choice_s3 = []
  exc_choice_s4 = []
  pairList = []
  for r1_x in dr1_x:
   for r1_y in dr1_y:
    for r2_x in dr2_x:
     for r2_y in dr2_y:
      if ((r2_x!=0) or (r2_y!=0)):
       x1 = abs(r1_x)
       x2 = abs(r2_x)
       y1 = abs(r1_y)
       y2 = abs(r2_y)

       if ((x1<=max1) and (x2<=max1) and (y1<=max1) and (y2<=max1) ):
        if ((x1>=min1) and (x2>=min1)):
         exc_choice_s4.append(NNs4)
         #print NN, '   ', r1_x,r1_y,r2_x,r2_y
         #print '%10d  %5d %6d %2d %6d %2d ' % (NN,4,r1_x,r1_y,r2_x,r2_y)
         NN+=1
       
       NNs4+=1

       #s4+= '%5d %6d %2d %6d %2d \n' % (4,r1_x,r1_y,r2_x,r2_y)
       if ((r2_x,r2_y,r1_x,r1_y) not in pairList):
        if ((x1<=max1) and (x2<=max1) and (y1<=max1) and (y2<=max1) ):
         if ((x1>=min1) and (x2>=min1)  ):
          exc_choice_s3.append(NNs3)
          #print NN, '   ', r1_x,r1_y,r2_x,r2_y
          #print '%10d  %5d %6d %2d %6d %2d ' % (NN,3,r1_x,r1_y,r2_x,r2_y)
          NN+=1
       
        #exc_choice.append(NN)
        #NN+=1#2
        NNs3+=1
        #s3+= '%5d %6d %2d %6d %2d \n' % (3,r1_x,r1_y,r2_x,r2_y)
        #s2+= '%5d %6d %2d %6d %2d \n' % (2,r1_x,r1_y,r2_x,r2_y)
        pairList.append((r1_x,r1_y,r2_x,r2_y))

  #print exc_choice
  s='exc_choice   0,1,'
  ii=0

  for num in exc_choice_s4:
    s+='%d,'%(num+2)
  for num in exc_choice_s3:
    s+='%d,'%(num+NNs4+2)
  print s[:-1]
  
  print NN
  exit()
  
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
      print 'your definition make no sense'
      print 'Please input numbers or pairs of numbers'
      print 'sepated by comma. Pairs must contains only one ":"'
      print 'example:'
      print '0:3,5,6,8:11'
      exit()
  return rangeOut

SelectExcitation()

