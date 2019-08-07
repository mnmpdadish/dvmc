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
  pairList = []
  for r1_x in dr1_x:
   for r1_y in dr1_y:
    for r2_x in dr2_x:
     for r2_y in dr2_y:
      if ((r2_x!=0) or (r2_y!=0)):
       NN+=1
       s4+= '%5d %6d %2d %6d %2d \n' % (4,r1_x,r1_y,r2_x,r2_y)
       if ((r2_x,r2_y,r1_x,r1_y) not in pairList):
        NN+=1#2
        s3+= '%5d %6d %2d %6d %2d \n' % (3,r1_x,r1_y,r2_x,r2_y)
        #s2+= '%5d %6d %2d %6d %2d \n' % (2,r1_x,r1_y,r2_x,r2_y)
        pairList.append((r1_x,r1_y,r2_x,r2_y))
  
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
  rangeOut = []
  range1 = inputStr.split(',')
  print range1
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

WriteExcitation()
