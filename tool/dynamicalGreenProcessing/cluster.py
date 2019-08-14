
import sys

if len(sys.argv)!=2:
  print("example:\n$ python cluster.py 12")
  sys.exit()


NN=int(sys.argv[1])

for ii in range(NN):
  print '#',
  for jj in range(NN):
    print ' %3d' % (NN*NN - (ii+1)*NN +jj),
  print ''
  
