

NN=12

for ii in range(NN):
  print '#',
  for jj in range(NN):
    print ' %3d' % (NN*NN - (ii+1)*NN +jj),
  print ''
  
