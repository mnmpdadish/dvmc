### greenone.def ###

def main():
  
  Ns = 16
  Np = Ns*Ns*2
  f = open('greenone.def','w')
  f.write(
    "===============================\n"+
    "NCisAjs         "+str(Np)+"\n"+
    "===============================\n"+
    "======== Green functions ======\n"+
    "===============================\n"
    )
  for ii in range(Ns):
    for ss in range(2):
      for jj in range(Ns):
        f.write(
          '%5s %6s %6s %6s \n' % (str(ii), str(ss), str(jj), str(ss)) )
  f.close()

main()
