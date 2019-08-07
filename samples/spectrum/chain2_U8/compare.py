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
#   misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.

import numpy as np
import os
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

def ReadFile(fileName):
  fileExist = os.path.isfile(fileName)
  if not fileExist:
    print('\nerror: file '+fileName+' does not exist.')   
    print('terminated.')
    exit()
   
  f = open(fileName, 'r')
  file1 = f.readlines()
  f.close
  
  file2 = []
  for lines in file1:
    file2.append(lines.strip())
  return file2

data1 = ReadFile("./output/green_comp_ave.dat")
data2 = ReadFile("./output/zvo_nCHAm_nAHCm_001.dat")

diff = 0.0

col = [6,7,8,9]
diffs = [0.0,0.0,0.0,0.0]
for column in col:
 for ii in range(len(data1)):
  if data1[ii][0] !='#':
    val1 = float((data1[ii].split())[column])
    val2 = float((data2[ii].split())[column-3])
    print ("%s   % 2.2f,        % 3.4f % 3.4f"% ( (data1[ii][0:14]), abs(val1-val2), val1, val2))
    diff += abs(val1 - val2 ) / (len(data1))
    diffs[column-6] += abs(val1 - val2 ) / (len(data1))
    
 print('')
print('SUMMARY:')
print('')
print('total error = %f' % diff)
print('')
print('<nCAm> error = %f' % diffs[0])
print('<nACm> error = %f' % diffs[1])
print('<nCHAm> error = %f' % diffs[2])
print('<nAHCm> error = %f' % diffs[3])

exit()




## analysis of the error as a funciton of sampling
y = np.transpose(np.array([[0.13169230769230786, 0.13169230769230786, 0.46769230769230896, 0.9796923076923081],
              [0.02326153846153849, 0.023261538461538434, 0.08123076923076919, 0.1668923076923077],
              [0.0020799999999999985, 0.002079999999999997, 0.007138461538461524, 0.014956307692307666],
              [0.002036923076923081, 0.0020369230769230803, 0.007560615384615395, 0.015212553846153875],
              [0.00023310769230768525, 0.0002331076923076852, 0.0008323692307692102, 0.001704492307692273]]))

x = np.array([50,500,5000,50000,500000])

#y = np.array([0.00722673893406,0.00316169828365,0.000943089430894,0.000237127371274,3.43270099368e-05])
fig = plt.figure()
ax = plt.axes()
plt.xlabel('log N')
plt.ylabel('log error')



plt.xlim([1,6])
#plt.ylim([0,0.008])
ax.plot(np.log(x)/np.log(10), np.log(y[0,:]) / np.log(10));
ax.plot(np.log(x)/np.log(10), np.log(y[1,:]) / np.log(10));
ax.plot(np.log(x)/np.log(10), np.log(y[2,:]) / np.log(10));
ax.plot(np.log(x)/np.log(10), np.log(y[3,:]) / np.log(10));

n = np.linspace(1, 6, 100)
x = 10**n
ax.plot(np.log(x)/np.log(10), np.log(0.25/np.sqrt(x)) / np.log(10));
ax.legend(['mVMC <n|CA|m>',
           'mVMC <n|AC|m>',
           'mVMC <n|CHA|m>',
           'mVMC <n|AHC|m>',
           '~1/sqrt(N)' ])

plt.savefig("figure_converge"+".pdf")
plt.show()
  
###############################################################
exit(0)

