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
#   misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.

import numpy as np
import numpy.linalg as la

def main():  
  np.set_printoptions(precision=3)
  
#  tc  = np.array([[ 0., -t , -tp, -t ],
#                  [-t ,  0., -t , -tp],
#                  [-tp, -t ,  0., -t ],
#                  [-t , -tp, -t ,  0.]])
  H = np.array([[1,2,6],
                [3,4,0],
                [1,2,1.]])

  H2 = np.array([[4,1,3],
                 [4,3,1],
                 [4,1,2.]])

  e,u = la.eig(H)
  
 
  print
  print 'H'
  print H
  print
  print 'multiplication'
  print np.dot(H,H2)
  print
  print 'inv'
  print la.inv(H)
  print
  print 'eig'
  print e
  print
  print u

  print u[0,1]

main()
