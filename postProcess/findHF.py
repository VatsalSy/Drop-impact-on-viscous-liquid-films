# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import sys

nGFS = 10000
ci = int(sys.argv[1])

name = "%d_getHF.dat" % ci

if os.path.exists(name):
    print("File %s found! New data will be appended to the file" % name)
for ti in range(nGFS):
    t = 0.01 * ti
    place = "intermediate/snapshot-%5.4f" % t
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
      print(("Doing %d of %d" % (ti, nGFS)))
      exe = "./findHF %s %s" % (place, name)
      os.system(exe)
            
