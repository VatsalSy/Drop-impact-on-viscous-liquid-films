# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import sys

nGFS = 5000
ci = int(sys.argv[1])
Rhor = float(sys.argv[2])
Ohd = float(sys.argv[3])
Ohf = float(sys.argv[4])
Bond  = float(sys.argv[5])

folder = 'bview'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

name = "%d_getEnergy.dat" % ci

if os.path.exists(name):
    print("File %s found! New data will be appended to the file" % name)
for ti in range(nGFS):
    t = 0.01 * ti
    place = "intermediate/snapshot-%5.4f" % t
    ImageName = "bview/%8.8d.png" % int(1e3*t)
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        if os.path.exists(ImageName):
            print("Image %s found!" % ImageName)
        else:
            exe = "./getVideo %s %s" % (place, ImageName)
            os.system(exe)
            print(("Doing %d of %d" % (ti, nGFS)))
            exe = "./getEnergyAxi %s %s %s %s %s %s" % (place, name, Rhor, Ohd, Ohf, Bond)
            os.system(exe)
            
