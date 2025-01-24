# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import sys

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'']

def gettingFacets(filename, tracer):
    print('Getting facets values')
    if tracer == 1:
        exe = ["./getFacet1", filename]
    else:
        exe = ["./getFacet2", filename]

    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e1):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2,z2)))
                    segs.append(((-r1, z1),(-r2,z2)))
                    skip = True
    print('Got facets values for tracer')
    return segs

def gettingfield(filename):
    print('Field values')
    exe = ["./getData_InsideDrop", filename, str(zmin), str(0), str(zmax), str(rmax), str(nz), str(Ohd), str(Ohf), str(hf)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    Rtemp, Ztemp, f1temp, f2temp, D2temp, Omegatemp, Utemp, Vtemp, veltemp = [] , [], [], [], [], [], [], [], []

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Ztemp.append(float(temp3[0]))
            Rtemp.append(float(temp3[1]))
            f1temp.append(float(temp3[2]))
            f2temp.append(float(temp3[3]))
            D2temp.append(float(temp3[4]))
            Omegatemp.append(float(temp3[5]))
            Utemp.append(float(temp3[6]))
            Vtemp.append(float(temp3[7]))
            veltemp.append(float(temp3[8]))

    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    f1 = np.asarray(f1temp)
    f2 = np.asarray(f2temp)
    D2 = np.asarray(D2temp)
    Omega = np.asarray(Omegatemp)
    U = np.asarray(Utemp)
    V = np.asarray(Vtemp)
    vel = np.asarray(veltemp)

    nr = int(len(R)/nz)

    print("nr is %d" % nr)
    R.resize((nz, nr))
    Z.resize((nz, nr))
    f1.resize((nz, nr))
    f2.resize((nz, nr))
    D2.resize((nz, nr))
    Omega.resize((nz, nr))
    U.resize((nz, nr))
    V.resize((nz, nr))
    vel.resize((nz, nr))

    print('Got Field values')
    return R, Z, f1, f2, D2, Omega, U, V, vel
# ------------------------------------------------------------------------------

nGFS = 2500
We = float(sys.argv[1])
Ohd= sys.argv[2]
Ohf = sys.argv[3]
hf = float(sys.argv[4])
print("Input Parameters: We = %f, Ohd = %s, Ohf = %s" % (We, Ohd, Ohf))

lw = 3

folder = 'TracerD2_InsideDrop'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

rmin, rmax, zmin, zmax = [-2.0, 2.0, -hf, 3.0-hf]

LEVEL = 10
nz = 2**(LEVEL)

for ti in range(nGFS):
    t = 0.01 * ti
    place = "intermediate/snapshot-%5.4f" % t
    name = "%s/%8.8d.png" %(folder, int(t*1e3))
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        if os.path.exists(name):
            print("Image %s found!" % name)
        else:
            facets1 = gettingFacets(place, 1)
            if (len(facets1)):
                facets2 = gettingFacets(place, 2)
                R, Z, f1, f2, D2, Omega, U, V, vel = gettingfield(place)
                speed = np.sqrt(U**2 + V**2)

                fig, ax = plt.subplots()
                fig.set_size_inches(19.20, 10.80)
                rc('axes', linewidth=2)
                # plt.xticks(fontsize=30)
                # plt.yticks(fontsize=30)

                ax.plot([0, 0], [zmin, zmax],'-.',color='grey',linewidth=lw)
                # ax.plot([rmin, rmax], [0, 0],'-',color='grey',linewidth=lw/2)
                ax.plot([rmin, rmin], [zmin, zmax],'-',color='black',linewidth=lw)
                ax.plot([rmin, rmax], [zmin, zmin],'-',color='black',linewidth=lw)
                ax.plot([rmin, rmax], [zmax, zmax],'-',color='black',linewidth=lw)
                ax.plot([rmax, rmax], [zmin, zmax],'-',color='black',linewidth=lw)

                ## omega
                cntrl1 = ax.imshow(D2, interpolation='bilinear', cmap="hot_r", origin='lower', extent=[0, -rmax, zmin, zmax], vmax = 1, vmin = -3)

                ## V
                maxs = speed.max()
                if maxs > 0.:
                    MaxVel = np.sqrt(We)
                    cntrl2 = ax.imshow(speed, interpolation='bilinear', cmap="Blues", origin='lower', extent=[0, rmax, zmin, zmax], vmax = MaxVel, vmin = 0.)

                line_segments1 = LineCollection(facets1, linewidths=4, colors='green', linestyle='solid')
                ax.add_collection(line_segments1)
                line_segments2 = LineCollection(facets2, linewidths=4, colors='green', linestyle='solid')
                ax.add_collection(line_segments2)

                ax.set_xlabel(r'$\mathcal{R}$', fontsize=40)
                ax.set_ylabel(r'$\mathcal{Z}$', fontsize=40)
                ax.set_aspect('equal')
                ax.set_xlim(rmin, rmax)
                ax.set_ylim(zmin, zmax)

                ax.set_title(r'$t/t_{\gamma} = %3.2f$' % t, fontsize=30)

                # l, b, w, h = ax.get_position().bounds
                # cb1 = fig.add_axes([l-0.16, b, 0.03, h])
                # c1 = plt.colorbar(cntrl1,cax=cb1,orientation='vertical')
                # c1.ax.tick_params(labelsize=30)
                # c1.set_label(r'$log_{10}\left(2\mathcal{O}h\|D_{ij}\|\right)$',fontsize=40)
                # c1.ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

                # cb2 = fig.add_axes([l+w+0.01, b, 0.03, h])
                # c2 = plt.colorbar(cntrl2,cax=cb2,orientation='vertical')
                # c2.ax.tick_params(labelsize=30)
                # c2.set_label(r'$\|V_i\|$',fontsize=40)
                # c2.ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
                # plt.show()


                ax.axis('off')

                plt.savefig(name, bbox_inches="tight")
                plt.close()
            else:
                print("Problem in the available file %s" % place)

    print(("Done %d of %d" % (ti+1, nGFS)))
