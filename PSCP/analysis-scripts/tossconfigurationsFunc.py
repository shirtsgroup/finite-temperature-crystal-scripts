# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 10:17:49 2015

@author: nps5kd
"""

from __future__ import print_function
import sys
import optparse
import numpy as np
import math
import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties as FP
import pdb


def readgro(filename, na, nm, nam):
    # Read in atom coordinate data (starts on 3rd line of the .gro file)
    # acoordA and acoordB are the atom coordinates for polymorphs A and B

    file = open(filename, 'r')
    lines = file.readlines()    
    #lines = filter(None, (line.rstrip() for line in file))
    file.close()
    ngro = (len(lines)) / (na + 3)
    line = 2
    boxv = np.zeros([3, ngro])
    acoord = np.zeros((ngro, nm, 3, nam))
    for gro in range(int(ngro)):
        mol = 0
        while mol < nm:
            acounter = 0
            while acounter < nam:
                acoord[gro, mol, 0, acounter] = float(lines[line].split()[3])
                acoord[gro, mol, 1, acounter] = float(lines[line].split()[4])
                acoord[gro, mol, 2, acounter] = float(lines[line].split()[5])
                line += 1
                acounter += 1
            mol += 1    
        boxv[:, gro] = np.array(lines[line].split()).astype(np.float)
        line += 3        
    return acoord, boxv, ngro
    
def find_centroids(acoord, ngro, boxv, nm):
    # ====================================================================
    # Calculate the centroid and idealized centroid of each molecule
    # ====================================================================

    all_centroids = np.zeros([ngro, nm, 3, 1])
    ideal_centroids = np.zeros([nm, 3, 1])

    for mol in range(int(nm)):
        for gro in range(int(ngro)):
            all_centroids[gro, mol, :, 0] = np.average(acoord[gro, mol, :, :], axis=1)
        centroids = np.average(all_centroids, axis=0)

    # for now, we start with centroids are the averages over all molecules (i.e. the centroids)
    for mol in range(int(nm)):
        ideal_centroids[mol, :, 0] = centroids[mol, :, 0].copy()
    return all_centroids, ideal_centroids



def rmscalc(coord, coordPrime, gronum, nm, nam):
    # Quantify error of mapping (least root mean square)
    diff = np.zeros((nm, 3, nam))
    lrms = np.zeros((nm, nam))
    for mol in range(int(nm)):
        for mole in range(int(nam)):
            for co in range(3):
                diff[mol, co, mole] = abs(coord[gronum, mol, co, mole] - coordPrime[gronum, mol, co, mole])
    for mol in range(int(nm)):
        for mole in range(int(nam)):
            lrms[mol, mole] = math.sqrt(diff[mol, 0, mole] ** 2 + diff[mol, 1, mole] ** 2 + diff[mol, 2, mole] ** 2)
    return lrms

def writegrofile(output, lines, numframe, nm, nam, coord, boxv):
    outgrofile = open(output, 'w')
    # Write output .gro file
    # tacount is the total atom count
    # mcount is the total molecule count
    x = np.zeros(3)
    na = nm * nam
    molname = 'BNZ'
    for gro in range(int(numframe)):
#        outgrofile.write('Converted from ' + sys.argv[1] + ' by ' + sys.argv[2]+', molecule ' + str(gro) + '\n')
        outgrofile.write('Frame number' + str(gro) + '\n')    
        outgrofile.write('  ' + str(na) + '\n')
        tacount = 0
        mcount = 0
        for mol in range(int(nm)):
            for atom in range(int(nam)):
                # mcount is the current molecule tracker
                mcount = mol + 1
                tacount += 1
                xstr = []
                for d in range(3):
    #                x[d] = round(acoordBprime[mol,d,atom],8)                    
                    x[d] = round(coord[gro, mol, d, atom], 8)
                    xstr.append("%13.8f" % (x[d]))
                line = str(mcount).rjust(5) + molname + lines[tacount + 1].split()[1].rjust(7) + \
                       str(tacount).rjust(5) + xstr[0] + xstr[1] + xstr[2] + '\n'
                outgrofile.write(line)
        boxvline = ''    
        for d in range(len(boxv)):
            bline = ("%13.8f" % (round(boxv[d], 8)))
            boxvline += bline
        outgrofile.write(boxvline+'\n')  # make sure the .gro ends on a newline
    outgrofile.close()

def normalcompare(gro1, gro2):
    angle = np.zeros([4, ])
    for x in range(4):
        C11 = gro1[x, :, 0]
        C31 = gro1[x, :, 4]
        C51 = gro1[x, :, 8]
        C12 = gro2[x, :, 0]
        C32 = gro2[x, :, 4]
        C52 = gro2[x, :, 8]
        A13 = C31 - C11
        A15 = C51 - C11
        B13 = C32 - C12
        B15 = C52 - C12
        gro1vect = np.cross(A13, A15)
        gro2vect = np.cross(B13, B15)
        gro1vectnorm = np.linalg.norm(gro1vect)
        gro2vectnorm = np.linalg.norm(gro2vect)
        cost = np.dot(gro1vect, gro2vect) / (gro1vectnorm * gro2vectnorm)
        angle[x] = np.arccos(cost)
    maxangle = angle
    return maxangle


# MAIN BODY
def tossconfigs(inputname, outputname, restraint):
    file = open(inputname, 'r')
    lines = file.readlines()
    #lines = filter(None, (line.rstrip() for line in file))
    file.close()

    nam = 12
    nm = 72
    # Read in the total number of atoms in the system (should be the only item on the second line of the .gro file)
    print(lines[1])
    if len(lines[1].split()) == 1:
        # na is the number of total atoms
        na = int(lines[1].split()[0])
        # nm is the number of total molecules
        nm = na / nam
    else:
        sys.exit('Unexpected .gro file format')

    coord, boxv, ngro = readgro(inputname, na, nm, nam)
    coord2, boxv, ngro2 = readgro(restraint, na, nm, nam)
    cshape = np.shape(coord)
    icoord = np.zeros(cshape)

    for x in range(0, ngro):
        icoord[x, :, :, :] = coord2[0, :, :, :]

    #todel1= []    
    todel2 = []

    allcent, idealcent = find_centroids(coord, ngro, boxv, nm)

    for x in range(ngro):
        for mol in range(nam):
             coord[x, :, :, mol] = coord[x, :, :, mol] - allcent[x, :, :, 0]
             icoord[x, :, :, mol] = icoord[x, :, :, mol] - allcent[0, :, :, 0]

    maxangles = np.zeros([ngro, 4])
    for x in range(ngro):
        maxangles[x, :] = normalcompare(coord[x, :, :, :], icoord[x, :, :, :]) * (180 / np.pi)

    """
    fig=plt.figure(0)
    plt.hist(maxangles[:ngro,0],40,alpha=0.3)
    plt.show()
    pdb.set_trace()	
    plt.hist(maxangles[:ngro,1],40,alpha=0.3)
    plt.show()
    pdb.set_trace()
    plt.hist(maxangles[:ngro,2],40,alpha=0.3)
    plt.show()
    pdb.set_trace()
    plt.hist(maxangles[:ngro,3],40,alpha=0.3)
    plt.show()
    """

    for x in range(ngro):
        for y in range(4):
            if maxangles[x, y] < 135 and maxangles[x, y] > 45:
                if x not in todel2:
                    todel2 = np.append(todel2, x)

    for x in range(ngro):
        for mol in range(nam):
            coord[x, :, :, mol] = coord[x, :, :, mol] + allcent[x, :, :, 0]
            icoord[x, :, :, mol] = icoord[x, :, :, mol] + allcent[0, :, :, 0]

    newgro = np.delete(coord, todel2, 0)

    frames = np.shape(newgro)[0]
    writegrofile(outputname, lines, frames, nm, nam, newgro, boxv)
    print(todel2)
    return todel2
