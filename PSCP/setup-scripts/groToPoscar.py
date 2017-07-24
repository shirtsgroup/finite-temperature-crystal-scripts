# -*- coding: utf-8 -*-

#!/usr/bin/python
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import sys
import optparse

parser = optparse.OptionParser()

parser.add_option("-I", "--inputgro", dest="filename",default='formam_p1_4.gro',
                  help = "initial gro file")
parser.add_option("-D", "--decscription", dest="descrip",default='benzene_p1',
                  help = "file description for POSCAR")
                  
(options, args) = parser.parse_args()

filename = options.filename
descrip = options.descrip

def readgro(filename):
    file = open(filename, 'r')
    lines = file.readlines() 
    nummol = int(lines[1])
    coords = np.zeros((nummol,3))
    identity = []
    types = []
    for x in range(2,2+nummol):
        line = lines[x].split()
        coords[x-2,0] = float(line[3])
        coords[x-2,1] = float(line[4])
        coords[x-2,2] = float(line[5])
        identity.append(str(line[1]))
        if str(line[1]) not in types:
            types.append(str(line[1]))
            
    amounts = np.zeros((len(types),), dtype = int)
    for x in range(len(types)):
        for y in range(len(identity)):
            if identity[y] == types[x]:
                amounts[x] = amounts[x]+1   
    boxvect = lines[nummol+2].split()
    boxv = [float(boxvect[0]),float(boxvect[1]),float(boxvect[2])]
            
    return coords, identity, types, amounts, boxv
    
def writeposcar(name, coords, identity, types, amounts, boxv):   
    outgrofile = open('POSCAR','w')
    outgrofile.write(name+'\n')
    outgrofile.write('1.0'+'\n')
    outgrofile.write(str(boxv[0])+' 0.0 0.0'+'\n')
    outgrofile.write('0.0 '+str(boxv[1])+' 0.0'+ '\n')
    outgrofile.write('0.0 0.0 '+str(boxv[1])+ '\n')
    boxiline=''
    boxtline=''
    for d in range(len(types)):
        iline = (types[d]+ ' ')
        tline = (str(amounts[d])+' ')
        boxiline += iline
        boxtline += tline
    outgrofile.write(boxiline+'\n') 
    outgrofile.write(boxtline+'\n')
    
    for d in range(len(identity)):
        line = str(coords[d,0])+' '+str(coords[d,1])+' '+str(coords[d,2])+' '+identity[d]+'\n'
        outgrofile.write(line)

def writeespresso(name, identity, coords, boxv):
    outespfile = open('FORMtest', 'w')
    for x in range(len(identity)):
        line = identity[x]+' '+str(coords[x,0]/boxv[0])+' '+str(coords[x,1]/boxv[1])+' '+str(coords[x,2]/boxv[2])+'\n'
        outespfile.write(line)
        
    
coords, identity, types, amounts, boxv = readgro(filename)
print(coords)
print(boxv)
writeespresso('BNZtest', identity, coords, boxv)
#writeposcar(descrip, coords, identity, types, amounts, boxv)
