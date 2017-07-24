# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 09:31:24 2016

@author: nasc5274
"""
import numpy as np
import sys
import optparse

parser = optparse.OptionParser()

parser.add_option("-I", "--inputgro", dest="inname",default='1BNZtsest.pos',
                  help = "initial tinker file")
parser.add_option("-B", "--inputboxv", dest="boxvname",default='traj.xsf',
                  help = "initial boxv file")
parser.add_option("-O", "--decscription", dest="outname",default='ben_p1_CPMD.gro',
                  help = "gro file output name")
parser.add_option("-F", "--frames", dest="frames",default=30,
                  help = "number of frames")
parser.add_option("-P", "--moleculesize", dest="mpa",default=12,
                  help = "atoms per molecule")
parser.add_option("-M", "--molecules", dest="mole",default=4,
                  help = "number of molecules")
parser.add_option("-T", "--type", dest="type",default='BNZ',
                  help = "type of molecule")
parser.add_option("-D", "--workdir", dest="workdir",default='/home/nasc5274/finishedjobs/QEtesting/QEtesting/',
                  help = "initial boxv file")
                  
(options, args) = parser.parse_args()

filename = options.inname
outname = options.outname
framenum = options.frames
mpa = options.mpa
mole = options.mole
typemol = options.type
boxvfile = options.boxvname
workdir = options.workdir

filename = workdir+filename
outname = workdir+outname
boxvfile = workdir+boxvfile

numtype1 = 6
numtype2 = 6
numtype3 = 0
numtype4 = 0

def readQE(filename, boxvfile, frames, atpermol, molnum):
    file = open(filename, 'r')
    fileb = open(boxvfile, 'r')
    linesb = fileb.readlines()
    lines = file.readlines()
    boxvtemp = np.zeros((3,))
    boxvtemp[0] = linesb[2].split()[0]    
    boxvtemp[1] = linesb[3].split()[1] 
    boxvtemp[2] = linesb[4].split()[2]
    
    numatom = atpermol*molnum
    coords = np.zeros((numatom,3,frames))
    boxv = np.zeros((frames,3))
    linetype = []
    it = 0   

    for x in range(frames):
        for g in range(0,numtype1):            
            linetype.append('H')
        for g in range(0,numtype2):            
            linetype.append('C')
        for g in range(0,numtype3):            
            linetype.append('N')      
        for g in range(0,numtype4):            
            linetype.append('O')
        boxv[x,:] = boxvtemp[0:3]

        for y in range(numatom):
            it = it+1
            coordtemp = lines[it].split()
            coords[y,:,x] = coordtemp[:]

        it = it+1
    return boxv,coords, linetype
        
           
    

def writegrofile(output,acoordBprime, frames, na, nam, nm, molname, boxvB, lines):
    outgrofile = open(output,'w')
    x = np.zeros(3)
    for gro in range(int(frames)):

        outgrofile.write('Converted from Tinker to GRO ' + str(gro) + '\n')    
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
                    x[d] = round(acoordBprime[tacount-1,d,gro],8)
                    xstr.append("%13.8f" % (x[d]))
                line = str(mcount).rjust(5) + molname + lines[tacount+1].rjust(7) + str(tacount).rjust(5) + xstr[0] + xstr[1] + xstr[2] + '\n'

                outgrofile.write(line)
        boxvline = ''    
        for d in range(3):

	    bline = ("%13.8f" % (round(boxvB[gro,d],8)))
            boxvline += bline
        outgrofile.write(boxvline+'\n')      # make sure the .gro ends on a newline
    outgrofile.close()


boxv,coords, lines = readQE(filename, boxvfile,framenum,mpa, mole)

writegrofile(outname,coords, framenum,mpa*mole, mpa, mole, typemol, boxv, lines)
