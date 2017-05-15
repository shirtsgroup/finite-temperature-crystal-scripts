#!/bin/python

#
# Replicate the unit cell of a polymorph into a supercell and expand its contents if desired
#
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import os
import pdb

parser = OptionParser()
parser.add_option('-g', dest = 'grofile', help = 'Name of the gro file to expand', default = 'None')
parser.add_option('-n', dest = 'polymorph', help = 'Polymorphs to analyze', default = 'p1')
parser.add_option('-M', dest = 'molecule', help = 'Name of the molecule', default = 'benzene')
parser.add_option('-N', dest = 'nmol', help = 'Number of molecules in the unit cell', default = 4)
parser.add_option('-l', dest = 'length', help = 'Supercell length', default = 1)
parser.add_option('-w', dest = 'width', help = 'Supercell width', default = 1)
parser.add_option('-z', dest = 'height', help = 'Supercell height', default = 1)
parser.add_option('-H', dest = 'Hinge', help = 'Optional String Hinge', default = '')
parser.add_option('-o', dest = 'outfile', help = 'Name of the output gro file', default = 'None')
nndx = 15 # number of ndx entries per line
#Ensure that user inputs are sensible
(options, args) = parser.parse_args()
poly = options.polymorph
molecule=options.molecule
nmol=int(options.nmol)
Length = int(options.length)
Width = int(options.width)
Height = int(options.height)
Hinge = options.Hinge
if Hinge != "":
    Hinge="_" + Hinge

#Read in the unit cell
if options.grofile == 'None':
    fname = molecule + "_" + poly +"_" + str(nmol) + Hinge + ".gro"
else:
    fname = options.grofile
infile = open(fname, 'r')
lines = infile.readlines()
infile.close()
print "loading " + fname

#Determine the unit dimensions from the last line of the .gro file
v1 = numpy.zeros(3,float)
v2 = numpy.zeros(3,float)
v3 = numpy.zeros(3,float)
tokens = lines[len(lines)-1].split()
stopper = ""	#Need to use this later during the loop
if len(tokens)==3:
    stopper = tokens[0]
    v1[0] = float(tokens[0])
    v2[1] = float(tokens[1])
    v3[2] = float(tokens[2])
    boxvect_length=3
else:
    stopper = tokens[0]
    v1[0] = float(tokens[0])
    v2[1] = float(tokens[1])
    v3[2] = float(tokens[2])
    v1[1] = float(tokens[3])
    v1[2] = float(tokens[4])
    v2[0] = float(tokens[5])
    v2[2] = float(tokens[6])
    v3[0] = float(tokens[7])
    v3[1] = float(tokens[8])
    boxvect_length=9

#Determine the number of molecules in the unit cell
tokens = lines[1].split()
#molecule_number = int(tokens[0])/apermol
apermol=int(tokens[0])/nmol

#Open the output gro file
if options.outfile == "None":
    numout = nmol*Length*Width*Height
    fname = molecule + "_" + poly + "_" + str(numout) + Hinge + ".gro"
    fname_ind = molecule + "_" + poly + "_" + str(numout) + "_" + str(nmol) + "ind" + Hinge + ".gro"
else:
    numout = nmol*Length*Width*Height
    fname = options.outfile
    fname_ind = options.outfile
open(fname, 'w').close() 	#delete any previous contents
outfile = open(fname, 'w')

#Open the output index file
fname_index = molecule + "_" + poly + "_" + str(numout) + "_" + str(nmol) + "ind" + ".ndx"
open(fname_index, 'w').close()        #delete any previous contents
outfile_index = open(fname_index, 'w')

#Open the output symmetry groups file
fname_symmetry = molecule + "_" + str(nmol) + "ind_symmetry_groups.txt"
open(fname_symmetry, 'w').close()        #delete any previous contents
outfile_symmetry = open(fname_symmetry, 'w')

#Write the first two lines to the new gro file
outfile.write(lines[0])
tokens = lines[1].split()
tokens[0] = str(int(tokens[0])*Length*Width*Height)		
outfile.write(" " + tokens[0])
outfile.write('\n')

#Initialize the atom header container and the dictionary of atom numbers for each header
index_headers = [];
atom_indicies = dict();
atom_types = dict();

atom_types['C']=0
atom_types['H']=0

atom_count = 0
total_atom_count = 0
linecount = 0
moleculecount = 1
for i in range(0,Length):
    for k in range(0,Height):
	for j in range(0,Width):
	    for l,line in enumerate(lines):
		tokens = line.split()
		linecount+=1
		if tokens[0] == "benzene" or tokens[0] == stopper:
        	    continue
    		if l == 0 or l == 1:
		    continue

		if (float(tokens[3])+v1[0]*i+v2[0]*j+v3[0]*k) < 0.0:
		    tokens[3] = "  " + str("%.8f" % (float(tokens[3])+v1[0]*i+v2[0]*j+v3[0]*k))
		else:
		    tokens[3] = "   " + str("%.8f" % (float(tokens[3])+v1[0]*i+v2[0]*j+v3[0]*k))
		if (float(tokens[4])+v1[1]*i+v2[1]*j+v3[1]*k) < 0.0:
                    tokens[4] = "  " + str("%.8f" % (float(tokens[4])+v1[1]*i+v2[1]*j+v3[1]*k))
                else:
                    tokens[4] = "   " + str("%.8f" % (float(tokens[4])+v1[1]*i+v2[1]*j+v3[1]*k))
		if (float(tokens[5])+v1[2]*i+v2[2]*j+v3[2]*k) < 0.0:
                    tokens[5] = "  " + str("%.8f" % (float(tokens[5])+v1[2]*i+v2[2]*j+v3[2]*k))
                else:
                    tokens[5] = "   " + str("%.8f" % (float(tokens[5])+v1[2]*i+v2[2]*j+v3[2]*k))

		atom_count+=1
		total_atom_count+=1 
		if atom_count > apermol:	#We have reached another molecule
		    atom_count = 1;
		    moleculecount+=1

		#Now format the outputs
		
		#If this is the first time encountering this atom, create a header for the index file
                if i==0 and j==0 and k==0:
		    if tokens[1] not in atom_types:
			atom_types[tokens[1]] = 0
                    #header = '[ ' + tokens[1] + '    ' + str((int(atom_types[tokens[1]]) + 1)) + ' ]'
                    header = '[ ' + tokens[1] + str((int(atom_types[tokens[1]]) + 1)) + ' ]'
                    atom_types[tokens[1]] = atom_types[tokens[1]] + 1
                    index_headers.append(header)
                    atom_indicies[header] = str(total_atom_count)
                else:
                    atom_indicies[index_headers[l-2]] = atom_indicies[index_headers[l-2]] + '    ' + str(total_atom_count)

		if moleculecount < 10:
		    tokens[0] = "    " + str(moleculecount) + "BNZ"
		elif moleculecount < 100:
		    tokens[0] = "   " + str(moleculecount) + "BNZ"
		elif moleculecount < 1000:
                    tokens[0] = "  " + str(moleculecount) + "BNZ"
		else:
		    tokens[0] = " " + str(moleculecount) + "BNZ"
	
		while len(tokens[1]) < 7:	
		    tokens[1] = " " + tokens[1]

		if total_atom_count < 10:
                    tokens[2] = "    " + str(total_atom_count)
                elif total_atom_count < 100:
                    tokens[2] = "   " + str(total_atom_count)
                elif total_atom_count < 1000:
                    tokens[2] = "  " + str(total_atom_count)
                else:
                    tokens[2] = " " + str(total_atom_count)
		outfile.write(''.join(tokens))
		outfile.write('\n')

#write the last line to the output file
tokens = lines[len(lines)-1].split()
if len(tokens)==3:
    tokens[0] = str(v1[0]*Length)
    tokens[1] = str(v2[1]*Width)
    tokens[2] = str(v3[2]*Height)
else:
    tokens[0] = str(v1[0]*Length)
    tokens[1] = str(v2[1]*Width)
    tokens[2] = str(v3[2]*Height)
    tokens[3] = str(v1[1]*Length)
    tokens[4] = str(v1[2]*Length)
    tokens[5] = str(v2[0]*Width)
    tokens[6] = str(v2[2]*Width)
    tokens[7] = str(v3[0]*Height)
    tokens[8] = str(v3[1]*Height)
outfile.write("   " + '   '.join(tokens))
outfile.write('\n')
outfile.close()
print 'wrote final gro file to: ' + fname
cpcommand = 'cp ' + fname + ' ' + fname_ind
os.system(cpcommand)
print 'also wrote final gro file to: ' + fname_ind

os.system('python ~/crystals/NMA/scripts/removeCOM.py -f ' + fname)
os.system('python ~/crystals/NMA/scripts/removeCOM.py -f ' + fname_ind)


#Write to the index file
# write [ system ] group
outfile_index.write('[ System ]\n')
for i in range(1,numout*apermol+1):
    outfile_index.write("%5d" % (i))
    if (i % nndx) == 0:
        outfile_index.write('\n')
outfile_index.write('\n')
for header in index_headers:
    outfile_index.write(header)
    outfile_index.write('\n')
    outfile_index.write(atom_indicies[header])
    outfile_index.write('\n')
outfile_index.close()
print 'wrote final ndx file to: ' + fname_index
print 'put this section in the mdp:'
print '-------------------'
print 'symmetry-averaging = yes'
outfile_symmetry.write('ave-ngroups = %d' % (len(index_headers)))
outfile_symmetry.write('\n')
for i,header in enumerate(index_headers):
    grpname = header[2:-2]
    outfile_symmetry.write("ave-group%d = %s" % (i,grpname))
    outfile_symmetry.write('\n')
outfile_symmetry.close()

#Write the expansion file
fname_expansion = molecule + '_' + poly + '_' + str(nmol) + '_expansion.txt'
os.system('cp expansion.txt '  + fname_expansion)
sedcommand = "sed -i s/XXXX/" + str(Length) + "/g " + fname_expansion
os.system(sedcommand)
sedcommand = "sed -i s/YYYY/" + str(Width) + "/g " + fname_expansion
os.system(sedcommand)
sedcommand = "sed -i s/ZZZZ/" + str(Height) + "/g " + fname_expansion
os.system(sedcommand)
#print 'ave-ngroups = %d' % (len(index_headers))
#for i,header in enumerate(index_headers):
#    grpname = header[2:-2]
#    print "ave-group%d = %s" % (i,grpname)
#print '-------------------'
