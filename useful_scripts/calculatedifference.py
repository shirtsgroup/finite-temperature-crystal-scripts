#/usr/bin/python
#PYTHON SCRIPT TO TAKE IN TWO FILES AND CALCULATE THE DIFFERENCE AT EACH STEP
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-f', dest = 'infile1', help = 'Input file 1', default="DefaultFile1")
parser.add_option('-g', dest = 'infile2', help = 'Input file 2', default="DefaultFile2")
parser.add_option('-o', dest = 'outfile', help = 'Output file', default="DefaultOutput")
(options, args) = parser.parse_args()

if options.infile1 == "DefaultFile1" or options.infile2 == "DefaultFile2" or options.outfile == "DefaultOutput":
    print "Inputs Wrong!"
    sys.exit()

ignore_symbols = ['#', '@', '@TYPE', 'STEP', '====================='];
fname = str(options.infile1)
infile = open(fname,'r')
lines1 = infile.readlines()
infile.close()
print "loading " + str(fname)

fname = str(options.infile2)
infile = open(fname,'r')
lines2 = infile.readlines()
infile.close()
print "loading " + str(fname)

fname = str(options.outfile)
outfile = open(fname, 'w')

file2_linenum = 1
for line in lines1:
    tokens1 = line.split()
    if tokens1[0] in ignore_symbols:
	continue
    tokens2 = lines2[file2_linenum].split()
    while tokens2[0] in ignore_symbols:
	file2_linenum+=1
	tokens2 = lines2[file2_linenum].split()	
    outfile.write(str(tokens1[0]) + ' ' + str(float(tokens1[1])-float(tokens2[1])) + '\n')
    file2_linenum+=1
outfile.close()
