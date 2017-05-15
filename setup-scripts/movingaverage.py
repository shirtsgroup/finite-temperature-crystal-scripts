#/usr/bin/python
#PYTHON SCRIPT TO TAKE IN AN XVG FILE AND CREATE A MOVING AVERAGE OF THE SAME DATA
import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-f', dest = 'infile', help = 'Input xvg file')
(options, args) = parser.parse_args()

ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================', '=====================', 'TIME(ps)'];
fname = str(options.infile)
infile = open(fname,'r')
lines = infile.readlines()
infile.close()
print "loading " + str(fname)
fname = fname[0:len(fname)-4]
fname = fname + '_avg.xvg'
outfile = open(fname, 'w')

counter = 0
movingavg = numpy.zeros(20,float)
for line in lines:
    tokens = line.split()
    if tokens[0] in ignore_symbols:
	outfile.write(line)
	continue
    for i,token in enumerate(tokens):
	if i==0:
	    continue
        movingavg[i] = movingavg[i]*float(counter)/(counter+1)+1.0/(counter+1)*float(token)
    outfile.write(str(tokens[0]) + ' '.join([str(j) for j in movingavg[:len(tokens)]]) + '\n')
    counter+=1
