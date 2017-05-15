#/usr/bin/python
#PYTHON SCRIPT TO TAKE IN AN XVG FILE AND FIND THE CONFIGURATION THAT BEST MATCHES THE AVERAGE VALUE
import numpy # numerical array library
from optparse import OptionParser # for parsing command-line options
import subprocess
import os
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-f', dest = 'infile', help = 'Input xvg file', default = 'volume.xvg')
parser.add_option('-m', dest = 'mdpfile', help = 'mdp file used to run the simulation', default = 'benzene_production.mdp')
parser.add_option('-i', dest = 'ignoreframes', help = 'Equilibration frames to ignore', default = '100')
parser.add_option('-n', dest = 'column', help = 'Column of the data in the file', default = '2')
parser.add_option('-V', dest = 'Volume', help = 'Optional volume to match', default = '-1')

GRO_LOC="/home/common/gromacs-forceaverage-3/install/bin"

(options, args) = parser.parse_args()
ignoreframes = int(options.ignoreframes)
column = int(options.column)
volume = float(options.Volume)
ignore_symbols = ['#', '@', '@TYPE', 'STEP', '====================='];
fname = str(options.infile)
infile = open(fname,'r')
lines = infile.readlines()
infile.close()
print "loading " + str(fname)

if volume == -1:
    #Loop through and calculate the average
    ignorecounter = 0
    counter = 0
    average = 0
    for line in lines:
	tokens = line.split()
	if tokens[0] in ignore_symbols:
	    continue

	if ignorecounter < ignoreframes:
	    ignorecounter += 1
	    continue

	average = float(average)*float(counter)/(counter+1)+1.0/(counter+1)*float(tokens[column-1])
	counter+=1
else:
    average=volume


#Now determine which frames were actually saved

commandstring = "less " + options.mdpfile + " | grep nstlog | grep -v compressed | awk '{print $3}'"
edr_steps = int(subprocess.Popen(commandstring, shell=True, stdout=subprocess.PIPE).communicate()[0])
commandstring = "less " + options.mdpfile + " | grep nstxout | grep -v compressed | awk '{print $3}'"
trj_steps = int(subprocess.Popen(commandstring, shell=True, stdout=subprocess.PIPE).communicate()[0])
frequency = int(trj_steps/edr_steps)
print edr_steps
print trj_steps

print frequency

#Now loop through again and find the index of the closest configuration to the average
ignorecounter = 0
counter = 0
closest_index = int(0)
closest_timestep = float(0.0)
closest_value = float(-1234.0)
for line in lines:
    tokens = line.split()
    if tokens[0] in ignore_symbols:
        continue

    if ignorecounter < ignoreframes:
        ignorecounter += 1
	counter+=1
        continue

    if counter % frequency > 0:
	counter+=1
	continue
    print float(tokens[0])
    if abs(float(tokens[column-1]) - float(average)) < abs(closest_value - average):
	closest_index = int(counter)
	closest_timestep = float(tokens[0])+0.00001
	closest_value = float(tokens[column-1])
    counter+=1

#Create the new gro file
commandstring = "echo 0 | " + GRO_LOC + "/trjconv_d -f benzene_PROD.trr -s benzene_PROD.tpr -o avgconfig.gro -pbc whole -ndec 8 -dump " + str(closest_timestep)
os.system(commandstring)

#Print out the results
print "The average is: "
print str(average)
print "The closest value in the trajectory is: "
print str(closest_value)
print "The deviation is: "
print str(abs(closest_value - average))
print "The index is: "
print str(closest_index)
print "The timestep is: "
print str(closest_timestep)
