#/usr/bin/python
#PYTHON SCRIPT TO TAKE IN AN XVG FILE AND CREATE A HISTOGRAM OF DATA
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import pdb
import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties as FP

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-f', dest = 'infiles', help = 'Input xvg file')
parser.add_option('-i', dest = 'ignoreframes', help = 'Equilibration frames to ignore', default = '200')
parser.add_option('-c', dest = 'columns', help = 'Column of the data in the file', default = '2')
parser.add_option('-l', dest = 'labels', help = 'Label of each dataset', default = 'Set1') 
parser.add_option('-b', dest = 'bins', help = 'Number of bins to use in the histogram', default = '10')
(options, args) = parser.parse_args()

#Read in the list of files, columns, and labels
infiles=[]
columns=[]
labels=[]
for i,token in enumerate(options.infiles.split()):
    infiles.append(token)

for i,token in enumerate(options.columns.split()):
    columns.append(token)

for i,token in enumerate(options.labels.split()):
    labels.append(token)

ignoreframes = int(options.ignoreframes)
bins = int(options.bins)
ignore_symbols = ['#', '@', '@TYPE', 'STEP', '====================='];
data_vector = numpy.zeros([len(infiles),100000],float)
N_k = numpy.zeros(len(infiles),int)
Colors=['b','g','r', 'c', 'm'];

#Read in the data
for i,file in enumerate(infiles):
    infile = open(file,'r')
    lines = infile.readlines()
    infile.close()
    print "loading " + str(file)

    ignorecounter = 0
    counter = 0
    for line in lines:
	tokens = line.split()
	if tokens[0] in ignore_symbols:
	    continue
	if ignorecounter < ignoreframes:
	    ignorecounter += 1
	    continue

	data_vector[i,counter] = float(tokens[int(columns[i])-1])
	counter+=1
    N_k[i]=counter

#Now histogram the data
title="Average: "
for i in range(len(infiles)):
    plt.hist(data_vector[i,0:N_k[i]],bins,alpha=0.3,facecolor=Colors[i],label=labels[i])
    title=title+str(str(numpy.average(data_vector[i,0:N_k[i]])) + " ")
#plt.hist(data_vector[0,0:N_k[0]]-data_vector[1,0:N_k[1]]-data_vector[2,0:N_k[2]],bins,alpha=0.3,facecolor=Colors[3],label="CHOnly")
plt.title(title)
if options.labels != 'Set1':
    plt.legend(loc='upper right')
plt.show()
