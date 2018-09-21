
#
# Conducts bootstrap sampling on a u_kln Matrix using MBAR
# 
# Copyright Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import pdb
import random
import os.path
import sys

#=============================================================================================
# INPUTS FROM THE USER
#=============================================================================================
#u_kln 		- The u_kln matrix that houses the raw data and that will be processed with MBAR
#N_k 		- The N_k matrix that determine the number of data points in each state (Either 1xK or NxK)
#beta_k		- Vector of beta values for each state
#N		- Number of molecules in the system
#indexVect	- 2xN Vector of the state free energy difference to store after the MBAR calculation
#datafile	- The name of the file which will store the bootstrap data
#stdfile	- The name of the file which will store the standard devation of the data
#n		- The total number of bootstrap iterations to conduct
#Header		- Optional string header to put at the top of the datafiles
#Entropy	- Optional flag to track the average energy and entropy

def runMBARBootstrap(u_kln_orig, N_k, beta_k, N, indexVect, datafile, stdfile, n, Header="None", Entropy='No'):
    i=0;
    j=0;
    k=0;
    previous_lines=0
    dG_Gas=dict()
    dG_Solv=dict()
    dG_Gas_all = 0
    dG_Solv_all = 0
    iteration_max = n
    u_kln=u_kln_orig.copy()
    K = len(u_kln[:,0,0])
    num_indicies = len(indexVect[0,:])
    results = numpy.zeros([num_indicies, iteration_max],float) #The free energy differences at each bootstrap iteration
    uncertainties = numpy.zeros([num_indicies, iteration_max], float) #The rolling standard deviation
    resultSlice = numpy.zeros(num_indicies+1,float)
    uncertaintySlice = numpy.zeros(num_indicies+1,float)
    if Entropy=='yes':
	results_U = numpy.zeros([num_indicies, iteration_max],float) #The average energy difference at each bootstrap iteration
    	uncertainties_U = numpy.zeros([num_indicies, iteration_max], float) #The rolling standard deviation
	resultSlice_U = numpy.zeros(num_indicies+1,float)
	uncertaintySlice_U = numpy.zeros(num_indicies+1,float)
	results_S = numpy.zeros([num_indicies, iteration_max],float) #The entropy differences at each bootstrap iteration
        uncertainties_S = numpy.zeros([num_indicies, iteration_max], float) #The rolling standard deviation
	resultSlice_S = numpy.zeros(num_indicies+1,float)
        uncertaintySlice_S = numpy.zeros(num_indicies+1,float)
	
    #=============================================================================================
    # ENSURE THE INPUTS ARE CORRECT
    #=============================================================================================
    if len(numpy.shape(N_k))==1 and len(N_k) != len(u_kln[:,0,0]):
	print "Number of states in the U_kln matrix: " + str(len(u_kln[:,0,0])) + " and in the N_k matrix: " + str(len(N_k)) + " are not the same!"
	sys.exit()
    elif int(len(indexVect[:,0]) != int(2)):
	print "indexVect shape: " + str(numpy.shape(indexVect)) + " is not valid"
	sys.exit()

    if len(numpy.shape(N_k)) != 1 and len(N_k[:,1]) != len(indexVect[0,:]):
	print "N_k shape: " + str(numpy.shape(N_k)) + " and indexVect shape: " + str(numpy.shape(indexVect)) + " is not valid"


    #Read in any previously generated data in the data file
    fnameData = datafile + '.txt'
    j=0
    if os.path.isfile(fnameData):
	infile = open(fnameData, 'r')
	lines = infile.readlines()
	infile.close()
	print "loading " + fnameData
	for line in lines:
	    token=line.split()
	    if token[0] == '#':
		continue
	    for i in range(num_indicies):
		results[i,j] = token[i+1]
	    j=j+1

    #Read in any previously generated data in the std file
    fnameStd = stdfile + '.txt'
    j=0
    if os.path.isfile(fnameStd):
	infile = open(fnameStd, 'r')
	lines = infile.readlines()
	infile.close()
	print "loading " + fnameStd
	for line in lines:
	    token=line.split()
	    if token[0] == '#':
		continue
	    for i in range(num_indicies):
		uncertainties[i,j] = token[i+1]
	    j=j+1
    previous_lines=j

    # Constants.
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184) # Boltzmann constant in kcal/mol

    # all data loaded from the three sets
    u_kln_save = u_kln_orig.copy() 

    #Now run the iterative bootstrapping algorithm on this dataset
    for j in range (previous_lines,iteration_max):
	print "BOOTSTRAP ITERATION: " + str(j+1)

	if len(numpy.shape(N_k))==1:
	    N_ks=1
	else:
	    N_ks=len(indexVect[0,:])

	for s in range(N_ks):

	    #Grab the new N_k matrix
	    if N_ks > 1:
		N_k_vect=N_k[s,:]
	    else:
		N_k_vect=N_k

	    #Generate the new random u_kln matrix
	    u_kln.fill(0) # u_kln[k,l,n] is the reduced potential energy of snapshot n from state k in state l
	    for k in range(K):
		for n in range(N_k_vect[k]):
		    randomnum = int(random.random()*N_k_vect[k])
		    #randomnum= n #FIX FIX FIX FIX FIX
		    u_kln[k,:,n] = u_kln_save[k,:,randomnum]
	
	    # Initialize MBAR.
	    print "Running MBAR..."

	    # generate the weights
	    if j>previous_lines:
		mbar = pymbar.MBAR(u_kln, N_k_vect, verbose = True, initial_f_k = guess_f_k, method = 'adaptive', use_optimized=False)
	    else:
		mbar = pymbar.MBAR(u_kln, N_k_vect, verbose = True, method = 'adaptive', use_optimized=False)	
	    (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()

	    # extract self-consistent weights and uncertainties
	    guess_f_k = mbar.f_k

	    #convert to kcal/mol
	    for k in range(K):
		df_i[k] /= (beta_k[k]*float(N))
		ddf_i[k] /= (beta_k[k]*float(N))
	
	    if N_ks > 1:
		resultSlice[s+1]=(df_i[indexVect[1,s],indexVect[0,s]])

	#Now extract the final values and store them in the output file
	
	#If this is the first iteration, pretty print a title
	if j==0:
	    with open(fnameData, "a") as myfile:
		myfile.write('# MBAR Bootstrap Data')
		myfile.write('\n')
		if Header != "None":
		    myfile.write(Header)
		    myfile.write('\n')

	#Store the data in the results file
	resultSlice[0]=j+1
	if N_ks==1:
	    for i in range(num_indicies):
	        resultSlice[i+1]=(df_i[indexVect[1,i],indexVect[0,i]])
	with open(fnameData, "a") as myfile:
	    myfile.write(' '.join(str(cell) for cell in resultSlice))
	    myfile.write('\n')
	results[:,j]=resultSlice[1:]

	##Compute the average energy difference if the entropy flag is specified
        #    if Entropy=='yes':
        #        for i in range(num_indicies):
        #            resultSlice[i] = numpy.average(u_kln[indexVect[1,i],indexVect[1,i],:]) - numpy.average(u_kln[indexVect[0,i],indexVect[0,i],:])




	#If this is the first iteration, pretty print a title
	if j==0:
	    with open(fnameStd, "a") as myfile:
		myfile.write('# MBAR Bootstrap Uncertainty')
		myfile.write('\n')
		if Header != "None":
                    myfile.write(Header)
                    myfile.write('\n')

	
	#Calculate the variance and store it in the std file
	uncertaintySlice[0]=j+1
	for i in range(num_indicies):
	    uncertainties[i,j]=numpy.std(results[i,:(j+1)])
	uncertaintySlice[1:]=uncertainties[:,j]
	with open(fnameStd, "a") as myfile:
	    myfile.write(' '.join(str(cell) for cell in uncertaintySlice))
	    myfile.write('\n')

    #Finally, calculate the average bootstrap value and report it at the end of the Data file
    Average=numpy.zeros(len(results[:,0]))
    for i in range(num_indicies):
	Average[i]=numpy.average(results[i,:])
    with open(fnameData, "a") as myfile:
	myfile.write('Avg:  ' + ' '.join(str(cell) for cell in Average))
	myfile.write('\n')


#fname			- The name of a bootstrapped std file to read in
#@Return		- The vector of bootstrapped uncertainties

def ExtractBootstrap(fname):

    if not os.path.isfile(fname):
        print "File: " + fname + " not found!"
	os.exit()
    
    infile = open(fname, 'r')
    lines = infile.readlines()
    infile.close()

    tokens = lines[len(lines)-1].split()
    ddG=numpy.zeros(len(tokens)-1,float)

    for i,token in enumerate(tokens):
	if i==0:
	    continue
	else:
	    ddG[i-1]=float(token)

    return ddG


#fname                  - The name of a bootstrapped data file to read in
#linenum		- The line in the data file to pull results from

#@Return                - The vector of bootstrapped data points

def ExtractBootstrapLine(fname, linenum):

    if not os.path.isfile(fname):
        print "File: " + fname + " not found!"
        os.exit()

    infile = open(fname, 'r')
    lines = infile.readlines()
    infile.close()

    for i,line in enumerate(lines):
	tokens = line.split()[0]
	try:	
	    value = int(tokens[0])
	except ValueError:
	    continue

	if value < int(linenum):
	    continue
    
	break

    tokens = line.split()
    dG=numpy.zeros(len(tokens)-1,float)

    for t,token in enumerate(tokens):
        if t==0:
            continue
        else:
            dG[t-1]=float(token)

    return dG
















