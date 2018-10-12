from __future__ import print_function
#
# Computing the free energy difference of an organic crystal polymorph at different temperatures using
# the constant volume variant of the Gibbs-Helmholtz Equation
# 
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-p', '--plot', dest = 'plot', help = 'Plot output (default false)', default=True, action = 'store_true')
parser.add_option('-C', dest = 'MinTemp', help = 'Minimum temperature to examine', default = 200)
parser.add_option('-H', dest = 'MaxTemp', help = 'Maximum temperature to examine', default = 400)
parser.add_option('-T', dest = 'Tspacing', help = 'Spacing between temperatures', default = 20)
parser.add_option('-G', dest = 'Gamma', help = 'Gamma point during the temperature change', default = 100)
parser.add_option('-n', dest = 'polymorphs', help = 'Polymorphs to analyze', default = 'all')
parser.add_option('-N', dest = 'molecules', help = 'number of supercell molecules', default = 72)
parser.add_option('-L', dest = 'lam', help = 'LAMBDA Point', default = 100)
parser.add_option('-d', '--hinge', dest = 'hinge', help = 'Optional string at end of jobs', default = 'DefaultHinge')

(options, args) = parser.parse_args()
MinT = float(options.MinTemp)
MaxT = float(options.MaxTemp)
dT = float(options.Tspacing)
GAMMA=int(options.Gamma)
Molecules = options.molecules
LAMBDA = options.lam
hinge = options.hinge

if (options.plot):
    import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.font_manager import FontProperties as FP

#=============================================================================================
# ENSURE THAT USER INPUTS ARE SENSIBLE
#=============================================================================================
#TEMPERATURE
if (MinT == -1 ) and (MaxT == -1) and (dT == -1):
    print("Using default values!")
    T_k = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] #The temperatures sampled
    temperature_name = ['100K', '200K', '300K', '400K', '500K', '600K', '700K', '800K', '900K', '1000K']
elif MinT < 0 or MaxT < 0 or dT < 0 or MinT > MaxT:
    print("Invalid Temperature Specifications")
    sys.exit()
else:
    temp = MinT
    i=0
    T_k = numpy.ones(numpy.ceil((MaxT-MinT)/dT)+1,float)
    temperature_name=[];
    while temp < MaxT:
        T_k[i]=float(temp)
	temperature_name.append(str(int(temp))+'K')
	temp=temp+dT
	i=i+1
    #Catch the final temperature
    T_k[len(T_k)-1]=MaxT
    temperature_name.append(str(int(MaxT))+'K')

#LAMBDA
if LAMBDA < 0 or LAMBDA > 100:
    print("Invalid Lambda Point: " + str(LAMBDA))
    sys.exit()

#GAMMA
if GAMMA < 0 or GAMMA > 100:
    print("Invalid Gamma Point: " + str(GAMMA))
    sys.exit()

#POLYMORPH
if (options.polymorphs == 'all'):
    polymorph = ['Polymorph1', 'Polymorph2', 'Polymorph3']
    polymorph_short = ['p1', 'p2', 'p3']
elif (options.polymorphs == 'p1'):
    polymorph = ['Polymorph1']
    polymorph_short = ['p1']
elif (options.polymorphs == 'p2'):
    polymorph = ['Polymorph2']
    polymorph_short = ['p2']
elif (options.polymorphs == 'p3'):
    polymorph = ['Polymorph3']
    polymorph_short = ['p3']
else:
    print("Polymorph Inputs Wrong")
    sys.exit()

#=============================================================================================
# FORMAT INPUTS
#=============================================================================================
# LAMBDA
if LAMBDA < 10:
    Lname="00"+str(LAMBDA)+"L"
elif LAMBDA < 100:
    Lname="0"+str(LAMBDA)+"L"
else:
    Lname="100L"

#GAMMA POINT
Gname = ""
if GAMMA < 10:
    Gname="00" + str(int(GAMMA))+"G"
elif GAMMA < 100:
    Gname="0" + str(int(GAMMA))+"G"
else:
    Gname=str(int(GAMMA))+"G"

# OPTIONAL HINGE
if hinge == "DefaultHinge":
    hinge = ""
else:
    hinge = "_" + hinge


#=============================================================================================
# READ IN RAW DATA
#=============================================================================================
# Constants.
kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184) # Boltzmann constant in kcal/mol

molecule = 'Benzene' #The crystalline molecule
ignoreframes = 500 #Ignore the first set of frames where the system is approaching equilibrium
omitT = ['4','5','6','7','8','9','10','11'] #Temperatures to be omitted from the analysis

# Parameters
T_k = T_k*numpy.ones(len(T_k),float) #Convert temperatures to floats 
print(T_k)
print(temperature_name)
K = len(T_k)  # How many states?
Kbig = K+0 # total number of states examined; none are unsampled
N_max = 200000 # maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
dA = numpy.zeros([3,(K+1)/2],float)
ddA = numpy.zeros([3,(K+1)/2],float)
avgenergy = numpy.zeros([3,len(T_k)],float)
#seeds = [201]; #The random seed used (not included at the moment)
g_k = numpy.zeros([len(T_k)],float)
convert_units = (0.2390057)*numpy.ones(len(T_k),float) #Convert all energies to kcal/mol


for i,poly in enumerate(polymorph):
    # Allocate storage for simulation data
    N_k = numpy.zeros([Kbig],numpy.int32) # N_k[k] is the total number of snapshots from alchemical state k

    u_kln = numpy.zeros([K,Kbig,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k 
            
    for k in range(K):
    	n = 0
	linenum = 0
        # cycle through all the input data
        fname = '../finishedJobs/' + poly + '/benzene_' + polymorph_short[i] + '_' + str(Molecules) + '_' + temperature_name[k] + '_' + Lname + '_' + Gname + hinge + '/energy.xvg'
        if k not in omitT:
            infile = open(fname, 'r')
            lines = infile.readlines()
            infile.close()
            print("loading " + fname)
            n = 0
            for line in lines:
		linenum = linenum + 1
                tokens = line.split()
                if (tokens[0]=='#' or tokens[0]=='@' or tokens[0]=='@TYPE'):
                    continue
		if (linenum<ignoreframes):
		    continue
                for l in range(Kbig):
                    # the energy of every configuration from each state evaluated at each state
                    u_kln[k,l,n] = (float(tokens[1]) * beta_k[l])*convert_units[k]  
                n+=1
            N_k[k] = n 

    u_kln_save = u_kln.copy()
    g_k = numpy.zeros([K])
    """        
    for k in range(K):
        #subsample correlated data - for now, use energy from current state
        if k not in omitT:
            g_k[k] = timeseries.statisticalInefficiency(u_kln[k,k,0:N_k[k]]) 
            print "Correlation time for sampled state %d is %10.3f" % (k,g_k[k])
            # subsample the data to get statistically uncorrelated data
            indices = numpy.array(timeseries.subsampleCorrelatedData(u_kln[k, k, 0:N_k[k]], g=g_k[k]))  # subsample
            N_k[k] = len(indices)
            u_kln[k,:,0:N_k[k]] = u_kln_save[k,:,indices].transpose()  # not sure why we have to transpose
    """
    print("Number of retained samples")
    print(N_k)

    # convert to nondimensional units from kcal/mol (already done)

    #u_kln *=  beta_k[0]
    # all data loaded from the three sets 

    #=============================================================================================
    # INTEGRATE POINTS USING SIMPSON INTEGRATION
    #=============================================================================================
    #Now divide by negative NRT^2 to have the matrix in the form of -u/RT^2
    for k in range(K):
	u_kln[k,:,:]=u_kln[k,:,:]*(-1/T_k[k])*(1/beta_k[k])/(float(Molecules))

    #Now run along the diagonal of the matrix and find the average
    for k in range(K):
	avgenergy[i,k]=numpy.sum(u_kln[k,k,:])/N_k[k]

    dT1 = T_k[1]-T_k[0]
    dT2 = T_k[2]-T_k[1]
    if len(T_k)%2 != 1:
	print("Intermediate points must be an even number for simpsons rule to work!")
	sys.exit()
    if dT2 != dT1:
	print("Intermediate points must be evenly spaced for simpsons rule to work!")
        sys.exit()
    simp = numpy.zeros((K-1)/2,float)
    T_xaxis = numpy.zeros((K+1)/2,float)
    T_xaxis[0]=T_k[0]
    for k in range((K-1)/2):
	if k==0:
	    simp[k] = dT1/3*(avgenergy[i,k]+4*avgenergy[i,k+1]+avgenergy[i,k+2])
        else:
	    simp[k] = simp[k-1]+dT1/3*(avgenergy[i,k]+4*avgenergy[i,k+1]+avgenergy[i,k+2])
	dA[i,k+1] = simp[k]*beta_k[2*k+2]
	ddA[i,k+1] = 0 #How do we calculate error bars with the simpson rule? Reweighting to get derivative?
        T_xaxis[k+1]=T_k[2*k+2]

    print("Free Energies Optained...")

    #=============================================================================================
    # PLOT THE FINAL DATA
    #=============================================================================================
    if (options.plot):    
        # now plot the free energy change as a function of temperature
        plt.figure(2*i)
        xlabel = 'Temperature'
        ylabel = 'Relative Free Energy'
        plt.title(poly)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        Xaxis = T_xaxis
        print('Xaxis:')
        print(Xaxis)
        print('YAxis:')
        print(dA[i,:])
        plt.plot(Xaxis,dA[i,:])
        filename = poly + '_' + str(Molecules) + '_dAvsT.png'
        plt.savefig(filename, bbox_inches='tight')

	# Also plot the raw inputs of u/RT^2 vs T
	plt.figure(2*i+1)
        xlabel = 'Temperature'
        ylabel = 'u/RT^2'
        plt.title(poly)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        Xaxis = T_k
        print('Xaxis:')
        print(Xaxis)
        print('YAxis:')
        print(avgenergy[i,:])
        plt.plot(Xaxis,avgenergy[i,:])
        filename = poly + '_' + str(Molecules) + '_uvsT.png'
        plt.savefig(filename, bbox_inches='tight')

if (options.plot) and options.polymorphs == 'all':
    # now plot the free energy change as a function of temperature
    plt.figure(7)
    xlabel = 'Temperature'
    ylabel = 'Relative Free Energy (kcal/mol)'
    plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    Xaxis = T_xaxis
    print('Xaxis:')
    print(Xaxis)
    print('YAxis:')
    print(dA)
    plt.plot(Xaxis,dA[0,:],'b',Xaxis,dA[1,:],'g',Xaxis,dA[2,:],'r')
    filename = 'ALL_' + str(Molecules) + '_dAvsT.png'
    plt.savefig(filename, bbox_inches='tight')

    # Also plot the raw inputs of u/RT^2 vs T
    plt.figure(8)
    xlabel = 'Temperature'
    ylabel = 'u/RT^2 (K-1)'
    plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    Xaxis = T_k
    print('Xaxis:')
    print(Xaxis)
    print('YAxis:')
    print(avgenergy)
    plt.plot(Xaxis,avgenergy[0,:],'b',Xaxis,avgenergy[1,:],'g',Xaxis,avgenergy[2,:],'r')
    filename = 'ALL_' + str(Molecules) + '_uvsT.png'
    plt.savefig(filename, bbox_inches='tight')
