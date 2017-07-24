#
# Computing the free energy difference of an organic crystal polymorph at different lambda values
# 
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
import MBARBootstrap # Bootstrapping algorithm
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import MBARBootstrap
import os.path
import pdb
import tossconfigurationsFunc

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-p', '--plot', dest = 'plot', help = 'Plot output (default false)', default=True, action = 'store_true')
parser.add_option('-U', dest = 'MinLambda', help = 'Minimum lambda point sampled', default = 0)
parser.add_option('-R', dest = 'MaxLambda', help = 'Maximum lambda point sampled', default = 100)
parser.add_option('-L', dest = 'spacing', help = 'Spacing between lambda points', default = 5)
parser.add_option('-G', dest = 'Gamma', help = 'Gamma point during the alchemical transformation', default = 100)
parser.add_option('-f', dest = 'exponent', help = 'functional form (exponent) of the spacing between the lambdas', default = 4)
parser.add_option('-M', dest = 'molecule', help = 'name of the molecule', default = 'benzene')
parser.add_option('-n', dest = 'polymorphs', help = 'Polymorphs to analyze', default = 'p1 p2')
parser.add_option('-N', dest = 'molecules', help = 'number of supercell molecules', default = 72)
parser.add_option('-I', dest = 'independent', help = 'number of independent molecules', default = 4)
parser.add_option('-T', dest = 'Temperature', help = 'Temperature', default = 200)
parser.add_option('-P', dest = 'Pressure', help = 'Pressure', default = 1)
parser.add_option('-k', dest = 'forceConstant', help = 'Harmonic Force Constant', default = 1000)
parser.add_option('-i', dest = 'ignoreframes', help = 'Initial frames to ignore', default = 2000)
parser.add_option('-j', dest = 'includeframes', help = 'Number of frames to include', default = 100000)
parser.add_option('-u', dest = 'potential', help = 'Potential used in the simulation', default = 'oplsaa')
parser.add_option('-H', '--hinge', dest = 'hinge', help = 'Optional string at end of jobs', default = 'DefaultHinge')

(options, args) = parser.parse_args()
MinL = float(options.MinLambda)
MaxL = float(options.MaxLambda)
dL = float(options.spacing)
GAMMA=int(options.Gamma)
exponent = int(options.exponent)
Temp = float(options.Temperature)
Pressure = int(options.Pressure)
Molecule = options.molecule
Molecules = int(options.molecules)
Independent = int(options.independent)
k = float(options.forceConstant)
ignoreframes = int(options.ignoreframes)
includeframes = int(options.includeframes)
potential = str(options.potential)
hinge = options.hinge

if (options.plot):
    import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.font_manager import FontProperties as FP
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 16}
    matplotlib.rc('font', **font)

#=============================================================================================
# ENSURE THAT USER INPUTS ARE SENSIBLE
#=============================================================================================
#TEMPERATURE
if Temp < 0:
    print "Invalid Temperature: " + str(Temp)
    sys.exit()

if Pressure < 0:
    print "Invalid Pressure: " + str(Pressure)
    sys.exit()

#LAMBDA
if (MinL == -1 ) and (MaxL == -1) and (dL == -1) and (exponent == 1):
    print "Using default values!"
    Lambdas = ['000L','010L','020L','030L','040L','050L','060L','070L','080L','090L','100L'] #The Lambda points sampled
elif MinL < 0 or MaxL < 0 or dL < 0 or MinL > MaxL:
    print "Invalid Lambda Specifications"
    sys.exit()
else:
    RawLambda = 0
    Lambdas = [];
    Lambda_names = [];
    Lambda_indicies=[];
    index=0
    while RawLambda < MaxL:
	if RawLambda >= MinL:
	    Lambda_indicies.append(index)
	    index+=1
	else:
	    index+=1
	    RawLambda=RawLambda+dL
	    continue
	Lambda=int(100*float(RawLambda**exponent)/float(MaxL**exponent))
	Lambdas.append(Lambda)
	#Format the lambda point name
	if RawLambda < 10:
            Lambda_names.append('00' + str(int(RawLambda)) + 'L')
	elif RawLambda < 100:
            Lambda_names.append('0' + str(int(RawLambda)) + 'L')
	else:
            Lambda_names.append('100L')
	RawLambda=RawLambda+dL
    #Catch the final lambda point
    Lambdas.append(MaxL)
    Lambda_indicies.append(index)
    if MaxL < 10:
        Lambda_names.append('00' + str(int(MaxL)) + 'L')
    elif MaxL < 100:
        Lambda_names.append('0' + str(int(MaxL)) + 'L')
    else:
        Lambda_names.append('100L')

#GAMMA
if GAMMA < 0 or GAMMA > 100:
    print "Invalid Gamma Point: " + str(GAMMA)
    sys.exit()

#POLYMORPH
polymorphs=options.polymorphs.split()
polymorph = []
polymorph_short = []
for i,token in enumerate(polymorphs):
    polymorph.append('Polymorph '+ str(token))
    polymorph_short.append(token)
#if (options.polymorphs == 'all'):
#    polymorph = ['Polymorph1', 'Polymorph2', 'Polymorph3']
#    polymorph_short = ['p1', 'p2', 'p3']
#    if Molecule != "benzene" and Molecule != "glycin":
#        polymorph = ['Polymorph1', 'Polymorph2']
#        polymorph_short = ['p1', 'p2']
#elif (options.polymorphs == 'p1'):
#    polymorph = ['Polymorph1']
#    polymorph_short = ['p1']
#elif (options.polymorphs == 'p2'):
#    polymorph = ['Polymorph2']
#    polymorph_short = ['p2']
#elif (options.polymorphs == 'p3'):
#    polymorph = ['Polymorph3']
#    polymorph_short = ['p3']
#elif (options.polymorphs == 'p4'):
#    polymorph = ['Polymorph4']
#    polymorph_short = ['p4']
#else:
#    print "Polymorph Inputs Wrong"
#    sys.exit()

#POTENTIAL
if potential != "oplsaa" and potential != "gromos" and potential != "designeda" and potential != "oplsaafakeg" and potential != "oplsaafakea":
    print "Invalid Potential"
    print "Supported potentials: oplsaa gromos designeda oplsaafakeg oplsaafakea"
    sys.exit()


#=============================================================================================
# FORMAT INPUTS
#=============================================================================================
#TEMPERATURE
Tname = ""
if Temp < 10:
    Tname="00" + str(int(Temp))+"K"
elif Temp < 100:
    Tname="0" + str(int(Temp))+"K"
else:
    Tname=str(int(Temp))+"K"

#PRESSURE
Pname = ""
if Pressure < 10:
    Pname="00" + str(int(Pressure))+"P"
elif Pressure < 100:
    Pname="0" + str(int(Pressure))+"P"
else:
    Pname=str(int(Pressure))+"P"

#GAMMA POINT
Gname = ""
if GAMMA < 10:
    Gname="00" + str(int(GAMMA))+"G"
elif GAMMA < 100:
    Gname="0" + str(int(GAMMA))+"G"
else:
    Gname=str(int(GAMMA))+"G"

#NUMBER OF MOLECULES
Molname = ""
if Molecules == Independent:
    Molname = str(Molecules) + '_'
else:
    Molname = str(Molecules) + '_' + str(Independent) + 'ind_'


#POTENTIAL
PotNAME = ""
if potential == "oplsaa":
    PotNAME = "OPLS"
elif potential == "gromos":
    PotNAME = "GROM"
elif potential == "designeda":
    PotNAME = "DESA"
elif potential == "oplsaafakeg":
    PotNAME = "FAKEG"
elif potential == "oplsaafakea":
    PotNAME = "FAKEA"

#CHARGE AND SIGMA HINGE
if potential == "oplsaa":
    ChargeHinge=""
elif potential == "gromos":
    ChargeHinge=""
elif potential == "designeda":
    ChargeHinge=""
elif potential == "oplsaafakeg":
    ChargeHinge="_C01150"
elif potential == "oplsaafakea":
    ChargeHinge="_C01150"

# OPTIONAL HINGE
if str(GAMMA) == "100":
    hingeLetter="L"
else:
    hingeLetter="R"

if hinge == "DefaultHinge":
    hinges = ["_" + hingeLetter];
else:
    # Read in each job
    hinges=[];
    hingevect = options.hinge.split();
    for i,token in enumerate(hingevect):
        hinges.append("_" + hingeLetter + "_" + str(token));
#=============================================================================================
# READ IN RAW DATA
#=============================================================================================
# Constants.
kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184) # Boltzmann constant in kcal/mol

omitK = []
#omitK = [0,1,2,3,4,5,6,7,8]
#pdb.set_trace()
#omitK = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18] #Temperatures to be omitted from the analysis
# Parameters
T_k = Temp*numpy.ones(len(Lambdas),float) #Convert temperatures to floats 
print T_k
print Lambdas
#seeds = [201]; #The random seed used (not included at the moment)
g_k = numpy.zeros([len(Lambdas)],float)
K = len(Lambdas)  # How many states?
Kbig = K+0 # total number of states examined; none are unsampled
N_max = 200000 # maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
dA = numpy.zeros([len(polymorph),len(Lambdas)],float)
ddA = numpy.zeros([len(polymorph),len(Lambdas)],float)
convert_units = (0.2390057)*numpy.ones(len(Lambdas),float) #Convert all energies to kcal/mol
#convert_units = (1.0)*numpy.ones(len(Lambdas),float) #Convert all energies to kcal/mol
ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================']; #Lines to ignore when reading in energies

for i,poly in enumerate(polymorph):
    # Allocate storage for simulation data
    N_k = numpy.zeros([Kbig],numpy.int32) # N_k[k] is the total number of snapshots from alchemical state k
    N_ksj = numpy.zeros([Kbig,len(hinges),100],numpy.int32) # N_k_s[k,s] is the total number of snapshots from alchemical state k from seed s in 'unflipped segment j'

    u_kln = numpy.zeros([K,Kbig,N_max], numpy.float64) # u_kln[k,l,n] is the adjusted energy of snapshot n from simulation k
    dhdl_kln = numpy.zeros([K,Kbig,N_max], numpy.float64) # dhdl_kln[k,l,n] is the restraint energy value of snapshop n from simulation k
    dhdl_kn = numpy.zeros([K,N_max], numpy.float64) # dhdl_kn[k,n] is the derivative of energy with respect to lambda of snapshot n from simulation k 

    #Load in the data for each run
    for k in range(K):
        n=0     
	for s,hinge in enumerate(hinges):
	    tossconfigs=[] #The index of each configuration to toss from the MBAR analysis
	    keepconfigs=numpy.arange(N_max) #The index of each configuration to keep in the MBAR analysis
	    linenum_energy = 0
	    linenum_dhdl = 0
	    # cycle through all the input total energy data
	    #dirpath='../finishedJobs/' + poly + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lambda_names[k] + '_' + Gname + '_' + Pname + hinge
	    dirpath='/oldhome/ecd4bd/finishedJobs_archive/PSCP/' + Molecule + '/' + Molecule + '_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lambda_names[k] + '_' + Gname + '_' + Pname + hinge
	    fname = dirpath + '/potenergy.xvg'
	    dhdlname = dirpath + '/dhdl_PROD.xvg'
	    groname = dirpath + '/Traj.gro'
	    outname = dirpath + '/Flipless.gro'
	    restname = dirpath + '/restraint.gro'
	    if k not in omitK:
		infile = open(fname, 'r')
		lines = infile.readlines()
		infile.close()
		print "loading " + fname
		infile = open(dhdlname, 'r')
		lines_dhdl = infile.readlines()
		infile.close()
		print "loading " + dhdlname

	        """ 		
		#Toss configurations if they are 'flipped'
		if os.path.isfile(dirpath+'/tossconfigs.npy'):
		    tossconfigs=numpy.load(dirpath+'/tossconfigs.npy')
		    print tossconfigs
		elif polymorph_short[i]=='p3':
		    print "Polymorph 3!"
		else:
		    tossconfigs = tossconfigurationsFunc.tossconfigs(groname, outname, restname)
		    if tossconfigs != []:
		        tossconfigs = (tossconfigs-1)
		    numpy.save(dirpath+'/tossconfigs',numpy.array(tossconfigs))
			

		#Remove tossed configurations from the keep list	
		keepconfigs=[j for j in keepconfigs if j not in tossconfigs];
		"""

		ignorecounter=0
		symbolcounter=0
		for counter,line in enumerate(lines):
		    tokens_energy = line.split()
		    if tokens_energy[0] in ignore_symbols:
			symbolcounter+=1
			continue
		    #ignore the first set of frames
		    if ignorecounter <= ignoreframes:
			ignorecounter+=1
			continue

		    #ignore the frames after the include frames
		    if counter > includeframes:
			continue
	
		    #ignore frames that are not in the trr file
		    #if counter % 200 != 0:
			#continue

		    #ignore frames that are flipped and we are tossing
		    if counter in tossconfigs:
			continue

	
		    #Grab the dhdl information (if possible)
		    tokens_dhdl = lines_dhdl[linenum_dhdl].split()
		    while tokens_dhdl[0] in ignore_symbols:
			linenum_dhdl+=1
			tokens_dhdl = lines_dhdl[linenum_dhdl].split()
		    #if not tokens_dhdl[0].isdigit(): 
			#continue
		    while float(tokens_energy[0]) != float(tokens_dhdl[0]) and (linenum_dhdl+1) < len(lines_dhdl) and linenum_dhdl < 1000000:
			linenum_dhdl+=1
			tokens_dhdl = lines_dhdl[linenum_dhdl].split()
		    #print tokens_dhdl
		    if float(tokens_energy[0]) != float(tokens_dhdl[0]):
			#print "Steps not equal for energy and dhdl!!"
			#print "Energy Step: " + tokens[0]
			#print "DHDL Step: " + tokens_dhdl[0]
			#sys.exit()
			continue
		    # the energy of every configuration from each state evaluated at its sampled state
		    u_kln[k,:,n] = (float(Independent)/Molecules)*(float(tokens_energy[1]) + numpy.asarray(numpy.array(tokens_dhdl)[5+numpy.array(Lambda_indicies)],float))*convert_units[k]
		    dhdl_kln[k,:,n] = numpy.asarray(numpy.array(tokens_dhdl)[5+numpy.array(Lambda_indicies)],float)*convert_units[k]
		    dhdl_kn[k,n] = (float(Independent)/Molecules)*float(tokens_dhdl[4])*convert_units[k]  

		    n+=1
		
		#Truncate the kept configuration list to be less than n
		keepconfigs = [j for j in keepconfigs if j < (counter-symbolcounter) and j >= ignoreframes]
		
		#Split up the retained configurations into connected segments
		j=0
		a=0
		for a in range(len(keepconfigs)):
		    if a==0:
			continue
		    elif int(keepconfigs[a-1])+1 != int(keepconfigs[a]):
			N_ksj[k,s,j] = a-(sum(N_ksj[k,s,0:j]))
			j+=1
		#Catch the final segment
		N_ksj[k,s,j] = len(keepconfigs)-sum(N_ksj[k,s,0:j])
		j+=1		

		
		#if s==0:
		    #N_k_s[k,s]=n
		#else:
		    #N_k_s[k,s] = n - sum(N_k_s[k,0:s])
        N_k[k] = n

    # convert to nondimensional units from kcal/mol

    u_kln *=  beta_k[0]
    # all data loaded from the three sets

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

    #Ignore the first state due to jumping

    print "Number of retained samples"
    print N_k

    #=============================================================================================
    # COMPUTE FREE ENERGY DIFFERENCE USING MBAR
    #=============================================================================================

    # Initialize MBAR.
    print "Running MBAR..."

    # generate the weights of each of the umbrella set
    #mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive', use_optimized=False)
    mbar = pymbar.MBAR(u_kln, N_k, verbose = True, subsampling_protocol=[{'method':'L-BFGS-B'}])

    print "MBAR Converged..."
    # testing

    for k in range(Kbig):
        w = numpy.exp(mbar.Log_W_nk[:,k])
        print "max weight in state %d is %12.7f" % (k,numpy.max(w))
        # using Kish (1965) formula.
        # effective # of samples =  (\sum_{i=1}^N w_i)^2 / \sum_{i=1}^N w_i^2
        #                        =  (\sum_{i=1}^N w_i^2)^-1
        neff = 1/numpy.sum(w**2)
        print "Effective number of sample in state %d is %10.3f" % (k,neff)
        print "Efficiency for state %d is %d/%d = %10.4f" % (k,neff,len(w),neff/len(w))

    # extract self-consistent weights and uncertainties
    (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()

    print "Free Energies Optained..."

    #convert PMF to kcal/mol and normalize by the number of molecules
    df_i /= (beta_k[0]*float(Independent))
    ddf_i /= (beta_k[0]*float(Independent))

    dA[i,:] = df_i[0]
    #ddA[i,:] = ddf_i[0]
    
    #=============================================================================================
    # COMPUTE UNCERTAINTY USING THE UNCORRELATED DATA
    #=============================================================================================

    for k in range(K):	#For each restraint state
	N_k[k]=0
        n_old=0
        if k not in omitK:
	    for s in range(len(hinges)):	#For each independent trajectory of this restraint state
		for j in range(100):					#For each untossed segment of each independent trajectory of this restraint state
               	    if N_ksj[k,s,j] == 0:
			continue
		    
		    #g_k[k] = timeseries.statisticalInefficiency(u_kln[k,k,0:N_k[k]])
		    g_k[k] = timeseries.statisticalInefficiency(dhdl_kn[k,n_old:(n_old+N_ksj[k,s,j])]) #Feed in the segment and calculate correlation time
		    print "Correlation time for sampled state %d is %10.3f" % (k,g_k[k])
		    # subsample the data to get statistically uncorrelated data
		    indices = numpy.array(timeseries.subsampleCorrelatedData(u_kln[k, k, n_old:(n_old+N_ksj[k,s,j])], g=g_k[k]))  # subsample indices within the segment
		    u_kln[k,:,N_k[k]:(N_k[k]+len(indices))] = u_kln_save[k,:,(indices+n_old)].transpose()  # Apphend the uncorrelated configurations in the segment to the u_kln matrix
		    N_k[k] = N_k[k]+len(indices)
		    n_old+=N_ksj[k,s,j]

    print "Number of retained samples"
    print N_k
    print "Number of retained samples from each seed"
    print N_ksj
	
    # generate the weights of each of the umbrella set
    #mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive', use_optimized=False)
    mbar = pymbar.MBAR(u_kln, N_k, verbose = True, subsampling_protocol=[{'method':'L-BFGS-B'}])

    print "MBAR Converged..."
    # testing

    # extract self-consistent weights and uncertainties
    (df_u, ddf_u, theta_i) = mbar.getFreeEnergyDifferences()

    print "Free Energies Optained..."

    #convert PMF to kcal/mol and normalize by the number of molecules
    df_u /= (beta_k[0]*float(Independent))
    ddf_u /= (beta_k[0]*float(Independent))

    ddA[i,:] = ddf_u[0]

    # Write out free energy differences
    print "Free Energy Difference (in units of kcal/mol)"
    for k in range(Kbig):
        print "%8.3f %8.3f" % (-df_i[k,0], ddf_u[k,0])
    
    ##Calculate the uncertainties using bootstrapping
    # 
    #indexVect = numpy.zeros([2,Kbig-1], numpy.float)
    ##indexVect[:,0] = [0, Kbig-1];
    ##indexVect[:,1] = [0, 10];
    ##indexVect[:,2] = [0, 1];
    #indexVect[0,:] = 0
    #indexVect[1,:] = numpy.arange(Kbig-1)+1
    #if len(hinges) > 1:
    #    datafile = 'BootstrapData/BootstrapData_' + Molecule + '_' + polymorph_short[i] + '_' + Molname + Tname + '_' + Pname + '_' + PotNAME + '_dAvsL_All'
    #    stdfile = 'BootstrapData/BootstrapStd_' + Molecule + '_' + polymorph_short[i] + '_' + Molname + Tname + '_' + Pname + '_' + PotNAME + '_dAvsL_All' 
    #else:
    #    datafile = 'BootstrapData/BootstrapData_' + Molecule + '_' + polymorph_short[i] + '_' + Molname + Tname + '_' + Pname + '_' + PotNAME + hinge
    #    stdfile = 'BootstrapData/BootstrapStd_' + Molecule + '_' + polymorph_short[i] + '_' + Molname + Tname + '_' + Pname + '_' + PotNAME + hinge
  
    #if not os.path.isfile(stdfile+'.txt'): 
    #    MBARBootstrap.runMBARBootstrap(u_kln, N_k, beta_k, Independent, indexVect, datafile, stdfile, 200) 
    

    #=============================================================================================
    # PLOT THE FINAL DATA
    #=============================================================================================
    """
    if (options.plot):    
        # now plot the free energy change as a function of temperature
        plt.figure(i)
        xlabel = 'Lambda Point'
        ylabel = 'Relative Free Energy'
        plt.title(poly)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        Xaxis = Lambdas
        print 'Xaxis:'
        print Xaxis
        print 'YAxis:'
        print dA[i,:]
        plt.errorbar(Xaxis,dA[i,:], ddA[i,:])
        filename = poly + '_' + str(Molecules) + '_' + Tname + '_dAvsL.png'
        plt.savefig(filename, bbox_inches='tight')
	plt.show()
    """    

#=============================================================================================
# PRINT THE FINAL DATA
#=============================================================================================
for i,poly in enumerate(polymorph):
    print poly + ": " + "%8.3f %8.3f" % (-dA[i,Kbig-1], ddA[i,Kbig-1])

#=============================================================================================
# PLOT THE FINAL DATA
#=============================================================================================
if (options.plot) and options.polymorphs == 'all':
    # now plot the free energy change as a function of temperature
    fig=plt.figure(4)
    ax=fig.add_subplot(111)
    #xlabel = 'Restraint Strength (%)'
    xlabel = 'Restraint Strength, $\lambda$'
    ylabel = 'Relative Free Energy (kcal/mol)'
    #if len(hinges) > 1:
#	plt.title('All Polymorphs Combined Runs')
#    else:
#	plt.title('All Polymorphs'+hinge)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    Xaxis = [float(j/100.0) for j in Lambdas]
    
    if os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph1_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsL_All'):
	ddA[0,:]=MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph1_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsL_All')
    elif os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsL_All'):
        ddA[1,:]=MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsL_All')
    elif os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsL_All'):
        ddA[2,:]=MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph3_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsL_All')
 
    ax.errorbar(Xaxis,dA[0,:], color='b', yerr=ddA[0,:], label='Benzene I')
    ax.errorbar(Xaxis,dA[1,:], color='g', yerr=ddA[1,:], label='Benzene II')
    ax.errorbar(Xaxis,dA[2,:], color='r', yerr=ddA[2,:], label='Benzene III')
    plt.legend(loc='upper left')
    #ax.set_ylim([0,3.0])
    #plt.errorbar(Xaxis,dA[0,:],ddA[0,:],Xaxis,dA[1,:],ddA[1,:],Xaxis,dA[2,:],ddA[2,:])
    #plt.plot(Xaxis,dA[0,:],'b',Xaxis,dA[1,:],'r',Xaxis,dA[2,:],'g')
    #plt.errorbar(Xaxis,dA[0,:], color='b', yerr=ddA[0,:],Xaxis,dA[1,:],color='g', yerr=ddA[1,:],Xaxis,dA[2,:],color='r',yerr=ddA[2,:])
    if len(hinges) > 1:
        filename = PotNAME + '_' + str(Molecules) + '_' + Tname + '_dAvsL.pdf'
    else:
        filename = PotNAME + '_' + str(Molecules) + '_' + Tname + hinge+ '_dAvsL.pdf'
    #plt.savefig(filename, bbox_inches='tight')
    plt.show()


