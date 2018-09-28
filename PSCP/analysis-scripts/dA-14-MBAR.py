from __future__ import print_function
#
# Computing the free energy difference of an organic crystal polymorphs as 14 interactions are removed
# 
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
import MBARBootstrap # Bootstrapping algorithm
import os.path
from optparse import OptionParser # for parsing command-line options
import Harvist #Hamiltonian Reweighting Visualization Toolkit
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-p', '--plot', dest = 'plot', help = 'Plot output (default false)', default=True, action = 'store_true')
parser.add_option('-L', dest = 'Lambda', help = 'Lambda point during the alchemical transformation', default = 100)
parser.add_option('-G', dest = 'Gamma', help = 'Gamma point during the alchemical transformation', default = 000)
parser.add_option('-f', dest = 'exponent', help = 'functional form (exponent) of the spacing between the gammas', default = 2)
parser.add_option('-n', dest = 'polymorphs', help = 'Polymorphs to analyze', default = 'p1 p2')
parser.add_option('-M', dest = 'molecule', help = 'name of the molecule', default = 'benzene')
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
LAMBDA=int(options.Lambda)
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

Colors=['b', 'r', 'g'];

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
    print("Invalid Temperature: " + str(Temp))
    sys.exit()

#LAMBDA
if LAMBDA < 0 or LAMBDA > 100:
    print("Invalid Lambda Point: " + str(LAMBDA))
    sys.exit()

#GAMMA
if GAMMA < 0 or GAMMA > 100:
    print("Invalid Gamma Point: " + str(GAMMA))
    sys.exit()

#POLYMORPH
polymorphs=options.polymorphs.split()
polymorph = []
polymorph_short = []
for i,token in enumerate(polymorphs):
    polymorph.append('Polymorph '+ str(token))
    polymorph_short.append(token)

#POTENTIAL
if potential != "oplsaa" and potential != "gromos" and potential != "designeda" and potential != "oplsaafakeg" and potential != "oplsaafakea":
    print("Invalid Potential")
    print("Supported potentials: oplsaa gromos designeda oplsaafakeg oplsaafakea")
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

#LAMBDA
Lname = ""
if LAMBDA < 10:
    Lname="00" + str(LAMBDA)+"L"
elif LAMBDA < 100:
    Lname="0" + str(LAMBDA)+"L"
else:
    Lname=str(LAMBDA)+"L"

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
if hinge == "DefaultHinge":
    hinges = ['_G']
else:
    # Read in each job
    hinges=[];
    hingevect = options.hinge.split();
    for i,token in enumerate(hingevect):
        hinges.append("_G_" + str(token));
#=============================================================================================
# READ IN RAW DATA
#=============================================================================================
# Constants.
kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184) # Boltzmann constant in kcal/mol

molecule = 'Benzene' #The crystalline molecule
omitT = [] #Temperatures to be omitted from the analysis

# Parameters
T_k = float(Temp)
K = 1  # How many states?
Kbig=K+1
N_max = 200000 # maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
dA = numpy.zeros([len(polymorph),Kbig],float)
ddA = numpy.zeros([len(polymorph),Kbig],float)
convert_units = (0.2390057)*numpy.ones(Kbig,float) #Convert all energies to kcal/mol
#convert_units = (1.0)*numpy.ones(len(Gammas),float) #Convert all energies to kcal/mol
ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================']; #Lines to ignore when reading in energies

for i,poly in enumerate(polymorph):
    # Allocate storage for simulation data
    N_k = numpy.zeros([Kbig],numpy.int32) # N_k[k] is the total number of snapshots from alchemical state k
    N_k_s = numpy.zeros([Kbig,len(hinges)],numpy.int32) # N_k_s[k,s] is the total number of snapshots from alchemical state k from seed s
    u_kln = numpy.zeros([K,Kbig,N_max], numpy.float64) # u_kln[k,l,n] is the adjusted energy of snapshot n from simulation k 
    dhdl_kn = numpy.zeros([K,N_max], numpy.float64) # dhdl_kn[k,n] is the derivative of energy with respect to lambda of snapshot n from simulation k

    #Load in the data for each run            
    for k in range(K):
	n = 0
	for s,hinge in enumerate(hinges):
	    linenum_energy = 0
	    linenum_dhdl = 0
	    # cycle through all the input total energy data
	    #dirpath = '../finishedJobs/' + poly + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lname + '_' + Gamma_names[k] + '_' + Pname + hinge
	    dirpath = '/oldhome/ecd4bd/finishedJobs_archive/PSCP/' + Molecule + '/' + Molecule  + '_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lname + '_' + Gname + '_' + Pname + hinge
	    #fname = '../finishedJobs/' + poly + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lname + '_' + Gamma_names[k] + '_' + Pname + hinge + '/potenergy.xvg'
	    #dhdlname = '../finishedJobs/' + poly + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lname + '_' + Gamma_names[k] + '_' + Pname + hinge + '/benzene_dhdl_PROD.xvg'
	    #fname = '/oldhome/ecd4bd/finishedJobs_archive/' + poly + '_PSCP72/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lname + '_' + Gamma_names[k] + '_' + Pname + hinge + '/potenergy.xvg'
            #dhdlname = '/oldhome/ecd4bd/finishedJobs_archive/' + poly + '_PSCP72/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lname + '_' + Gamma_names[k] + '_' + Pname + hinge + '/benzene_dhdl_PROD.xvg'
	    #fname = '/oldhome/ecd4bd/finishedJobs_archive/' + poly + '_QuadraticGammas/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + str(Molecules) + '_4ind_' + Tname + '_' + Lname + '_' + Gamma_names[k] + '__' + Pname + hinge + '/potenergy.xvg'
	    #dhdlname = '/oldhome/ecd4bd/finishedJobs_archive/' + poly + '_QuadraticGammas/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + str(Molecules) + '_' + Tname + '_' + Lname + '_' + Gamma_names[k] + '__' + Pname + hinge + '/benzene_dhdl_PROD.xvg'
	    fname = dirpath + '/potenergy.xvg'
	    fname14 = dirpath + '/potenergy14.xvg'
            groname = dirpath + '/Traj.gro'
            restname = dirpath + '/restraint.gro'
	    if k not in omitT:
		infile = open(fname, 'r')
		lines = infile.readlines()
		infile.close()
		print("loading " + fname)
		
		ignorecounter=0
		for counter,line in enumerate(lines):
		    tokens_energy = line.split()
		    if tokens_energy[0] in ignore_symbols:
			continue
		    #ignore the first set of frames
		    if ignorecounter <= ignoreframes:
			ignorecounter+=1
			continue

		    #ignore the frames after the include frames
                    if counter > includeframes:
                        continue

		    u_kln[k,:K,n] = (float(Independent)/Molecules)*(float(tokens_energy[1]))*convert_units[k]
		    n+=1
	    if s==0:
                N_k_s[k,s]=n
	    else:
                N_k_s[k,s] = n - sum(N_k_s[k,0:s])
	    #Fill in the additional state with LJ14 interactions removed
	    energies14 = Harvist.GrabTerms(fname14,['LJ-14', 'Coulomb-14'], ignoreframes=ignoreframes)[0]
	    u_kln[k,Kbig-1,:] = u_kln[k,K-1,:]
	    u_kln[k,Kbig-1,:N_k_s[k,s]+1] -= energies14[:,0]*(float(Independent)/Molecules)*convert_units[k]
	    u_kln[k,Kbig-1,:N_k_s[k,s]+1] -= energies14[:,1]*(float(Independent)/Molecules)*convert_units[k]
	N_k[k] = n 

    #data_vector = u_kln[Kbig-12,Kbig-12,0:N_k[Kbig-12]]-u_kln[Kbig-12,Kbig-1,0:N_k[Kbig-12]];
    #plt.hist(data_vector,20,facecolor='b')
    #plt.title('Average: ' + str(numpy.average(data_vector[0:counter])))
    #plt.show()

    # convert to nondimensional units from kcal/mol

    u_kln *=  beta_k
    # all data loaded from the three sets

    u_kln_save = u_kln.copy()
    N_k_save = N_k.copy()
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
    print("Number of retained samples from each seed")
    print(N_k_s)

    #=============================================================================================
    # COMPUTE FREE ENERGY DIFFERENCE USING MBAR
    #=============================================================================================
    
    # Initialize MBAR.
    print("Running MBAR...")

    # generate the weights of each of the umbrella set
    #mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive', use_optimized=False)
    mbar = pymbar.MBAR(u_kln, N_k, verbose = True, subsampling_protocol=[{'method':'L-BFGS-B'}])

    print("MBAR Converged...")
    # testing
    
    for k in range(Kbig):
        w = numpy.exp(mbar.Log_W_nk[:,k])
        print("max weight in state %d is %12.7f" % (k,numpy.max(w)))
        # using Kish (1965) formula.
        # effective # of samples =  (\sum_{i=1}^N w_i)^2 / \sum_{i=1}^N w_i^2
        #                        =  (\sum_{i=1}^N w_i^2)^-1
        neff = 1/numpy.sum(w**2)
        print("Effective number of sample in state %d is %10.3f" % (k,neff))
        print("Efficiency for state %d is %d/%d = %10.4f" % (k,neff,len(w),neff/len(w)))

    # extract self-consistent weights and uncertainties
    (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()

    print("Free Energies Optained...")

    #convert PMF to kcal/mol and normalize by the number of molecules
    df_i /= (beta_k*float(Independent)) 
    ddf_i /= (beta_k*float(Independent))

    dA[i,:] = df_i[0]
    #ddA[i,:] = ddf_i[0]

        
    #=============================================================================================
    # COMPUTE UNCERTAINTY USING THE UNCORRELATED DATA
    #=============================================================================================
    
    #for k in range(K):
    #    N_k[k]=0
    #    n_old=0
    #    if k not in omitT:
    #        for s in range(len(hinges)):
    #    	#g_k[k] = timeseries.statisticalInefficiency(u_kln[k,k,0:N_k[k]])
    #    	g_k[k] = timeseries.statisticalInefficiency(dhdl_kn[k,n_old:(n_old+N_k_s[k,s])])
    #    	print "Correlation time for sampled state %d is %10.3f" % (k,g_k[k])
    #    	# subsample the data to get statistically uncorrelated data
    #    	indices = numpy.array(timeseries.subsampleCorrelatedData(u_kln[k, k, n_old:(n_old+N_k_s[k,s])], g=g_k[k]))  # subsample
    #    	u_kln[k,:,N_k[k]:(N_k[k]+len(indices))] = u_kln_save[k,:,(indices+n_old)].transpose()  # not sure why we have to transpose
    #    	N_k[k] = N_k[k]+len(indices)
    #    	n_old+=N_k_s[k,s]
    #    
    #print "Number of retained samples"
    #print N_k
    #print "Number of retained samples from each seed"
    #print N_k_s

    # generate the weights of each of the umbrella set
    #mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive', use_optimized=False)
    mbar = pymbar.MBAR(u_kln, N_k, verbose = True, subsampling_protocol=[{'method':'L-BFGS-B'}])

    print("MBAR Converged...")
    # testing

    # extract self-consistent weights and uncertainties
    (df_u, ddf_u, theta_i) = mbar.getFreeEnergyDifferences()

    print("Free Energies Optained...")

    #convert PMF to kcal/mol and normalize by the number of molecules
    df_u /= (beta_k*float(Independent))
    ddf_u /= (beta_k*float(Independent))

    ddA[i,:] = ddf_u[0]
    
    # Write out free energy differences
    print("Free Energy Difference (in units of kcal/mol)")
    for k in range(Kbig):
        print("%8.3f %8.3f" % (-df_i[k,0], ddf_u[k,0]))
      
    """ 
    #Calculate the uncertainties using bootstrapping
    indexVect = numpy.zeros([2,Kbig-1], numpy.float)
    #indexVect[:,0] = [0, Kbig-1];
    #indexVect[:,1] = [0, 10];
    #indexVect[:,2] = [0, 1];
    indexVect[0,:] = 0
    indexVect[1,:] = numpy.arange(Kbig-1)+1
    if len(hinges) > 1:
	datafile = 'BootstrapData_' + PotNAME + '_' + poly + '_' + Molname + Tname + '_' + Pname + '_dAvsG_All'
        stdfile = 'BootstrapStd_' + PotNAME + '_' + poly + '_' + Molname + Tname + '_' + Pname + '_dAvsG_All'
    else:
	datafile = 'BootstrapData_' + PotNAME + '_' + poly + '_' + Molname + Tname + '_' + Pname + '_dAvsG' + hinge
        stdfile = 'BootstrapStd_' + PotNAME + '_' + poly + '_' + Molname + Tname + '_' + Pname + '_dAvsG' + hinge
    MBARBootstrap.runMBARBootstrap(u_kln, N_k, beta_k, Independent, indexVect, datafile, stdfile, 200)
    """


#=============================================================================================
# PRINT THE FINAL DATA
#=============================================================================================
for i,poly in enumerate(polymorph):
    print(poly + ": " + "%8.3f %8.3f" % (dA[i,Kbig-1], ddA[i,Kbig-1]))


#=============================================================================================
# PLOT THE FINAL DATA
#=============================================================================================
"""
for k in range(K):
    #if k==0:
    #    fig=plt.figure(k+2)
    #    xlabel = 'Energy Difference (kcal/mol)'
    #    ylabel = 'Frequency'
    #    plt.title(polymorph_short[0] + ' dU Full Intramolecular and no Intramolecular')
    #    #plt.title(polymorph_short[0]+' sampled in ' + potentials[k])
    #    plt.xlabel(xlabel)
    #    plt.ylabel(ylabel)
    ##if omitk[k]==0:
    #    #continue
    ##Energy_Differences = (u_kln_save[k,11,ignoreframes:N_k_save[k]-1] - u_kln_save[k,10,ignoreframes:N_k_save[k]-1])/float(Molecules)
    ##print Energy_Differences
    #print 'Xaxis:'
    #print Xaxis_OPLS
    #print 'OPLS Energy:'
    #print Energy_OPLS
    #iplt.hist(Energy_Differences,10,facecolor=Colors[0],alpha=0.2,label=str(k))
    ##n,bins,patches = plt.hist(Energy_Differences,10,facecolor='g',alpha=0.5,label=potentials_short[k])
    #plt.hold(True)
    #plt.legend(loc='upper right')
    #plt.set_xlim([0,5])
    #Histogram_data, bins = numpy.histogram(Energy_Differences,10)
    #print 'Amoeba Energy:'
    #print Energy_Amoeba
    #pdb.set_trace()
    #ax.errorbar(Xaxis,dA[0,:], color='b', yerr=ddA[0,:])
    #ax.errorbar(Xaxis,dA[1,:], color='g', yerr=ddA[1,:])
    #ax.errorbar(Xaxis,dA[2,:], color='r', yerr=ddA[2,:])
    #plt.errorbar(Xaxis,dA[0,:],ddA[0,:],Xaxis,dA[1,:],ddA[1,:],Xaxis,dA[2,:],ddA[2,:])
    #plt.plot(bins,Histogram_data,'b')
    #plt.errorbar(Xaxis,dA[0,:], color='b', yerr=ddA[0,:],Xaxis,dA[1,:],color='g', yerr=ddA[1,:],Xaxis,dA[2,:],color='r',yerr=ddA[2,:])
    #filename = polymorph_short[0] + ' ' + potentials[k] + ' dU Histogram.png'
    #plt.savefig(filename, bbox_inches='tight')
plt.legend(loc='upper right')
plt.show()
"""




if (options.plot) and options.polymorphs == 'all':
    # now plot the free energy change as a function of temperature
    fig = plt.figure(4)
    ax = fig.add_subplot(111)
    #xlabel = 'Gamma Point'
    #ylabel = 'Relative Free Energy'
    xlabel = 'Interaction Strength, $\gamma$'
    ylabel = 'Relative Free Energy (kcal/mol)'
    #if len(hinges) > 1:
#	plt.title('All Polymorphs Combined Runs')
#    else:
#        plt.title('All Polymorphs ' + hinge)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    Xaxis = [float(j/100.0) for j in Gammas]
   
    if os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph1_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsG_All'):
        ddA[0,:]=MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph1_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsG_All')
    elif os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsG_All'):
        ddA[1,:]=MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsG_All')
    elif os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsG_All'):
        ddA[2,:]=MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph3_' + str(Molecules) + '_' + Tname + '_' + Pname + '_dAvsG_All') 

    ax.errorbar(Xaxis,dA[0,:], color='b', yerr=ddA[0,:],label='Benzene I')
    ax.errorbar(Xaxis,dA[1,:], color='g', yerr=ddA[1,:],label='Benzene II')
    ax.errorbar(Xaxis,dA[2,:], color='r', yerr=ddA[2,:],label='Benzene III')
    plt.legend(loc='upper right')
    #plt.errorbar(Xaxis,dA[0,:],ddA[0,:],Xaxis,dA[1,:],ddA[1,:],Xaxis,dA[2,:],ddA[2,:])
    #plt.plot(Xaxis,dA[0,:],'b',Xaxis,dA[1,:],'r',Xaxis,dA[2,:],'g')
    #plt.errorbar(Xaxis,dA[0,:], color='b', yerr=0.1,Xaxis,dA[1,:],color='g', yerr=0.1,Xaxis,dA[2,:],color='r',yerr=0.1)
    if len(hinges) > 1:
	filename = PotNAME + '_' + str(Molecules) + '_' + Tname + '_dAvsG.pdf'
    else:
   	filename = PotNAME + '_' + str(Molecules) + '_' + Tname + hinge + '_dAvsG.pdf'
    plt.savefig(filename, bbox_inches='tight')
    #plt.show()
