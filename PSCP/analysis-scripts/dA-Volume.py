from __future__ import print_function
#
# Computing the free energy difference of an organic crystal polymorph at different volumes
# 
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
import calculate_gro_volume
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import pdb

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-p', '--plot', dest = 'plot', help = 'Plot output (default false)', default=True, action = 'store_true')
parser.add_option('-S', dest = 'MinVolume', help = 'Minimum box vector scaling', default = 840)
parser.add_option('-B', dest = 'MaxVolume', help = 'Maximum box vector scaling', default = 890)
parser.add_option('-s', dest = 'spacing', help = 'Spacing between lambda points', default = 5)
parser.add_option('-G', dest = 'Gamma', help = 'Gamma Point at which  volume was changed', default = 100)
parser.add_option('-n', dest = 'polymorphs', help = 'Polymorphs to analyze', default = 'all')
parser.add_option('-N', dest = 'molecules', help = 'number of supercell molecules', default = 72)
parser.add_option('-I', dest = 'independent', help = 'number of independent molecules', default = 4)
parser.add_option('-i', dest = 'ignoreframes', help = 'Initial frames to ignore', default = 1000)
parser.add_option('-j', dest = 'includeframes', help = 'Number of frames to include', default = 100000)
parser.add_option('-u', dest = 'potential', help = 'Potential used in the simulation', default = 'oplsaa')
parser.add_option('-T', dest = 'Temperature', help = 'Temperature', default = 200)
parser.add_option('-P', dest = 'Pressure', help = 'Pressure', default = 1)
parser.add_option('-d', dest = 'directory', help = 'Parent directory of the volume change directories', default = 'Volume_Coarse')
parser.add_option('-H', '--hinge', dest = 'hinge', help = 'Optional string at end of jobs', default = 'DefaultHinge')

(options, args) = parser.parse_args()
MinV = float(options.MinVolume)
MaxV = float(options.MaxVolume)
dV = float(options.spacing)
GAMMA = int(options.Gamma)
Temp = float(options.Temperature)
Pressure = int(options.Pressure)
Molecules = int(options.molecules)
Independent = int(options.independent)
ignoreframes = int(options.ignoreframes)
includeframes = int(options.includeframes)
potential = str(options.potential)
directory = options.directory
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

if Temp < 0:
    print("Invalid Temperature: " + str(Temp))
    sys.exit()

if (MinV == -1 ) and (MaxV == -1) and (dV == -1):
    print("Using default values!")
    Volumes = ['v100','v102','v104','v106','v108','v110','v112','v114','v116','v118','v120'] #The scaling parameters sampled
elif MinV < 0 or MaxV < 0 or dV < 0 or MinV > MaxV:
    print("Invalid Volume Specifications")
    sys.exit()
else:
    Volume = MinV
    Volumes = [];
    Volume_names = [];
    while Volume < MaxV:
        Volumes.append(Volume)
        #Format the  name
        if Volume < 10:
            Volume_names.append('000' + str(int(Volume))+'V')
        elif Volume < 100:
            Volume_names.append('00' + str(int(Volume))+'V')
	#elif Volume < 1000:
	    #Volume_names.append('0' + str(int(Volume))+'V')
        else:
            Volume_names.append(str(int(Volume))+'V')
        Volume=Volume+dV
    #Catch the final Volume point
    Volumes.append(MaxV)
    if MaxV < 10:
        Volume_names.append('000' + str(int(Volume))+'V')
    elif MaxV < 100:
        Volume_names.append('00' + str(int(Volume))+'V')
    #elif MaxV < 1000:
        #Volume_names.append('0' + str(int(Volume))+'V')
    else:
        Volume_names.append(str(int(Volume))+'V')

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

#POTENTIAL
if potential != "oplsaa" and potential != "gromos" and potential != "oplsaafakeg" and potential != "oplsaafakea":
    print("Invalid Potential")
    print("Supported potentials: oplsaa gromos oplsaafakeg oplsaafakea")
    sys.exit()

#=============================================================================================
# FORMAT INPUTS
#=============================================================================================
#POTENTIAL
PotNAME = ""
if potential == "oplsaa":
    PotNAME = "OPLS"
    Chargehinge=""
elif potential == "gromos":
    PotNAME = "GROM"
    Chargehinge=""
elif potential == "oplsaafakeg":
    PotNAME = "FAKEG"
    Chargehinge="_C01150_C100H100"
elif potential == "oplsaafakea":
    PotNAME = "FAKEA"
    Chargehinge="_C01150_C100H100"

#NUMBER OF MOLECULES
Molname = ""
if Molecules == Independent:
    Molname = str(Molecules) + '_'
else:
    Molname = str(Molecules) + '_' + str(Independent) + 'ind_'

#PRESSURE
Pname = ""
if Pressure < 10:
    Pname="00" + str(int(Pressure))+"P"
elif Pressure < 100:
    Pname="0" + str(int(Pressure))+"P"
else:
    Pname=str(int(Pressure))+"P"

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

gname = "benzene_pre_EM.gro"
molecule = 'Benzene' #The crystalline molecule
ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================']; #Lines to ignore when reading in energies

# Parameters
T_k = Temp*numpy.ones(len(Volumes),float) #Convert temperatures to floats
P = numpy.zeros([3,len(Volumes)],float) #Convert Pressures to floats
V = Temp*numpy.ones(len(Volumes),float) #Convert Volumes to floats
K = len(T_k)  # How many states?
Kbig = K+0 # total number of states examined; none are unsampled
N_max = 200000 # maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
dA = numpy.zeros([3,(K+1)/2],float)
ddA = numpy.zeros([3,(K+1)/2],float) 
#seeds = [201]; #The random seed used (not included at the moment)
g_k = numpy.zeros([len(T_k)],float)

for i,poly in enumerate(polymorph):
    # Allocate storage for simulation data
    N_k = numpy.zeros([Kbig],numpy.int32) # N_k[k] is the total number of snapshots from alchemical state k

    u_kln = numpy.zeros([K,Kbig,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k 
            
    for k in range(K):
    	n = 0
	linenum = 0
	#dirname='../finishedJobs/' + directory + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + str(int(Temp))+'K' + Chargehinge + '_000L_100G_' + Pname + '_' + Volume_names[k]+hinge
	dirname='/oldhome/ecd4bd/finishedJobs_archive/' + directory + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + str(int(Temp))+'K' + Chargehinge + '_000L_100G_' + Pname + '_' + Volume_names[k]+hinge
        # cycle through all the pressure data
	fname=dirname+'/pressure.xvg'
	infile = open(fname, 'r')
	lines = infile.readlines()
	infile.close()
	print("loading " + fname)
	ignorecounter=0
	counter=0
        for line in lines:
	    tokens_pressure = line.split()
            if tokens_pressure[0] in ignore_symbols:
                continue
            #ignore the first set of frames
            if ignorecounter < ignoreframes:
                ignorecounter+=1
                continue
	    P[i,k] = P[i,k] * (float(counter)/(counter+1)) +  float(tokens_pressure[1])*(float(1)/(counter+1)) #Moving average of the pressure
	    counter+=1
        fname=dirname+'/'+gname
	print("loading " + fname)
	V[k] = numpy.round(calculate_gro_volume.Volume(fname),3)

        
    #============================================================================================= 
    # INTEGRATE POINTS USING SIMPSON INTEGRATION
    #=============================================================================================
    dV1 = numpy.round(V[1]-V[0],3)
    dV2 = numpy.round(V[2]-V[1],3)
    nm_to_M = 1.0e-09	#Conversion from angstroms into meters
    Bar_to_Pa = 100000	#Conversion from bar to pascals
    J_to_kcal = 0.2390057*0.001 #Conversion from kJ to kcal 
    Na = 6.022*10**23	#Avogadros numbers    

    if len(V)%2 != 1:
	print("Intermediate points must be an even number for simpsons rule to work!")
	sys.exit()
    if dV1 != dV2:
	print("Intermediate points must be evenly spaced for simpsons rule to work!")
	print("V1: " + str(dV1))
        print("V2: " + str(dV2))
        sys.exit()
    simp = numpy.zeros((K-1)/2,float)
    V_xaxis = numpy.zeros((K+1)/2,float)
    V_xaxis[0]=V[0]
    for k in range((K-1)/2):
	if k==0:
	    simp[k] = dV1/3*(P[i,k]+4*P[i,k+1]+P[i,k+2])
        else:
	    simp[k] = simp[k-1]+dV1/3*(P[i,k]+4*P[i,k+1]+P[i,k+2])
	dA[i,k+1] = simp[k]*(nm_to_M)**3 * Bar_to_Pa * J_to_kcal * Na / Molecules
	ddA[i,k+1] = 0 #How do we calculate error bars with the simpson rule? Reweighting to get derivative?
        V_xaxis[k+1]=V[2*k+2]

    print("Free Energies Optained...")
    
    #=============================================================================================
    # PLOT THE FINAL DATA
    #=============================================================================================
    """
    if (options.plot):    
        # now plot the pressure as a function of volume
        plt.figure(2*i)
        xlabel = 'Volume'
        ylabel = 'Average Pressure'
        plt.title(poly)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        print 'Xaxis:'
        print V
        print 'YAxis:'
        print P[i,:]
        plt.errorbar(V,P[i,:],numpy.zeros(len(V),float))
        #filename = poly + '_' + str(Molecules) + '_' + str(int(Temp)) + 'K_'+Gname+'_PvsV.png'
        plt.show()
	#plt.savefig(filename, bbox_inches='tight')

    	
	# Also plot the free energy vs V
	plt.figure(2*i+1)
        xlabel = 'Volume'
        ylabel = 'Free Energy (kcal/mol)'
        plt.title(poly)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        Xaxis = V_xaxis
        print 'Xaxis:'
        print Xaxis
        print 'YAxis:'
        print dA[i,:]
        plt.errorbar(Xaxis,dA[i,:])
        filename = poly + '_' + str(Molecules) + '_' + str(int(Temp)) + 'K_'+Gname+'_dAvsV.png'
        plt.savefig(filename, bbox_inches='tight')
    """
# now plot all 3 lines of the pressure as a function of volume
if (options.plot) and options.polymorphs == 'all':
    fig=plt.figure(7)
    ax=fig.add_subplot(111)
    xlabel = 'Volume (nm^3)'
    ylabel = 'Average Pressure (Bar)'
    #plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    print('Xaxis:')
    print(V)
    print('YAxis:')
    print(P[i,:])
    ax.plot(V,P[0,:], color='b', label='Polymorph 1')
    ax.plot(V,P[0,:], color='g', label='Polymorph 2')
    ax.plot(V,P[0,:], color='r', label='Polymorph 3')
    plt.plot(V,P[0,:],'b',V,P[1,:],'g',V,P[2,:],'r')
    filename =  'ALL_' + str(Molecules) + '_' + str(int(Temp)) + 'K_'+gname+'_PvsV.png'
    plt.legend(loc='upper right')
    plt.savefig(filename, bbox_inches='tight')
    plt.show()

    # Also plot all three lines of the free energy vs V
    plt.figure(8)
    xlabel = 'Volume (nm^3)'
    ylabel = 'Free Energy (kcal/mol)'
    #plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    Xaxis = V_xaxis
    print('Xaxis:')
    print(Xaxis)
    print('YAxis:')
    print(dA)
    plt.plot(Xaxis,dA[0,:],'b',Xaxis,dA[1,:],'g',Xaxis,dA[2,:],'r')
    #filename = 'All_' + str(Molecules) + '_' + str(int(Temp)) + 'K_'+Gname+'_dAvsV.png'
    plt.show()
    #plt.savefig(filename, bbox_inches='tight')

"""
#Plot the ideal gas value as well
    plt.figure(9)
    xlabel = 'Volume'
    ylabel = 'Average Pressure (Bar)'
    plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    print 'Xaxis:'
    print V
    print 'YAxis:'
    P[0,:] = Molecules*(1.38*10**(-23.0))*float(Temp)*(nm_to_M**(-3))/(Bar_to_Pa)/V
    print P[0,:]
    plt.plot(V,P[0,:])
    filename =  'ALL_' + str(Molecules) + '_' + str(int(Temp)) + 'K_Ideal_PvsV.png'
    plt.savefig(filename, bbox_inches='tight')

    dA = numpy.zeros(len(V));
    for i in range(len(V)):
	dA[i] = -(8.31)*float(Temp)*numpy.log(V[i]/V[0])*J_to_kcal
    plt.figure(10)
    xlabel = 'Volume'
    ylabel = 'Free Energy (kcal/mol)'
    plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    print 'Xaxis:'
    print V
    print 'YAxis:'
    print dA
    plt.plot(V,dA)
    filename = 'All_' + str(Molecules) + '_' + str(int(Temp)) + 'K_Ideal_dAvsV.png'
    #plt.savefig(filename, bbox_inches='tight')
"""


