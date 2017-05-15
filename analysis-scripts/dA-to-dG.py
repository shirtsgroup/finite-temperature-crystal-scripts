#
# Computing the Gibbs free energy difference between a set of polymorphs by integrating the helmholtz free energy at all volumes
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
parser.add_option('-S', dest = 'MinVolume', help = 'Minimum box vector scaling', default = 770)
parser.add_option('-B', dest = 'MaxVolume', help = 'Maximum box vector scaling', default = 970)
parser.add_option('-s', dest = 'spacing', help = 'Spacing between lambda points', default = 5)
parser.add_option('-n', dest = 'polymorphs', help = 'Polymorphs to analyze', default = 'all')
parser.add_option('-N', dest = 'molecules', help = 'number of supercell molecules', default = 72)
parser.add_option('-I', dest = 'independent', help = 'number of independent molecules', default = 4)
parser.add_option('-i', dest = 'ignoreframes', help = 'Initial frames to ignore', default = 2000)
parser.add_option('-j', dest = 'includeframes', help = 'Number of frames to include', default = 100000)
parser.add_option('-u', dest = 'potential', help = 'Potential used in the simulation', default = 'oplsaa')
parser.add_option('-T', dest = 'Temperature', help = 'Temperature', default = 200)
parser.add_option('-P', dest = 'Pressure', help = 'Pressure', default = 1)
#parser.add_option('-V', dest = 'Volume', help = 'Reference volume for the helmholtz free energy calculations', default = "8.755 8.679 8.667") #Gromos 72ind
#parser.add_option('-A', dest = 'Helmholtz', help = 'Helmholtz free energy for each polymorph at the reference volume', default = "-11.014 -10.796 -10.835") #Gromos 72ind
#parser.add_option('-V', dest = 'Volume', help = 'Reference volume for the helmholtz free energy calculations', default = "8.773 8.651 8.645") #OPLSAAFAKEG 72ind
#parser.add_option('-A', dest = 'Helmholtz', help = 'Helmholtz free energy for each polymorph at the reference volume', default = "-11.378 -11.213 -11.248") #OPLSAAFAKEG 72ind
#parser.add_option('-V', dest = 'Volume', help = 'Reference volume for the helmholtz free energy calculations', default = "8.755 8.679 8.667") #Gromos 72ind
#parser.add_option('-A', dest = 'Helmholtz', help = 'Helmholtz free energy for each polymorph at the reference volume', default = "-11.014 -10.796 -10.835") #Gromos 72ind
#parser.add_option('-V', dest = 'Volume', help = 'Reference volume for the helmholtz free energy calculations', default = "0.484 0.477 0.481") #GROMOS 4ind
#parser.add_option('-A', dest = 'Helmholtz', help = 'Helmholtz free energy for each polymorph at the reference volume', default = "-11.347 -11.201 -10.953") #OPLSAAFAKEG 4ind
parser.add_option('-d', dest = 'directory', help = 'Parent directory of the volume change directories', default = 'Volume')
parser.add_option('-H', '--hinge', dest = 'hinge', help = 'Optional string at end of jobs', default = 'DefaultHinge')

(options, args) = parser.parse_args()
MinV = float(options.MinVolume)
MaxV = float(options.MaxVolume)
dV = float(options.spacing)
Temp = float(options.Temperature)
Pressure = int(options.Pressure)
Molecules = int(options.molecules)
Independent = int(options.independent)
ignoreframes = int(options.ignoreframes)
includeframes = int(options.includeframes)
potential = str(options.potential)
directory = options.directory
hinge = options.hinge

AvgVolume = dict()
Helmholtz = dict()

AvgVolume['gromos'] = [0.480, 0.475, 0.479]
AvgVolume['oplsaafakeg'] = [0.481, 0.474, 0.476]
AvgVolume['oplsaa'] = [0.486, 0.480, 0.478]
AvgVolume['oplsaafakea'] = [0.494, 0.474, 0.478]

#200K
Helmholtz['gromos'] = [-10.881, -10.626, -10.529]
Helmholtz['oplsaafakeg'] = [-11.221, -11.038, -10.918]
Helmholtz['oplsaa'] = [-11.077, -10.888, -10.770]
Helmholtz['oplsaafakea'] = [-10.618, -10.985, -10.911]

#100K
#Helmholtz['gromos'] = [-10,145, -9.867, -9.814]
#Helmholtz['oplsaa'] = [-10.402, -10.176, -10.115]


#72 Benzene
#Helmholtz['gromos'] = [-11.077, -10.888, -10.770]
#Helmholtz['oplsaa'] = [-11.261, -11.058, -11.061]

if (options.plot):
    import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.font_manager import FontProperties as FP
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
    matplotlib.rc('font', **font)
    colors = ('b', 'g', 'r')
    markers= ('o', 'o', 'o')
#=============================================================================================
# ENSURE THAT USER INPUTS ARE SENSIBLE
#=============================================================================================

if Temp < 0:
    print "Invalid Temperature: " + str(Temp)
    sys.exit()

if (MinV == -1 ) and (MaxV == -1) and (dV == -1):
    print "Using default values!"
    Volumes = ['v100','v102','v104','v106','v108','v110','v112','v114','v116','v118','v120'] #The scaling parameters sampled
elif MinV < 0 or MaxV < 0 or dV < 0 or MinV > MaxV:
    print "Invalid Volume Specifications"
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
    polymorph = ['Benzene I', 'Benzene II', 'Benzene III']
    polymorph_short = ['p1', 'p2', 'p3']
elif (options.polymorphs == 'p1'):
    polymorph = ['Benzene I']
    polymorph_short = ['p1']
elif (options.polymorphs == 'p2'):
    polymorph = ['Benzene II']
    polymorph_short = ['p2']
elif (options.polymorphs == 'p3'):
    polymorph = ['Benzene III']
    polymorph_short = ['p3']
else:
    print "Polymorph Inputs Wrong"
    sys.exit()

#POTENTIAL
if potential != "oplsaa" and potential != "gromos" and potential != "oplsaafakeg" and potential != "oplsaafakea":
    print "Invalid Potential"
    print "Supported potentials: oplsaa gromos oplsaafakeg oplsaafakea"
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
    Chargehinge="_C01150"
elif potential == "oplsaafakea":
    PotNAME = "FAKEA"
    Chargehinge="_C01150"

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

gname = "pre_EQ.gro"
molecule = 'Benzene' #The crystalline molecule
ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================']; #Lines to ignore when reading in energies


#Read in the reference volume and helmholtz free energy
Vref = numpy.zeros(len(polymorph),float)
Aref = numpy.zeros(len(polymorph),float)
A_shift = 0.0
for i,token in enumerate(AvgVolume[potential]):
    Vref[i]=float(token)
for i,token in enumerate(Helmholtz[potential]):
    Aref[i]=int(Independent)*float(token) + A_shift


# Parameters
T_k = Temp*numpy.ones(len(Volumes),float) #Convert temperatures to floats
P = numpy.zeros([3,len(Volumes)],float) #Convert Pressures to floats
ddP = numpy.zeros([3,len(Volumes)],float) #Standard Deviation of the Mean Pressure at each volume
V = Temp*numpy.ones(len(Volumes),float) #Convert Volumes to floats
P_NVT = numpy.zeros([3,20000],float) #Pressure of each configuration in the NVT simulation
V_NPT = numpy.zeros([3,20000],float) #Volume of each configuration in the NPT simulation
N_k = numpy.zeros(3,int) #Number of configurations from the NPT simulation for each polymorph
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
    #N_k = numpy.zeros([Kbig],numpy.int32) # N_k[k] is the total number of snapshots from alchemical state k

    #u_kln = numpy.zeros([K,Kbig,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k 
            
    for k in range(K):
    	n = 0
	linenum = 0
	#dirname='../finishedJobs/' + directory + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + str(int(Temp))+'K' + Chargehinge + '_000L_100G_' + Pname
        dirname='/oldhome/ecd4bd/finishedJobs_archive/' + directory + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + str(int(Temp))+'K' + Chargehinge + '_000L_100G_' + Pname
        # cycle through all the pressure data
	fname=dirname + '_' + Volume_names[k]+hinge + '/pressure.xvg'
	infile = open(fname, 'r')
	lines = infile.readlines()
	infile.close()
	print "loading " + fname
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
	    P_NVT[i,counter] = float(tokens_pressure[1])
	    P[i,k] = P[i,k] * (float(counter)/(counter+1)) +  float(tokens_pressure[1])*(float(1)/(counter+1)) #Moving average of the pressure
	    counter+=1
	#Calculate the standard deviation of the pressure
	ddP[i,k] = numpy.std(P_NVT[i,:counter])/(counter**0.5)
	fname=dirname+ '_' + Volume_names[k]+hinge +'/'+gname
        print "loading " + fname
	V[k] = (float(Independent)/Molecules)*numpy.round(calculate_gro_volume.Volume(fname),3)

    ##Read in the actual volume distribution
    #fname_volume=dirname+'/volume.xvg'
    ##fname_volume = '../finishedJobs/' + directory + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + str(int(Temp))+'K' + Chargehinge + '_000L_100G_' + Pname + '_c/volume.xvg'
    #infile = open(fname_volume, 'r')
    #lines_volume = infile.readlines()
    #infile.close()
    #print "loading " + fname_volume
    #ignorecounter=0
    #n=0;
    #for j,line in enumerate(lines_volume):
    #    tokens_volume = line.split()
    #    if tokens_volume[0] in ignore_symbols:
    #        continue
    #    #ignore the first set of frames
    #    if ignorecounter < ignoreframes/10.0:
    #        ignorecounter+=1
    #        continue
    #    #time[n] = float(tokens_volume[0])/1000.0
    #    V_NPT[i,n] = float(tokens_volume[1]);
    #    n+=1
    #N_k[i]=n




        
#=====================================================================================================
# REGRESS A POLYNOMIAL FIT TO THE PV CURVES AND INTEGRATE TO FIND A(V) FOR EACH POLYMORPH
#=====================================================================================================
nm_to_M = 1.0e-09	#Conversion from angstroms into meters
Bar_to_Pa = 100000	#Conversion from bar to pascals
J_to_kcal = 0.2390057*0.001 #Conversion from kJ to kcal 
Na = 6.022*10**23	#Avogadros numbers    
degree = 4		#degree of the polynomial fit
dV=0.001		#Incremental volume for plotting the regressed fit
V_axis = numpy.arange(V[0],V[len(V)-1]+dV,dV) #V axis for plotting and numerically integrating
P_axis = numpy.zeros([len(polymorph),len(V_axis)],float)#Fitted pressure at each volume in V_axis for each polymorph
A = numpy.zeros([len(polymorph),len(V_axis)],float) 	#Helmholtz free energy estimate for each polymorph at each volume
Prob_V = numpy.zeros([len(polymorph),len(V_axis)],float)     #Probability of occupying a volume between V and V+dV for each polymorph
G = numpy.zeros([len(polymorph)],float)     #Gibbs free energy estimate for each polymorph (from interacting crystal to noninteracting ideal gas)

#Regress the polynomial fit
polyfit = numpy.transpose(numpy.polyfit(V,numpy.transpose(P),degree))
polyfit_integral = numpy.zeros([len(polymorph),degree+1],float)

for i in range(degree+1):
    polyfit_integral[:,i] = polyfit[:,i]*1.0/(degree-i+1)

#Calculate the regressed pressure and the helmholtz free energy
for i in range(len(polymorph)):
    A[i,:] = Aref[i]	#Add the reference free energy term
    for j in range(degree+1):
	P_axis[i,:] += polyfit[i,j]*numpy.power(V_axis[:],degree-j)
	A[i,:] -= (nm_to_M)**3 * Bar_to_Pa * J_to_kcal * Na *(polyfit_integral[i,j]*numpy.power(V_axis[:],degree-j+1) - polyfit_integral[i,j]*numpy.power(Vref[i],degree-j+1))

pdb.set_trace()

#Calculate the probability of being at each volume
for i in range(len(polymorph)):
    Prob_V[i,:] = numpy.exp(-1*beta_k[0]*A[i,:]) * numpy.exp(-((nm_to_M)**3 * Bar_to_Pa * J_to_kcal * Na)*beta_k[0]*Pressure*V_axis[:])*dV
    G[i] = -1.0/beta_k[0] * numpy.log(numpy.sum(Prob_V[i,:]))/float(Independent)
    Prob_V[i,:] /= numpy.sum(Prob_V[i,:])

#Calculate the Gibbs free energy difference between the different polymorphs
#for i in range(len(polymorph)):
#    G[i] = numpy.dot((A[i,:]-A_shift),Prob_V[i,:])# - numpy.dot(A[0,:],Prob_V[0,:]) 
#G /= float(Independent)

#Print the integrated free energy difference:
print "Helmholtz Reference:"
print (Aref-A_shift)/float(Independent)
print "PV Correction:"
print G-(Aref-A_shift)/float(Independent)
print "Gibbs Free Energy:"
print G
#=============================================================================================
# PLOT THE FINAL DATA
#=============================================================================================
# now plot all 3 lines of the pressure as a function of volume
if (options.plot) and options.polymorphs == 'all':
    fig=plt.figure(1)
    ax=fig.add_subplot(111)
    xlabel = 'Volume (nm^3)'
    ylabel = 'Average Pressure (Bar)'
    #plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #print 'Xaxis:'
    #print V
    #print 'YAxis:'
    #print P[i,:]
    for i in range(len(polymorph)):
        ax.errorbar(V,P[i,:],ddP[i,:],color=colors[i],marker=markers[i],linestyle='None') #Plot the actual data
	ax.plot(V_axis,P_axis[i,:],color=colors[i],label=polymorph[i]) #Plot the regressed fit
	#plt.hold(true)
    plt.legend(loc='upper right')
    filename='../../Pictures/' + PotNAME + '_' + str(Molecules) + '_' + str(Independent) + 'ind_' + str(int(Temp)) + 'K'+'_PvsV.pdf'
    plt.savefig(filename, bbox_inches='tight')
    plt.show()

    # Also plot all three lines of the free energy vs V
    fig=plt.figure(2)
    ax=fig.add_subplot(111)
    xlabel = 'Volume (nm^3)'
    ylabel = 'Free Energy (kcal/mol)'
    #plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #print V_axis
    #print A
    for i in range(len(polymorph)):
	ax.plot(V_axis,A[i,:],color=colors[i]) #Plot the free energy    
    #filename='All_' + PotNAME + '_' + str(Molecules) + '_' + str(int(Temp)) + 'K'+'_dAvsV.png'
    #plt.savefig(filename, bbox_inches='tight')
    plt.show()

    # Also plot all three histograms of volume
    fig=plt.figure(3)
    xlabel = 'Volume'
    ylabel = 'Probability'
    #plt.title('All Polymorphs')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #print V_axis
    #print Prob_V
    for i in range(len(polymorph)):
        plt.plot(V_axis,Prob_V[i,:]/(V_axis[len(V_axis)-1]-V_axis[0])/10,color=colors[i]) #Plot the free energy
	Volumes = V_NPT[i,0:N_k[i]-1]
	weights = numpy.ones_like(Volumes)/len(Volumes)
        plt.hist(Volumes,20,alpha=0.3,color=colors[i],label=polymorph[i],weights=weights)
    #filename = 'All_' + str(Molecules) + '_' + str(int(Temp)) + 'K'+'_dAvsV.png'
    plt.show()
    #plt.savefig(filename, bbox_inches='tight')
   


