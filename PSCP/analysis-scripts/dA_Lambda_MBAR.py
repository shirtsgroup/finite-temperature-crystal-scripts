from __future__ import print_function
#
# Computing the free energy difference of an organic crystal polymorph at different lambda values
# 
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy as np
import pymbar # multistate Bennett acceptance ratio
import MBARBootstrap # Bootstrapping algorithm
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import MBARBootstrap
import os.path
import sys
import panedr
import pdb
import tossconfigurationsFunc

def dA_Lambda_MBAR(plot_out=True, MinL=0, MaxL=100, dL=5, GAMMA=100, exponent=4, Molecule='benzene', polymorphs='p1 p2', 
                   Molecules=72, Independent=4, Temp=200, Pressure=1, k=1000, ignoreframes=500, includeframes=100000,
                   potential='oplsaa', hinge='DefaultHinge'):
    if (plot_out):
        import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from matplotlib.font_manager import FontProperties as FP
        font = {'family': 'normal',
                'weight': 'normal',
                'size': 16}
        matplotlib.rc('font', **font)
    
    
    # =============================================================================================
    # ENSURE THAT USER INPUTS ARE SENSIBLE
    # =============================================================================================
    # TEMPERATURE
    if Temp < 0:
        print("Invalid Temperature: " + str(Temp))
        sys.exit()
    
    if Pressure < 0:
        print("Invalid Pressure: " + str(Pressure))
        sys.exit()
    
    #LAMBDA
    if (MinL == -1) and (MaxL == -1) and (dL == -1) and (exponent == 1):
        print("Using default values!")
    
        # The Lambda points sampled
        Lambdas = ['000L', '010L', '020L', '030L', '040L', '050L', '060L', '070L', '080L', '090L', '100L']
    elif MinL < 0 or MaxL < 0 or dL < 0 or MinL > MaxL:
        print("Invalid Lambda Specifications")
        sys.exit()
    else:
        RawLambda = 0
        Lambdas = []
        lambda_names = np.arange(MinL, MaxL + dL, dL)

        Lambda_names = []
        Lambda_indicies=[]
        index = 0
        while RawLambda < MaxL:
            if RawLambda >= MinL:
                Lambda_indicies.append(index)
                index += 1
            else:
                index += 1
                RawLambda = RawLambda + dL
                continue
            Lambda = int(100 * float(RawLambda ** exponent) / float(MaxL ** exponent))
            Lambdas.append(Lambda)
            # Format the lambda point name
            if RawLambda < 10:
                Lambda_names.append('00' + str(int(RawLambda)) + 'L')
            elif RawLambda < 100:
                Lambda_names.append('0' + str(int(RawLambda)) + 'L')
            else:
                Lambda_names.append('100L')
            RawLambda=RawLambda+dL
        # Catch the final lambda point
        Lambdas.append(MaxL)
        Lambda_indicies.append(index)
        if MaxL < 10:
            Lambda_names.append('00' + str(int(MaxL)) + 'L')
        elif MaxL < 100:
            Lambda_names.append('0' + str(int(MaxL)) + 'L')
        else:
            Lambda_names.append('100L')
    
    # GAMMA
    if GAMMA < 0 or GAMMA > 100:
        print("Invalid Gamma Point: " + str(GAMMA))
        sys.exit()
    
    # POLYMORPH
    polymorphs = polymorphs.split()
    polymorph = []
    polymorph_short = []
    for i, token in enumerate(polymorphs):
        polymorph.append('Polymorph ' + str(token))
        polymorph_short.append(token)
    
    # POTENTIAL
    if potential not in ["oplsaa", "gromos", "designeda", "oplsaafakeg", "oplsaafakea"]:
        print("Invalid Potential")
        print("Supported potentials: oplsaa gromos designeda oplsaafakeg oplsaafakea")
        sys.exit()
    
    
    # =============================================================================================
    # FORMAT INPUTS
    # =============================================================================================
    # TEMPERATURE
#    Tname = ""
#    if Temp < 10:
#        Tname = "00" + str(int(Temp)) + "K"
#    elif Temp < 100:
#        Tname = "0" + str(int(Temp)) + "K"
#    else:
#        Tname = str(int(Temp)) + "K"
    
    # PRESSURE
#    Pname = ""
#    if Pressure < 10:
#        Pname = "00" + str(int(Pressure)) + "P"
#    elif Pressure < 100:
#        Pname = "0" + str(int(Pressure)) + "P"
#    else:
#        Pname = str(int(Pressure)) + "P"
    
    # GAMMA POINT
#    Gname = ""
#    if GAMMA < 10:
#        Gname = "00" + str(int(GAMMA)) + "G"
#    elif GAMMA < 100:
#        Gname = "0" + str(int(GAMMA)) + "G"
#    else:
#        Gname = str(int(GAMMA)) + "G"
    
    # NUMBER OF MOLECULES
#    Molname = ""
#    if Molecules == Independent:
#        Molname = str(Molecules) + '_'
#    else:
#        Molname = str(Molecules) + '_' + str(Independent) + 'ind_'
    
    
    # POTENTIAL
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
    
    # CHARGE AND SIGMA HINGE
#    if potential == "oplsaa":
#        ChargeHinge = ""
#    elif potential == "gromos":
#        ChargeHinge = ""
#    elif potential == "designeda":
#        ChargeHinge = ""
#    elif potential == "oplsaafakeg":
#        ChargeHinge = "_C01150"
#    elif potential == "oplsaafakea":
#        ChargeHinge = "_C01150"
    
    # OPTIONAL HINGE
    if str(GAMMA) == "100":
        hingeLetter = "L"
    else:
        hingeLetter = "R"
    
    if hinge == "DefaultHinge":
        hinges = ["_" + hingeLetter]
    else:
        # Read in each job
        hinges = []
        hingevect = hinge.split()
        for i, token in enumerate(hingevect):
            hinges.append("_" + hingeLetter + "_" + str(token))
    
    
    # =============================================================================================
    # READ IN RAW DATA
    # =============================================================================================
    # Constants.
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)  # Boltzmann constant in kcal/mol
    
    omitK = []
    
    # Parameters
    T_k = Temp*np.ones(len(Lambdas), float)  # Convert temperatures to floats
    
    g_k = np.zeros([len(Lambdas)], float)
    K = len(Lambdas)  # How many states?
    
    # total number of states examined; none are unsampled
    Kbig = K + 0
    
    # maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
    N_max = 200000
    
    # beta factor for the different temperatures
    beta_k = 1.0 / (kB * T_k)
    dA = np.zeros([len(polymorph), len(Lambdas)], float)
    ddA = np.zeros([len(polymorph), len(Lambdas)], float)
    convert_units = (0.2390057) * np.ones(len(Lambdas), float)  # Convert all energies to kcal/mol

    # Lines to ignore when reading in energies
    ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================']
    
    for i, poly in enumerate(polymorph):
        # Allocate storage for simulation data
        # N_k[k] is the total number of snapshots from alchemical state k
        N_k = np.zeros([Kbig], np.int32)
    
        # N_k_s[k,s] is the total number of snapshots from alchemical state k from seed s in 'unflipped segment j'
        N_ksj = np.zeros([Kbig, len(hinges), 100], np.int32)
    
        # u_kln[k,l,n] is the adjusted energy of snapshot n from simulation k
        u_kln = np.zeros([K, Kbig, N_max], np.float64)
    
        # dhdl_kln[k,l,n] is the restraint energy value of snapshop n from simulation k
        dhdl_kln = np.zeros([K, Kbig, N_max], np.float64)
    
        # dhdl_kn[k,n] is the derivative of energy with respect to lambda of snapshot n from simulation k
        dhdl_kn = np.zeros([K, N_max], np.float64)
    
        # Load in the data for each run
        for k in range(K):
            n = 0
            for s, hinge in enumerate(hinges):
#                tossconfigs = []  # The index of each configuration to toss from the MBAR analysis
                keepconfigs = np.arange(N_max)  # The index of each configuration to keep in the MBAR analysis
#                linenum_energy = 0
                linenum_dhdl = 0
                # cycle through all the input total energy data
                # dirpath='../finishedJobs/' + poly + '/benzene_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + '_' + Lambda_names[k] + '_' + Gname + '_' + Pname + hinge
                #dirpath = Molecule + '_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + \
                #          '_' + Lambda_names[k] + '_' + Gname + '_' + Pname + hinge
                dirpath = polymorph_short[i] + '/restraints/' + str(lambda_names[k])
                fname = dirpath + '/PROD.edr'
                dhdlname = dirpath + '/dhdl_PROD.xvg'
#                groname = dirpath + '/Traj.gro'
#                outname = dirpath + '/Flipless.gro'
#                restname = dirpath + '/restraint.gro'
                if k not in omitK:
#                    infile = open(fname, 'r')
#                    lines = infile.readlines()
#                    infile.close()
                    potential_energy = panedr.edr_to_df(fname)['Potential'].values
                    print("loading " + fname)

                    dhdl_energy = np.loadtxt(dhdlname, comments=['#', '$', '@', '!'])
                    print("loading " + dhdlname)

                    # Removing any non-equilibrated points of the simulation
                    [start_production, _, _] = timeseries.detectEquilibration(potential_energy)
                    potential_energy = potential_energy[start_production:]
                    dhdl_energy = dhdl_energy[start_production:]


                    # the energy of every configuration from each state evaluated at its sampled state
                    n = len(potential_energy)
                    u_kln[k, :, :n] = (float(Independent) / Molecules) * (potential_energy.reshape((n, 1)) +
                                                                          dhdl_energy[:, 5:]) * convert_units[k]
                    dhdl_kln[k, :, :n] = dhdl_energy[:, 5:] * convert_units[k]
                    dhdl_kn[k, :n] = (float(Independent) / Molecules) * dhdl_energy[:, 4] * convert_units[k]




                    ignorecounter=0
                    symbolcounter=0

#                    u_kln_hold = np.zeros([Kbig, N_max], np.float64)
#                    dhdl_kln_hold = np.zeros([Kbig, N_max], np.float64)
#                    dhdl_kn_hold = np.zeros([N_max], np.float64)

#                    for counter, line in enumerate(lines):
#                        tokens_energy = line.split()
#                        if tokens_energy[0] in ignore_symbols:
#                            symbolcounter += 1
#                            continue
#
#                        # ignore the first set of frames
#                        if ignorecounter <= ignoreframes:
#                            ignorecounter += 1
#                            continue
#
#                        # ignore the frames after the include frames
#                        if counter > includeframes:
#                            continue
#
#                        # Grab the dhdl information (if possible)
#                        tokens_dhdl = lines_dhdl[linenum_dhdl].split()
#                        while tokens_dhdl[0] in ignore_symbols:
#                            linenum_dhdl += 1
#                            tokens_dhdl = lines_dhdl[linenum_dhdl].split()
#
#                        while float(tokens_energy[0]) != float(tokens_dhdl[0]) and (linenum_dhdl+1) < len(lines_dhdl) \
#                                and linenum_dhdl < 1000000:
#                            linenum_dhdl += 1
#                            tokens_dhdl = lines_dhdl[linenum_dhdl].split()
#
#                        if float(tokens_energy[0]) != float(tokens_dhdl[0]):
#                            continue
#
#                        # the energy of every configuration from each state evaluated at its sampled state
#                        u_kln[k, :, n] = (float(Independent) / Molecules) * \
#                                         (float(tokens_energy[1]) +
#                                          np.asarray(np.array(tokens_dhdl)[5 + np.array(Lambda_indicies)], float)) \
#                                          * convert_units[k]
#                        dhdl_kln[k, :, n] = np.asarray(np.array(tokens_dhdl)[5 + np.array(Lambda_indicies)], float) \
#                                            * convert_units[k]
#                        dhdl_kn[k, n] = (float(Independent) / Molecules) * float(tokens_dhdl[4]) * convert_units[k]
#                        n += 1
                    

                    # Truncate the kept configuration list to be less than n
                    keepconfigs = [j for j in keepconfigs if j < (len(potential_energy)-symbolcounter) and j >= 0]
    
                    # Split up the retained configurations into connected segments
                    j = 0
                    a = 0
                    for a in range(len(keepconfigs)):
                        if a == 0:
                            continue
                        elif int(keepconfigs[a - 1]) + 1 != int(keepconfigs[a]):
                            N_ksj[k, s, j] = a - (sum(N_ksj[k, s, 0: j]))
                            j += 1
                    # Catch the final segment
                    N_ksj[k, s, j] = len(keepconfigs) - sum(N_ksj[k, s, 0: j])
                    j += 1

            N_k[k] = n
    
        # convert to nondimensional units from kcal/mol
        u_kln *= beta_k[0]
        # all data loaded from the three sets
        u_kln_save = u_kln.copy()
        g_k = np.zeros([K])
    
        # Ignore the first state due to jumping
        print("Number of retained samples")
        print(N_k)
    
        # =============================================================================================
        # COMPUTE FREE ENERGY DIFFERENCE USING MBAR
        # =============================================================================================
    
        # Initialize MBAR.
        print("Running MBAR...")
    
        # generate the weights of each of the umbrella set
        mbar = pymbar.MBAR(u_kln, N_k, verbose=True, subsampling_protocol=[{'method': 'L-BFGS-B'}])
    
        print("MBAR Converged...")
        # testing
    
        for k in range(Kbig):
            w = np.exp(mbar.Log_W_nk[:, k])
            print("max weight in state %d is %12.7f" % (k, np.max(w)))
            neff = 1 / np.sum(w ** 2)
            print("Effective number of sample in state %d is %10.3f" % (k, neff))
            print("Efficiency for state %d is %d/%d = %10.4f" % (k, neff, len(w), neff / len(w)))
    
        # extract self-consistent weights and uncertainties
        (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()
    
        print("Free Energies Optained...")
    
        #convert PMF to kcal/mol and normalize by the number of molecules
        df_i /= (beta_k[0] * float(Independent))
        ddf_i /= (beta_k[0] * float(Independent))
    
        dA[i, :] = df_i[0]


        # =============================================================================================
        # COMPUTE UNCERTAINTY USING THE UNCORRELATED DATA
        # =============================================================================================
    
        for k in range(K):  # For each restraint state
            N_k[k] = 0
            n_old = 0
            if k not in omitK:
                for s in range(len(hinges)):  # For each independent trajectory of this restraint state
                    for j in range(100):  # For each untossed segment of each independent trajectory of this restraint state
                        if N_ksj[k, s, j] == 0:
                            continue
                        #g_k[k] = timeseries.statisticalInefficiency(u_kln[k,k,0:N_k[k]])
                        # Feed in the segment and calculate correlation time
                        g_k[k] = timeseries.statisticalInefficiency(dhdl_kn[k, n_old:(n_old + N_ksj[k, s, j])])
                        print("Correlation time for sampled state %d is %10.3f" % (k, g_k[k]))
                        # subsample the data to get statistically uncorrelated data
                        # subsample indices within the segment
                        indices = np.array(timeseries.subsampleCorrelatedData(u_kln[k, k, n_old:(n_old + N_ksj[k, s, j])],
                                                                              g=g_k[k]))
                        # Apphend the uncorrelated configurations in the segment to the u_kln matrix
                        u_kln[k, :, N_k[k]: (N_k[k] + len(indices))] = u_kln_save[k, :, (indices+n_old)].transpose()
                        N_k[k] = N_k[k] + len(indices)
                        n_old += N_ksj[k, s, j]
    
        print("Number of retained samples")
        print(N_k)
        print("Number of retained samples from each seed")
        print(N_ksj)
    
        # generate the weights of each of the umbrella set
        mbar = pymbar.MBAR(u_kln, N_k, verbose=True, subsampling_protocol=[{'method': 'L-BFGS-B'}])
    
        print("MBAR Converged...")
        # testing
    
        # extract self-consistent weights and uncertainties
        (df_u, ddf_u, theta_i) = mbar.getFreeEnergyDifferences()
    
        print("Free Energies Optained...")
    
        # convert PMF to kcal/mol and normalize by the number of molecules
        df_u /= (beta_k[0] * float(Independent))
        ddf_u /= (beta_k[0] * float(Independent))
    
        ddA[i, :] = ddf_u[0]
    
        # Write out free energy differences
        print("Free Energy Difference (in units of kcal/mol)")
        for k in range(Kbig):
            print("%8.3f %8.3f" % (-df_i[k, 0], ddf_u[k, 0]))

    
    # =============================================================================================
    # PRINT THE FINAL DATA
    # =============================================================================================

    out_dA = np.zeros(len(polymorph))
    out_ddA = np.zeros(len(polymorph))
    for i, poly in enumerate(polymorph):
         out_dA[i] = -dA[i, Kbig-1]
         out_ddA[i] = ddA[i, Kbig - 1]

    
    # =============================================================================================
    # PLOT THE FINAL DATA
    # =============================================================================================

    if (plot_out) and polymorphs == 'all':
        # now plot the free energy change as a function of temperature
        fig = plt.figure(4)
        ax = fig.add_subplot(111)
        xlabel = 'Restraint Strength, $\lambda$'
        ylabel = 'Relative Free Energy (kcal/mol)'
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        Xaxis = [float(j / 100.0) for j in Lambdas]
        
        if os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph1_' + str(Molecules) + '_' + Tname + '_' + Pname +
                                  '_dAvsL_All'):
            ddA[0, :] = MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph1_' + str(Molecules) + '_' +
                                                       Tname + '_' + Pname + '_dAvsL_All')
        elif os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname +
                                    '_dAvsL_All'):
            ddA[1, :] = MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' +
                                                       Tname + '_' + Pname + '_dAvsL_All')
        elif os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname +
                                    '_dAvsL_All'):
            ddA[2, :] = MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph3_' + str(Molecules) + '_' +
                                                       Tname + '_' + Pname + '_dAvsL_All')
     
        ax.errorbar(Xaxis, dA[0, :], color='b', yerr=ddA[0, :], label='Benzene I')
        ax.errorbar(Xaxis, dA[1, :], color='g', yerr=ddA[1, :], label='Benzene II')
        ax.errorbar(Xaxis, dA[2, :], color='r', yerr=ddA[2, :], label='Benzene III')
        plt.legend(loc='upper left')

        if len(hinges) > 1:
            filename = PotNAME + '_' + str(Molecules) + '_' + Tname + '_dAvsL.pdf'
        else:
            filename = PotNAME + '_' + str(Molecules) + '_' + Tname + hinge + '_dAvsL.pdf'
        plt.show()
    return out_dA, out_ddA

if __name__ == '__main__':
    #=============================================================================================
    # READ IN USER INPUTS
    #=============================================================================================
    parser = OptionParser()
    parser.add_option('-p', '--plot', dest='plot', help='Plot output (default false)', default=True, action='store_true')
    parser.add_option('-U', dest='MinLambda', help='Minimum lambda point sampled', default=0)
    parser.add_option('-R', dest='MaxLambda', help='Maximum lambda point sampled', default=100)
    parser.add_option('-L', dest='spacing', help='Spacing between lambda points', default=5)
    parser.add_option('-G', dest='Gamma', help='Gamma point during the alchemical transformation', default=100)
    parser.add_option('-f', dest='exponent', help='functional form (exponent) of the spacing between the lambdas',
                      default=4)
    parser.add_option('-M', dest='molecule', help='name of the molecule', default='benzene')
    parser.add_option('-n', dest='polymorphs', help='Polymorphs to analyze', default='p1 p2')
    parser.add_option('-N', dest='molecules', help='number of supercell molecules', default=72)
    parser.add_option('-I', dest='independent', help='number of independent molecules', default=4)
    parser.add_option('-T', dest='Temperature', help='Temperature', default=200)
    parser.add_option('-P', dest='Pressure', help='Pressure', default=1)
    parser.add_option('-k', dest='forceConstant', help='Harmonic Force Constant', default=1000)
    parser.add_option('-i', dest='ignoreframes', help='Initial frames to ignore', default=2000)
    parser.add_option('-j', dest='includeframes', help='Number of frames to include', default=100000)
    parser.add_option('-u', dest='potential', help='Potential used in the simulation', default='oplsaa')
    parser.add_option('-H', '--hinge', dest='hinge', help='Optional string at end of jobs', default='DefaultHinge')
    
    (options, args) = parser.parse_args()
    plot_out = options.plot
    MinL = float(options.MinLambda)
    MaxL = float(options.MaxLambda)
    dL = float(options.spacing)
    GAMMA=int(options.Gamma)
    exponent = int(options.exponent)
    Molecule = options.molecule
    polymorphs = options.polymorphs
    Molecules = int(options.molecules)
    Independent = int(options.independent)
    Temp = float(options.Temperature)
    Pressure = int(options.Pressure)
    k = float(options.forceConstant)
    ignoreframes = int(options.ignoreframes)
    includeframes = int(options.includeframes)
    potential = str(options.potential)
    hinge = options.hinge
    
    dA, ddA = dA_Lambda_MBAR(plot_out=plot_out, MinL=MinL, MaxL=MaxL, dL=dL, GAMMA=GAMMA, exponent=exponent, 
                             Molecule=Molecule, polymorphs=polymorphs, Molecules=Molecules, Independent=Independent, 
                             Temp=Temp, Pressure=Pressure, k=k, ignoreframes=ignoreframes, includeframes=includeframes, 
                             potential=potential, hinge=hinge)
    for i, poly in enumerate(polymorphs):
        print(poly + ": " + "%8.3f %8.3f" % (dA[i], ddA[i]))
