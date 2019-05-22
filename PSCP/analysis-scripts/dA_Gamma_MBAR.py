from __future__ import print_function
#
# Computing the free energy difference of an organic crystal polymorph at different gamma values
# 
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#
import numpy as np
import pymbar  # multistate Bennett acceptance ratio
from pymbar import timeseries  # timeseries analysis
import MBARBootstrap  # Bootstrapping algorithm
import os.path
from optparse import OptionParser # for parsing command-line options
import sys
import Harvist #Hamiltonian Reweighting Visualization Toolkit
import pdb


def dA_Gamma_MBAR(plot_out=True, MINGAMMA=0, MAXGAMMA=100, GSPACING=10, LAMBDA=100, exponent=2, polymorphs='p1 p2', 
                  Molecule='benzene', Molecules=72, Independent=4, Temp=200, Pressure=1, k=1000, ignoreframes=100,
                  includeframes=100000, potential='oplsaa',bonds='no', hinge='DefaultHinge'):
    Colors = ['b', 'r', 'g']
    
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
    
    # GAMMA
    if (MINGAMMA == -1) and (MAXGAMMA == -1) and (GSPACING == -1) and (exponent == 1):
        print("Using default values!")
        # The Gamma points sampled
        Gammas = ['000L', '010L', '020L', '030L', '040L', '050L', '060L', '070L', '080L', '090L', '100L']
    elif MINGAMMA < 0 or MAXGAMMA < 0 or GSPACING < 0 or MINGAMMA > MAXGAMMA:
        print("Invalid Gamma Specifications")
        sys.exit()
    else:
        RawGamma = MINGAMMA
        Gammas = []
        Gamma_names = []
        gamma_names = np.arange(MINGAMMA, MAXGAMMA + GSPACING, GSPACING)
        while RawGamma < MAXGAMMA:
            Gamma = int(100 * float(RawGamma ** exponent) / float(MAXGAMMA ** exponent))
            Gammas.append(Gamma)
            # Format the gamma point name
            if RawGamma < 10:
                Gamma_names.append('00' + str(int(RawGamma)) + 'G')
            elif RawGamma < 100:
                Gamma_names.append('0' + str(int(RawGamma)) + 'G')
            else:
                Gamma_names.append('100G')
            RawGamma = RawGamma + GSPACING
        # Catch the final gamma point
        Gammas.append(int(MAXGAMMA))
        if MAXGAMMA < 10:
            Gamma_names.append('00' + str(int(MAXGAMMA)) + 'G')
        elif MAXGAMMA < 100:
            Gamma_names.append('0' + str(int(MAXGAMMA)) + 'G')
        else:
            Gamma_names.append('100G')
    
    # LAMBDA
    if LAMBDA < 0 or LAMBDA > 100:
        print("Invalid Lambda Point: " + str(LAMBDA))
        sys.exit()
    
    # POLYMORPH
    polymorphs = polymorphs.split()
    polymorph = []
    polymorph_short = []
    for i, token in enumerate(polymorphs):
        polymorph.append('Polymorph ' + str(token))
        polymorph_short.append(token)
    
    # POTENTIAL
    if potential != "oplsaa" and potential != "gromos" and potential != "designeda" and potential != "oplsaafakeg" and \
                    potential != "oplsaafakea":
        print("Invalid Potential")
        print("Supported potentials: oplsaa gromos designeda oplsaafakeg oplsaafakea")
        sys.exit()
    
    # =============================================================================================
    # FORMAT INPUTS
    # =============================================================================================
    # TEMPERATURE
    Tname = ""
    if Temp < 10:
        Tname = "00" + str(int(Temp)) + "K"
    elif Temp < 100:
        Tname = "0" + str(int(Temp)) + "K"
    else:
        Tname = str(int(Temp)) + "K"
    
    # PRESSURE
    Pname = ""
    if Pressure < 10:
        Pname = "00" + str(int(Pressure)) + "P"
    elif Pressure < 100:
        Pname = "0" + str(int(Pressure)) + "P"
    else:
        Pname = str(int(Pressure)) + "P"
    
    # LAMBDA
    Lname = ""
    if LAMBDA < 10:
        Lname = "00" + str(LAMBDA) + "L"
    elif LAMBDA < 100:
        Lname = "0" + str(LAMBDA) + "L"
    else:
        Lname = str(LAMBDA) + "L"
    
    # NUMBER OF MOLECULES
    Molname = ""
    if Molecules == Independent:
        Molname = str(Molecules) + '_'
    else:
        Molname = str(Molecules) + '_' + str(Independent) + 'ind_'
    
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
    if potential == "oplsaa":
        ChargeHinge = ""
    elif potential == "gromos":
        ChargeHinge = ""
    elif potential == "designeda":
        ChargeHinge = ""
    elif potential == "oplsaafakeg":
        ChargeHinge = "_C01150"
    elif potential == "oplsaafakea":
        ChargeHinge = "_C01150"
    
    # OPTIONAL HINGE
    if hinge == "DefaultHinge":
        hinges = ['_G']
    else:
        # Read in each job
        hinges = []
        hingevect = options.hinge.split()
        for i, token in enumerate(hingevect):
            hinges.append("_G_" + str(token))
    
    
    # =============================================================================================
    # READ IN RAW DATA
    # =============================================================================================
    # Constants.
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)  # Boltzmann constant in kcal/mol
    
    omitT = []  # Temperatures to be omitted from the analysis
    
    # Parameters
    T_k = Temp * np.ones(len(Gammas), float)  # Convert temperatures to floats
    print(T_k)
    print(Gammas)

    g_k = np.zeros([len(Gammas)], float)
    K = len(Gammas)  # How many states?
     
    # total number of states examined; 0 are unsampled if bonds are left on, 1 is unsampled if the bonds are removed
    if bonds != "no":
        Kbig = K + 1
    else:
        Kbig = K
    
    # maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
    N_max = 200000
    
    # beta factor for the different temperatures
    beta_k = 1.0 / (kB * T_k)
    dA = np.zeros([len(polymorph), Kbig], float)
    ddA = np.zeros([len(polymorph), Kbig], float)
    convert_units = 0.2390057 * np.ones(Kbig, float)  # Convert all energies to kcal/mol
    #convert_units = (1.0)*np.ones(len(Gammas),float)  # Convert all energies to kcal/mol
    ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================']  # Lines to ignore when reading in energies
    
    # Allocate storage for simulation data
    for i, poly in enumerate(polymorph):
        # N_k[k] is the total number of snapshots from alchemical state k
        N_k = np.zeros([Kbig], np.int32)
    
        # N_k_s[k,s] is the total number of snapshots from alchemical state k from seed s
        N_k_s = np.zeros([Kbig, len(hinges)], np.int32)
    
        # u_kln[k,l,n] is the adjusted energy of snapshot n from simulation k
        u_kln = np.zeros([K, Kbig, N_max], np.float64)
    
        # dhdl_kn[k,n] is the derivative of energy with respect to lambda of snapshot n from simulation k
        dhdl_kn = np.zeros([K, N_max], np.float64)
    
        #Load in the data for each run
        for k in range(K):
            n = 0
            for s, hinge in enumerate(hinges):
                linenum_energy = 0
                linenum_dhdl = 0
                # cycle through all the input total energy data
                dirpath = Molecule  + '_GRO_' + PotNAME + '_' + polymorph_short[i] + '_' + Molname + Tname + ChargeHinge + \
                          '_' + Lname + '_' + Gamma_names[k] + '_' + Pname + hinge
                dirpath = polymorph_short[i] + '/interactions/' + str(gamma_names[k])

                fname = dirpath + '/potenergy.xvg'
                dhdlname = dirpath + '/dhdl_PROD.xvg'
                fname14 = dirpath + '/potenergy14.xvg'
                groname = dirpath + '/Traj.gro'
                outname = dirpath + '/Flipless.gro'
                restname = dirpath + '/restraint.gro'
                if k not in omitT:
                    #if not os.path.isfile(fname):
                    #    continue
                    infile = open(fname, 'r')
                    lines = infile.readlines()
                    infile.close()
                    print("loading " + fname)
                    infile = open(dhdlname, 'r')
                    lines_dhdl = infile.readlines()
                    infile.close()
                    print("loading " + dhdlname)
                    ignorecounter = 0
                    for counter, line in enumerate(lines):
                        tokens_energy = line.split()
                        if tokens_energy[0] in ignore_symbols:
                            continue
                        # ignore the first set of frames
                        if ignorecounter <= ignoreframes:
                            ignorecounter += 1
                            continue
                        ## ignore the frames after the include frames
                        #if counter > includeframes:
                        #    continue
                        # Grab the dhdl information
                        tokens_dhdl = lines_dhdl[linenum_dhdl].split()
                        while tokens_dhdl[0] in ignore_symbols:
                            linenum_dhdl += 1
                            tokens_dhdl = lines_dhdl[linenum_dhdl].split()
                        while float(tokens_energy[0]) != float(tokens_dhdl[0]) and (linenum_dhdl + 1) < len(lines_dhdl) \
                                and linenum_dhdl < 1000000:
                            linenum_dhdl += 1
                            tokens_dhdl = lines_dhdl[linenum_dhdl].split()

                        # the energy of every configuration from each state evaluated at its sampled state
                        if float(tokens_energy[0]) != float(tokens_dhdl[0]):
                            continue
    
                        # Use this one if this is a normal PSCP
                        u_kln[k, :K, n] = (float(tokens_energy[1]) + np.asarray(tokens_dhdl[5:], float)) * convert_units[k]
                        dhdl_kn[k, n] = (float(Independent) / Molecules) * (float(tokens_dhdl[2]) + float(tokens_dhdl[3])
                                                                            + 0.0 * float(tokens_dhdl[4])) \
                                        * convert_units[k]
                        n += 1
                if s == 0:
                    N_k_s[k, s] = n
                else:
                    N_k_s[k, s] = n - sum(N_k_s[k, 0:s])
            N_k[k] = n

        # convert to nondimensional units from kcal/mol
    
        u_kln *= beta_k[0]
        # all data loaded from the three sets
    
        u_kln_save = u_kln.copy()
        N_k_save = N_k.copy()
        g_k = np.zeros([K])

        print("Number of retained samples")
        print(N_k)
        print("Number of retained samples from each seed")
        print(N_k_s)


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
            # using Kish (1965) formula.
            # effective # of samples =  (\sum_{i=1}^N w_i)^2 / \sum_{i=1}^N w_i^2
            #                        =  (\sum_{i=1}^N w_i^2)^-1
            neff = 1 / np.sum(w ** 2)
            print("Effective number of sample in state %d is %10.3f" % (k, neff))
            print("Efficiency for state %d is %d/%d = %10.4f" % (k, neff, len(w), neff / len(w)))
    
        # extract self-consistent weights and uncertainties
        (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()
    
        print("Free Energies Optained...")
    
        # convert PMF to kcal/mol and normalize by the number of molecules
        df_i /= (beta_k[0] * float(Independent))
        ddf_i /= (beta_k[0] * float(Independent))
    
        dA[i, :] = df_i[0]

    
        # =============================================================================================
        # COMPUTE UNCERTAINTY USING THE UNCORRELATED DATA
        # =============================================================================================
        
        for k in range(K):
            N_k[k] = 0
            n_old = 0
            if k not in omitT:
                for s in range(len(hinges)):
                    g_k[k] = timeseries.statisticalInefficiency(dhdl_kn[k, n_old: (n_old + N_k_s[k, s])])
                    print("Correlation time for sampled state %d is %10.3f" % (k, g_k[k]))
                    # subsample the data to get statistically uncorrelated data
                    indices = np.array(timeseries.subsampleCorrelatedData(u_kln[k, k, n_old:(n_old + N_k_s[k, s])],
                                                                          g=g_k[k]))  # subsample
    
                    # not sure why we have to transpose
                    u_kln[k, :, N_k[k]: (N_k[k] + len(indices))] = u_kln_save[k, :, (indices + n_old)].transpose()
                    N_k[k] = N_k[k] + len(indices)
                    n_old += N_k_s[k, s]
        print("Number of retained samples")
        print(N_k)
        print("Number of retained samples from each seed")
        print(N_k_s)
    
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
        out_dA[i] = dA[i, Kbig - 1]
        out_ddA[i] = ddA[i, Kbig - 1]


    # =============================================================================================
    # PLOT THE FINAL DATA
    # =============================================================================================
    
    if (plot_out) and polymorphs == 'all':
        # now plot the free energy change as a function of temperature
        fig = plt.figure(4)
        ax = fig.add_subplot(111)
        xlabel = 'Interaction Strength, $\gamma$'
        ylabel = 'Relative Free Energy (kcal/mol)'
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        Xaxis = [float(j / 100.0) for j in Gammas]
       
        if os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph1_' + str(Molecules) + '_' + Tname + '_' + Pname +
                                  '_dAvsG_All'):
            ddA[0, :] = MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph1_' + str(Molecules) + '_' +
                                                       Tname + '_' + Pname + '_dAvsG_All')
        elif os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname +
                                    '_dAvsG_All'):
            ddA[1, :] = MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' +
                                                       Tname + '_' + Pname + '_dAvsG_All')
        elif os.path.isfile('BootstrapStd_' + PotNAME + '_Polymorph2_' + str(Molecules) + '_' + Tname + '_' + Pname +
                                    '_dAvsG_All'):
            ddA[2, :] = MBARBootstrap.ExtractBootstrap('BootstrapStd_' + PotNAME + '_Polymorph3_' + str(Molecules) + '_' +
                                                       Tname + '_' + Pname + '_dAvsG_All')
    
        ax.errorbar(Xaxis, dA[0, :], color='b', yerr=ddA[0, :], label='Benzene I')
        ax.errorbar(Xaxis, dA[1, :], color='g', yerr=ddA[1, :], label='Benzene II')
        ax.errorbar(Xaxis, dA[2, :], color='r', yerr=ddA[2, :], label='Benzene III')
        plt.legend(loc='upper right')

        if len(hinges) > 1:
            filename = PotNAME + '_' + str(Molecules) + '_' + Tname + '_dAvsG.pdf'
        else:
            filename = PotNAME + '_' + str(Molecules) + '_' + Tname + hinge + '_dAvsG.pdf'
        plt.savefig(filename, bbox_inches='tight')
    return out_dA, out_ddA


if __name__ == '__main__':
    # =============================================================================================
    # READ IN USER INPUTS
    # =============================================================================================
    parser = OptionParser()
    parser.add_option('-p', '--plot', dest='plot', help='Plot output (default false)', default=True, action='store_true')
    parser.add_option('-A', dest='MaxGamma', help='Maximum gamma point sampled', default=100)
    parser.add_option('-B', dest='MinGamma', help='Minimum gamma point sampled', default=0)
    parser.add_option('-g', dest='spacing', help='Spacing between gamma points', default=10)
    parser.add_option('-L', dest='Lambda', help='Lambda point during the alchemical transformation', default=100)
    parser.add_option('-f', dest='exponent', help='functional form (exponent) of the spacing between the gammas', default=2)
    parser.add_option('-n', dest='polymorphs', help='Polymorphs to analyze', default='p1 p2')
    parser.add_option('-M', dest='molecule', help='name of the molecule', default='benzene')
    parser.add_option('-N', dest='molecules', help='number of supercell molecules', default=72)
    parser.add_option('-I', dest='independent', help='number of independent molecules', default=4)
    parser.add_option('-T', dest='Temperature', help='Temperature', default=200)
    parser.add_option('-P', dest='Pressure', help='Pressure', default=1)
    parser.add_option('-k', dest='forceConstant', help='Harmonic Force Constant', default=1000)
    parser.add_option('-i', dest='ignoreframes', help='Initial frames to ignore', default=2000)
    parser.add_option('-j', dest='includeframes', help='Number of frames to include', default=100000)
    parser.add_option('-u', dest='potential', help='Potential used in the simulation', default='oplsaa')
    parser.add_option('-b', dest='bonds', help='Bonds also removed (yes/no)', default='no')
    parser.add_option('-H', '--hinge', dest='hinge', help='Optional string at end of jobs', default='DefaultHinge')
    
    (options, args) = parser.parse_args()
    plot_out = options.plot
    MINGAMMA = float(options.MinGamma)
    MAXGAMMA = float(options.MaxGamma)
    GSPACING = float(options.spacing)
    LAMBDA=int(options.Lambda)
    exponent = int(options.exponent)
    polymorphs = options.polymorphs
    Molecule = options.molecule
    Molecules = int(options.molecules)
    Independent = int(options.independent)
    Temp = float(options.Temperature)
    Pressure = int(options.Pressure)
    k = float(options.forceConstant)
    ignoreframes = int(options.ignoreframes)
    includeframes = int(options.includeframes)
    potential = str(options.potential)
    bonds = str(options.bonds)
    hinge = options.hinge

    dA, ddA = dA_Gamma_MBAR(plot_out=plot_out, MINGAMMA=MINGAMMA, MAXGAMMA=MAXGAMMA, GSPACING=GSPACING, LAMBDA=LAMBDA, 
                            exponent=exponent, polymorphs=polymorphs, Molecule=Molecule, Molecules=Molecules, 
                            Independent=Independent, Temp=Temp, Pressure=Pressure, k=k, ignoreframes=ignoreframes, 
                            includeframes=includeframes, potential=potential,bonds=bonds, hinge=hinge)

    for i, poly in enumerate(polymorph):
        print(poly + ": " + "%8.3f %8.3f" % (dA, ddA))



