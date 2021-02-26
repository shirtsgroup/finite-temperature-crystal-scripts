from __future__ import print_function
#
# Computing the free energy difference of an organic crystal polymorph at different gamma values
# 
# Copyright Eric Dybeck and Michael R. Shirts, University of Virginia, 2014
#

import numpy as np
import pymbar  # multistate Bennett acceptance ratio
from pymbar import timeseries  # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import sys
import panedr
import os

def dA_MBAR(minimum=0, maximum=100, spacing=10, exponent=2, polymorphs='p1 p2', Molecules=72, Independent=4, Temp=200,
            bonds=False, primary_directory='.', added_directories=[]):
    # =============================================================================================
    # Setting up the values for gamma or lambda states
    # =============================================================================================
#    raw_value = minimum
#    values = []
    directory_names = np.arange(minimum, maximum + spacing, spacing)
    directory_names = np.sort(np.append(directory_names, added_directories)) 

#    while raw_value <= maximum:
#        if exponent >= 0:
#            value = int(100 * (float(raw_value) / float(maximum)) ** abs(exponent))
#        else:
#            value = int(100 * (1 - (float(maximum - raw_value) / float(maximum)) ** abs(exponent)))
#        values.append(value)
#        raw_value = raw_value + spacing
#    print(values)
#    print(directory_names)
#    exit()   
 
    # POLYMORPH
    polymorphs = polymorphs.split()

    # =============================================================================================
    # READ IN RAW DATA
    # =============================================================================================
    # Constants.
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)  # Boltzmann constant in kcal/mol
    
    # Parameters
    T_k = Temp * np.ones(len(directory_names), float)  # Convert temperatures to floats
    print(T_k)
  #  print(values)

    K = len(directory_names)  # How many states?
     
    # total number of states examined; 0 are unsampled if bonds are left on, 1 is unsampled if the bonds are removed
    Kbig = K
    
    # maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
    N_max = 5000
    
    # beta factor for the different temperatures
    beta_k = 1.0 / (kB * T_k)
    dA = np.zeros([len(polymorphs), Kbig], float)
    ddA = np.zeros([len(polymorphs), Kbig], float)
    convert_units = 0.2390057 * np.ones(Kbig, float)  # Convert all energies to kcal/mol
    
    # Allocate storage for simulation data
    for i, poly in enumerate(polymorphs):
        # N_k[k] is the total number of snapshots from alchemical state k
        N_k = np.zeros([Kbig], np.int32)
    
        # N_k_s[k,s] is the total number of snapshots from alchemical state k from seed s
        N_k_s = np.zeros([Kbig], np.int32)
    
        # u_kln[k,l,n] is the adjusted energy of snapshot n from simulation k
        u_kln = np.zeros([K, Kbig, N_max], np.float64)
    
        # dhdl_kn[k,n] is the derivative of energy with respect to lambda of snapshot n from simulation k
        dhdl_kn = np.zeros([K, N_max], np.float64)
    
        #Load in the data for each run
        for k in range(K):
            n = 0

            # cycle through all the input total energy data
            if directory_names[k] == int(directory_names[k]):
                dirpath = polymorphs[i] + '/' + primary_directory + '/' + str(int(directory_names[k]))
            else:
                dirpath = polymorphs[i] + '/' + primary_directory + '/' + str(directory_names[k])
            if os.path.isdir(dirpath):
                fname = dirpath + '/PROD.edr'
                dhdlname = dirpath + '/dhdl_PROD.xvg'

                potential_energy = panedr.edr_to_df(fname)['Potential'].values
                print("loading " + fname)

                dhdl_energy = np.loadtxt(dhdlname, comments=['#', '$', '@', '!'])
                print("loading " + dhdlname)

                # Removing any non-equilibrated points of the simulation
                [start_production, _, _] = timeseries.detectEquilibration(potential_energy)
                potential_energy = potential_energy[start_production:]
                dhdl_energy = dhdl_energy[start_production:,:]

                # Cutting points if they exceed N_max
                if len(potential_energy) > N_max:
                    potential_energy = potential_energy[len(potential_energy) - N_max:]
                    dhdl_energy = dhdl_energy[len(dhdl_energy) - N_max:,:]

                # the energy of every configuration from each state evaluated at its sampled state
                n = len(potential_energy)
                dhdl_placement = len(dhdl_energy[0, :]) - K
                u_kln[k, :K, :n] = (potential_energy.reshape((n, 1)) + dhdl_energy[:, dhdl_placement:]).T * convert_units[k]
                dhdl_kn[k, :n] = (float(Independent) / Molecules) * \
                                 np.sum(dhdl_energy[:, 2:dhdl_placement], axis=1) * convert_units[k]

                N_k_s[k] = n
                N_k[k] = n

        # convert to nondimensional units from kcal/mol
        u_kln *= beta_k[0]
    
        #u_kln_save = u_kln.copy()
        u_kln_save = u_kln[:]
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
            neff = 1 / np.sum(w ** 2)
            print("Effective number of sample in state %d is %10.3f" % (k, neff))
            print("Efficiency for state %d is %d/%d = %10.4f" % (k, neff, len(w), neff / len(w)))
    
        # extract self-consistent weights and uncertainties
        (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()
 
        print("Free Energies Optained...")
    
        # convert PMF to kcal/mol and normalize by the number of molecules
        df_i /= (beta_k[0] * float(Independent))
        ddf_i /= (beta_k[0] * float(Independent))
    
        dA[i, :] = df_i[-1]

        # =============================================================================================
        # COMPUTE UNCERTAINTY USING THE UNCORRELATED DATA
        # =============================================================================================
        
        for k in range(K):
            N_k[k] = 0
            n_old = 0

            g_k[k] = timeseries.statisticalInefficiency(dhdl_kn[k, n_old: (n_old + N_k_s[k])])
            print("Correlation time for sampled state %d is %10.3f" % (k, g_k[k]))
            # subsample the data to get statistically uncorrelated data
            indices = np.array(timeseries.subsampleCorrelatedData(u_kln[k, k, n_old:(n_old + N_k_s[k])],
                                                                          g=g_k[k]))  # subsample
    
            # not sure why we have to transpose
            if indices != []:
                u_kln[k, :, N_k[k]: (N_k[k] + len(indices))] = u_kln_save[k, :, (indices + n_old)].transpose()
                N_k[k] = N_k[k] + len(indices)
                n_old += N_k_s[k]

        print("Number of retained samples")
        print(N_k)
        print("Number of retained samples from each seed")
        print(N_k_s)
    
        # generate the weights of each of the umbrella set
        mbar = pymbar.MBAR(u_kln, N_k, verbose=True, subsampling_protocol=[{'method': 'L-BFGS-B'}])
    
        print("MBAR Converged...")
    
        # extract self-consistent weights and uncertainties
        try:
            (df_u, ddf_u, theta_i) = mbar.getFreeEnergyDifferences()
        except ValueError:
            pass
    
        print("Free Energies Optained...")
    
        # convert PMF to kcal/mol and normalize by the number of molecules
        df_u /= (beta_k[0] * float(Independent))
        ddf_u /= (beta_k[0] * float(Independent))
    
        ddA[i, :] = ddf_u[-1]
#        ddA[i, :] = ddf_i[-1]
        
        # Write out free energy differences
        print("Free Energy Difference (in units of kcal/mol)")
        print("  dA(Gamma) = A(Gamma) - A(Interactions Off)")
        for k in range(Kbig):
            print("%8.3f %8.3f" % (df_i[k, -1], ddf_u[k, -1]))

        del N_k
        del N_k_s
        del u_kln
        del dhdl_kn

    out_dA = np.zeros(len(polymorphs))
    out_ddA = np.zeros(len(polymorphs))
    for i, poly in enumerate(polymorphs):
        out_dA[i] = dA[i, 0]
        out_ddA[i] = ddA[i, 0]

    return out_dA, out_ddA

