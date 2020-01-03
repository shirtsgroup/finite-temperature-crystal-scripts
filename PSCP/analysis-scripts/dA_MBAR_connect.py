#!/bin/python
#
# Create a plot of free energy vs temperature for a polymorph
# 
# Copyright Michael R. Shirts, University of Virginia, 2014
#
from __future__ import print_function
import numpy as np
import pymbar  # multistate Bennett acceptance ratio
from pymbar import timeseries  # timeseries analysis
import panedr


def dA_MBAR_connect(end_point, start_point, Temperature=200, Molecules=72, Independent=0,
          phase='solid', Polymorphs=['p1', 'p2', 'p3']):

    if Independent == 0:
        Independent = Molecules

    # =============================================================================================
    # READ IN RAW DATA
    # =============================================================================================
    # Constants.
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)  # Boltzmann constant in kcal/mol/K
    
    # Parameters
    # How many states?
    K = 2
    
    #  maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
    N_max = 30000
    
    # beta factor for the different temperatures
    beta_k = 1.0 / (kB * Temperature)
    
    # Conversion from kJ to kcal
    kJ_to_kcal = 0.2390057

    # This is the sampling efficiency for each potential in each combination of potentials
    Efficiency = np.zeros(K, float)

    # N_k[k] is the total number of snapshots from alchemical state k
    N_k = np.zeros(K, np.int32)

    # dA[p,i,k] is the free energy between potential 0 and state k in polymorph p
    dA = np.zeros([len(Polymorphs), K], float)
    
    # ddA[p,i,k] is the uncertainty in the free energy between potential 0 and state k in polymorph p
    ddA = np.zeros([len(Polymorphs), K], float)
    
    # u_kn[k,n] is the reduced potential energy of configuration n from potential k
    u_kn = np.zeros([K, N_max], np.float64)

    # directory path
    dirpath = [start_point, end_point]

    # Cycle through all polymorphs
    for p, polymorph in enumerate(Polymorphs):
        count = 0
        for k in range(K):
            count += 1
            all_energy = panedr.edr_to_df(dirpath[k] + 'PROD.edr')
            [start_production, _, _] = timeseries.detectEquilibration(np.array(all_energy['Potential']))

            # Setting the end point of the simulation
            N = len(np.array(all_energy['Total Energy'])[start_production:])
            N_k[k] = N
    
            u_kn[k, :N] = np.array(all_energy['Potential'])[start_production:]

        print("Start1")
        # Convert all units to kcal
        u_kn *= kJ_to_kcal
        
        print("Start2")
        # Remove dependent molecules
        u_kn *= float(Independent) / Molecules
    
        print("Start3")
        # convert to nondimensional units from kcal/mol
        for k, beta in enumerate(beta_k):
            u_kn[k, :] *= beta
    
        u_kn_save = np.copy(u_kn)
    
        print("Number of retained samples")
        print(N_k)
    
        # =============================================================================================
        # COMPUTE FREE ENERGY DIFFERENCE USING MBAR FOR EACH SPACING
        # =============================================================================================
        for i in range(K):
            if i == 0:
                continue
            # Initialize MBAR
            print("Running MBAR...")

            # generate the weights of each of the umbrella set
            mbar = pymbar.MBAR(u_kn, N_k[i], verbose=True)
            print("MBAR Converged...")
       
            hold = mbar.computeEffectiveSampleNumber(verbose=True)
            print(hold)
             
            # extract self-consistent weights and uncertainties
            (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()
        
            # Store the dimensionless results in the dA container
            print(df_i)
            dA[p, i, :] = df_i[1]
        print(dA)
        
        # =============================================================================================
        # COMPUTE UNCERTAINTY USING MBAR
        # =============================================================================================
        g_k = np.zeros([K])
        for k in range(K):
            # subsample correlated data - for now, use energy from current state
            if N_k[k] > 0:
                print(N_k[k])
                g_k[k] = timeseries.statisticalInefficiency(u_kn_save[k, 0:100])
                print("Correlation time for phase (%s), sampled state %d is %10.3f" % (phase, k, g_k[k]))

                # subsample the data to get statistically uncorrelated data
                indices = np.array(timeseries.subsampleCorrelatedData(u_kn_save[k, 0:N_k[k]],
                                                                          g=g_k[k]))
                N_k[k] = len(indices)
                u_kn[k, :, 0:N_k[k]] = u_kn_save[k, indices].transpose()
    
        print("Number of retained samples")
        print(N_k)
    
        print("Running MBAR...")
    
        # generate the weights of each state
        mbar = pymbar.MBAR(u_kn, N_k[:], verbose=True)
        print("MBAR Converged...")
    
        # extract self-consistent weights and uncertainties
        (df_u, ddf_u, theta_u) = mbar.getFreeEnergyDifferences()
    
        # testing
        weights_in_gromos = np.zeros(K, float)
        for k in range(K):
            w = np.exp(mbar.Log_W_nk[:, k])
            print("max weight in state %d is %12.7f" % (k, np.max(w)))
            neff = 1 / np.sum(w ** 2)

            print("Effective number of sample in state %d is %10.3f" % (k, neff))
            print("Efficiency for state %d is %d/%d = %10.4f" % (k, neff, len(w), neff / len(w)))
            Efficiency[k] = neff / len(w)  # Store the efficiency
            w_0 = np.exp(mbar.Log_W_nk[:, 0])  # Weights in gromos
            initial_configs = np.sum(N_k[0:k])
            final_configs = np.sum(N_k[0:k + 1])

            print("Total weight in gromos " + str(np.sum(w_0[initial_configs:final_configs])))
            weights_in_gromos[k] = np.sum(w_0[initial_configs:final_configs])
    
        # Store the dimensionless results in the ddA container
        ddA[p, :] = ddf_u[1]
    return dA, ddA



