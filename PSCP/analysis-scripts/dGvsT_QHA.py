#!/bin/python
#
# Create a plot of free energy vs temperature for a polymorph
# 
# Copyright Michael R. Shirts, University of Virginia, 2014
#
from __future__ import print_function
import numpy as np
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
import mdtraj as md
import os
import panedr
import subprocess
import scipy as sp


def crystal_matrix_to_lattice_parameters(crystal_matrix):
    """
    This function takes any strained crystal lattice matrix and return the lattice parameters

    **Required Inputs
    crystal_matrix = crystal lattice matrix ([[Vxx,Vxy,Vxz],
                                              [Vyx,Vyy,Vyz],
                                              [Vzx,Vzy,Vzz]])
    """
    # Computing lattice parameters
    a = np.linalg.norm(crystal_matrix[:, 0])
    b = np.linalg.norm(crystal_matrix[:, 1])
    c = np.linalg.norm(crystal_matrix[:, 2])

    gamma = np.arccos(np.dot(np.squeeze(np.asarray(crystal_matrix[:, 0])), np.squeeze(np.asarray(crystal_matrix[:, 1])))
                      / (a * b)) * 180. / np.pi
    alpha = np.arccos(np.dot(np.squeeze(np.asarray(crystal_matrix[:, 1])), np.squeeze(np.asarray(crystal_matrix[:, 2])))
                      / (b * c)) * 180. / np.pi
    beta = np.arccos(np.dot(np.squeeze(np.asarray(crystal_matrix[:, 2])), np.squeeze(np.asarray(crystal_matrix[:, 0])))
                     / (c * a)) * 180. / np.pi

    # Creating an array of lattice parameters
    lattice_parameters = np.array([a, b, c, alpha, beta, gamma])
    return lattice_parameters

def compute_gradient(Tr, T, dA, T1, T2, G1, G2):
    dG = np.interp(Tr, T2, G2) - np.interp(Tr, T1, G1)
    dG_dTr = ((np.interp(Tr + 3, T2, G2) - np.interp(Tr + 3, T1, G1)) - (np.interp(Tr - 3, T2, G2) - np.interp(Tr - 3, T1, G1))) / 6
    df_dTr = (np.interp(Tr + 3, T, dA) - np.interp(Tr - 3, T, dA)) / 6
    return df_dTr, dG_dTr, dG

def compute_error_gradient(Tr, T, dA, T1, T2, G1, G2):
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)
    df_dTr, dG_dTr, dG = compute_gradient(Tr, T, dA, T1, T2, G1, G2)
    return np.absolute(-kB*df_dTr - 1 / Tr**2 * dG + 1 / Tr * dG_dTr)

def QHA_optimized_ref(T, dA, refT_files, refG_files):
    refT = np.ones(len(refT_files))
    refdG = np.zeros(len(refT_files))
    for i in range(1, len(refT)):
        T1 = np.load(refT_files[0])
        T2 = np.load(refT_files[i])
        G1 = np.load(refG_files[0])
        G2 = np.load(refG_files[i])
        x = sp.optimize.minimize(compute_error_gradient, 50.0, args=(T, dA[i], T1, T2, G1, G2),  tol=1e-8)
        print(x)
        refT[i] = x.x
        refdG[i] = np.interp(x.x, T2, G2) - np.interp(x.x, T1, G1)

    x = np.arange(10,300,1)
    y = np.zeros((len(refT) - 1, len(x)))
    for i in range(1, len(refT)):
        T1 = np.load(refT_files[0])
        T2 = np.load(refT_files[i])
        G1 = np.load(refG_files[0])
        G2 = np.load(refG_files[i])
        for j in range(len(x)):
            y[0, j] = compute_error_gradient(x[j], T,dA[i], T1, T2, G1, G2)
    np.save('x', x)
    np.save('y', y) 

    return refT, refdG


def compute_COV_dGref(refT_cov, Temperatures_MD, Molecules, Polymorphs):
    # setting a place to store the reference free energy differences
    refdG = np.zeros((len(refT_cov), len(Polymorphs)))

    # setting key variables for QHA
    natoms = len(md.load(Polymorphs[0] + '/temperature/0/pre_EQ.gro').xyz[0,:,0])
    nmodes = natoms * 3

    # boltzmann constant in kcal/(mol * K)
    kB = 0.0019872041

    # converting kcal to g*nm**2 / (ps**2)
    ekcal = 418.4

    # speed of light in cm / ps
    speed_of_light = 0.0299792458

    # Reduced planks constant
    h_bar = 2.520 * 10 ** (-38)

    # Avogadro's number
    Na = 6.022 * 10 ** 23

    for i, t in enumerate(refT_cov):
        # determining what directory to look into for this temperature
        directory = '/temperature/' + str(np.where(t == Temperatures_MD)[0][0])
        for j, p in enumerate(Polymorphs):
            path = p + directory
            edr = panedr.edr_to_df(path + '/PROD.edr')

            if not os.path.isfile(path + '/eigenvalues.xvg'):
                # Generating the eigenvalues from the covarience matrix
                c = subprocess.Popen(['echo', '0', ';', 'echo', '0'], stdout=subprocess.PIPE)
                output = subprocess.check_output(
                    ['gmx', 'covar', '-f', path + '/PROD.trr', '-s', path + '/PROD_0.tpr', '-o',
                     path + '/eigenvalues.xvg', '-mwa', 'yes', '-pbc', 'yes', '-last', str(nmodes)], stdin=c.stdout)
                c.wait()

                # Removing excess files that take up too much space
                subprocess.call(['rm', 'eigenvec.trr', 'covar.log', 'average.pdb'])

            # Loading in eigenvalues and converting them to wavenumbers
            wavenumbers = np.loadtxt(path + '/eigenvalues.xvg', comments=['#', '@'])[:, 1]
            #wavenumbers = kB * t / (np.absolute(wavenumbers[np.where(wavenumbers > 0.)])*100)
            wavenumbers = kB * t / (np.absolute(wavenumbers)*100)
            wavenumbers = np.sort(np.sqrt(wavenumbers[3:] * ekcal) / (2 * np.pi * speed_of_light))
            print(len(wavenumbers))
            # Getting the potential energy
            U = np.average(edr['Potential'].values) / 4.184

            # Computing the vbirational energy
            Av = kB * t * np.sum(np.log(Na * h_bar * wavenumbers * speed_of_light * 10 ** 12 / (kB * t)))
            print(U/Molecules, Av/Molecules, np.average(edr['Volume'].values) * Na * 0.024201 * 10 ** (-24) / Molecules)
            # Computing the free energy
            refdG[i, j] = (U + Av + np.average(edr['Volume'].values) * Na * 0.024201 * 10 ** (-24) ) / Molecules

    refdG -= refdG[:, 0]
    print(refdG, refT_cov)
    exit()
    return np.array(refT_cov), refdG

def dGvsT_QHA(Temperatures_MD=np.array([100,200,300]), Temperatures_unsampled=[], Molecules=72, molecule='benzene',
              Independent=0, potential='oplsaa', spacing=1, phase='solid', Polymorphs=['p1', 'p2', 'p3'],
              refdG_type='QHA',output_directory='output_QHA',
              refT_files=['', '', ''], refG_files=['', '', ''], refT_cov=[]):
    # Setting-up if the simulation is suppose to use QHA or covarience for dG ref
    if refdG_type == 'QHA':
        # Loading in the longest string of temperatures for refT
        refT = []
        for i in refT_files:
            temp_T = np.load(i)
            if len(temp_T) > len(refT):
                refT = np.load(i)

        # Cutting off any zero values form refT
        if refT[0] == 0.:
            refT = refT[1:]

        # Adding in the reference free energy differences for each temperature
        refdG = np.zeros((len(refT), len(Polymorphs)))
        for i in range(len(Polymorphs)):
            G0 = np.load(refG_files[0])
            G1 = np.load(refG_files[i])

            T0 = np.load(refT_files[0])
            T1 = np.load(refT_files[i])
            for j, t in enumerate(refT):
                placement_0 = np.where(T0 == t)
                placement_1 = np.where(T1 == t)
                if (len(placement_0[0]) == 1) and (len(placement_1[0]) == 1):
                    refdG[j, i] = G1[placement_1[0]] - G0[placement_0[0]]
                else:
                    refdG[j, i] = np.nan

    elif refdG_type == 'COV':
        if output_directory == 'output_QHA':
            output_directory = 'output_COV'
        refT, refdG = compute_COV_dGref(refT_cov, Temperatures_MD, Molecules, Polymorphs)

    else:
        print("ERROR: refdG_type " + refdG_type + " is not a valid input.")
        exit()

    if Independent == 0:
        Independent = Molecules
    # =============================================================================================
    # Load reference free energy differences
    # =============================================================================================


    # Hard set from old dictionary funciton
    refPot = 0
    ExtraPressures = []
    Temperatures = np.sort(np.append(Temperatures_MD, Temperatures_unsampled))
    Temperatures = np.sort(np.unique(np.append(Temperatures, refT)))
    Pressures = np.ones(len(Temperatures), int)
    Pressures[len(Pressures) - len(ExtraPressures): len(Pressures)] = ExtraPressures
    Potentials = [potential]
    
    # =============================================================================================
    # FORMAT INPUTS
    # =============================================================================================
    # TEMPERATURE
    refk = []
    for k, temp in enumerate(refT):
        refk.append(np.where(temp == Temperatures)[0][0])

    # =============================================================================================
    # READ IN RAW DATA
    # =============================================================================================
    # Constants.
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)  # Boltzmann constant in kcal/mol/K
    
    # Parameters
    # How many states?
    K = len(Potentials) * len(Temperatures)
    
    #  maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
    N_max = 30000
    
    # beta factor for the different temperatures
    beta_k = 1.0 / (kB * Temperatures)
    beta_k = np.tile(beta_k, (1, len(Potentials)))[0]
    
    # Conversion from kJ to kcal
    kJ_to_kcal = 0.2390057

    # This is the sampling efficiency for each potential in each combination of potentials
    Efficiency = np.zeros(K, float)
    
    # N_k[k] is the total number of snapshots from alchemical state k
    N_k = np.zeros(K, np.int32)
    
    # dA[p,i,k] is the p.interp(Tr, T1, G1)free energy between potential 0 and state k for spacing i in polymorph p
    dA = np.zeros([len(refT), len(Polymorphs), spacing + 1, K], float)
    
    # ddA[p,i,k] is the uncertainty in the free energy between potential 0 and state k for spacing i in polymorph p
    ddA = np.zeros([len(refT), len(Polymorphs), spacing + 1, K], float)

    run_dA_analysis = True

    if os.path.isdir(output_directory + '/dA_raw.npy') and os.path.isdir(output_directory + '/ddA_raw.npy'):
        hold_dA = np.load(output_directory + '/dA_raw.npy')
        if np.shape(hold_dA) == np.shape(dA):
            dA = np.load(output_directory + '/dA_raw.npy')
            ddA = np.load(output_directory + '/ddA_raw.npy')
            run_dA_analysis = False

    # dG[p,i,t] is the free energy between polymorph 1 and polymorph p for spacing i and temperature t
    dG = np.zeros([len(refT), len(Polymorphs), spacing + 1, len(Temperatures)])
    
    # ddG[p,i,t] is the uncertanity in the free energy between polymorph 1 and polymorph p for spacing i and temperature t
    ddG = np.zeros([len(Polymorphs), spacing + 1, len(Temperatures)])
    
    # dS[p,i,t] is the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
    dS = np.zeros([len(refT), len(Polymorphs), spacing + 1, len(Temperatures)])
    
    # ddS[p,i,t] is the uncertanity in the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
    ddS = np.zeros([len(Polymorphs), spacing + 1, len(Temperatures)])
    
    dS_mbar = np.zeros([len(Polymorphs), spacing + 1, len(Temperatures)])
    ddS_mbar = np.zeros([len(Polymorphs), spacing + 1, len(Temperatures)])
    dH_mbar = np.zeros([len(Polymorphs), spacing + 1, len(Temperatures)])
    
    # O_pij[p,i,j] is the overlap within polymorph p between temperature state i and temperature state j
    O_pij = np.zeros([len(Polymorphs), len(Temperatures), len(Temperatures)])
    dU = np.zeros([len(Polymorphs), len(Temperatures)])
    ddU = np.zeros([len(Polymorphs), len(Temperatures)])
    
    # u_kln[k,l,n] is the reduced potential energy of configuration n from potential k in potential l
    u_kln = np.zeros([K, K, N_max], np.float64)
    
    # V_pkn is the volume of configuration n of polymorph p at temperature k
    V_pkn = np.zeros([len(Polymorphs), len(Temperatures), N_max], float)

    # V_avg is the average volume of polymorph p at temperature k
    V_avg = np.zeros([len(Polymorphs), len(Temperatures)], float)
    
    # ddV_avg is the standard deviation of the volume of polymorph p at temperature k
    ddV_avg = np.zeros([len(Polymorphs), len(Temperatures)], float)

    # C_pkn is the lattice tensor of the polymorph p at temperature k
    box_place = np.matrix([[0, 0], [1, 1], [2, 2], [0, 1], [0, 2], [1, 2]])
    C_pkn = np.zeros([len(Polymorphs), len(Temperatures), N_max, 3, 3], float)

    # h_avg is the average lattice parameters of polymorph p at temperature k
    h_avg = np.zeros([len(Polymorphs), len(Temperatures), 6], float)
    
    # dh is the standard deviation of the lattice parameters of polymorph p at temperature k
    dh = np.zeros([len(Polymorphs), len(Temperatures), 6], float)

    # Cycle through all polymorphs
    if run_dA_analysis:
        for p, polymorph in enumerate(Polymorphs):
            # Cycle through all sampled potentials
            for i, potential_k in enumerate(Potentials):
                count = 0
                for t in range(len(Temperatures)):
                    k = len(Temperatures) * i + t
                    # Cycle through all evaluated potentials
                    for j, potential_l in enumerate(Potentials):
                        l = len(Temperatures) * j
        
                        dirpath = polymorph + '/temperature/' + str(count) + '/'
                        if os.path.isfile(dirpath + 'PROD.edr') and (Temperatures[t] in Temperatures_MD):
                            count += 1
                            print("loading " + dirpath + 'PROD.edr')
                            all_energy = panedr.edr_to_df(dirpath + 'PROD.edr')
                            if len(all_energy['Potential'].values) > N_max:
                                [start_production, _, _] = timeseries.detectEquilibration(all_energy['Potential'].values[::10])
                                start_production *= 10
                            else:
                                [start_production, _, _] = timeseries.detectEquilibration(all_energy['Potential'].values)

                            # Now read in the lattice tensor and average them
                            if 'Box-XX' in list(all_energy):
                                box_letters = ['XX', 'YY', 'ZZ', 'YX', 'ZX', 'ZY']
                            else:
                                box_letters = ['X', 'Y', 'Z']
        
                            for b in range(len(box_letters)):
                                if len(all_energy['Potential'].values) > N_max:
                                    [hold,_,_] = timeseries.detectEquilibration(all_energy['Box-' + box_letters[b]].values[::10])
                                    hold *= 10
                                else:
                                    [hold,_,_] = timeseries.detectEquilibration(all_energy['Box-' + box_letters[b]].values)

                                if hold > start_production:
                                    start_production = hold

                            if len(all_energy['Total Energy'].values[start_production:]) > N_max:
                                start_production = len(all_energy['Total Energy'].values) - N_max

                            # Setting the end point of the simulation
                            N = len(all_energy['Total Energy'].values[start_production:])
                            N_k[k] = N
        
                            u_kln[k, l, :N] = all_energy['Potential'].values[start_production:]
        
                            # Now set these energies over all temperatures
                            u_kln[k, l:(l + len(Temperatures)), :N] = u_kln[k, l, :N]
        
                            # Now read in the volumes and average them
                            V_pkn[p, t, :N] = all_energy['Volume'].values[start_production:]
                            V_avg[p, t] = np.average(V_pkn[p, t, :N]) / float(Independent)
                            ddV_avg[p, t] = np.std(V_pkn[p, t, :N]) / float(Independent)
        
                            # Making the lattice tensor all the correct sign with time    
                            if count == 1:
                                sign = np.sign(md.load(dirpath + 'pre_EQ.gro').unitcell_vectors[0].T)
                                for s in range(3):
                                    for j in range(3):
                                        if sign[s, j] == 0.:
                                            # Correcting for the sign of the lattice parameters
                                            sign[s, j] = 1.

                            for b in range(len(box_letters)):
                                C_pkn[p, t, :N, box_place[b, 0], box_place[b, 1]] = np.absolute(all_energy['Box-' + box_letters[b]].values[start_production:]) * \
                                        sign[box_place[b, 0], box_place[b, 1]] * 10
                            C_avg = np.average(C_pkn[p, t, :N], axis=0)
                            dC = np.std(C_pkn[p, t, :N], axis=0)
                            h_avg[p, t] = crystal_matrix_to_lattice_parameters(C_avg) 
                            dh[p, t] = np.absolute(crystal_matrix_to_lattice_parameters(C_avg + dC) - h_avg[p, t])
                        else:
                            N_k[k] = 0
                            V_avg[p, t] = np.nan
                            ddV_avg[p, t] = np.nan
                            h_avg[p, t] = np.nan
                            dh[p, t] = np.nan

            print("Start1")
            # Convert all units to kcal
            #u_pklnT[p, :, :, :] *= kJ_to_kcal
            u_kln *= kJ_to_kcal
            
            print("Start2")
            # If this was already in kcal or already fully independent, revert
            for j in range(len(Potentials)):
                if Potentials[j][:6] == "amoeba":
                    #u_pklnT[p, :, j * len(Temperatures):(j + 1) * len(Temperatures), :, :] /= kJ_to_kcal
                    u_kln[:, j * len(Temperatures):(j + 1) * len(Temperatures), :] /= kJ_to_kcal
            
            print("Start3")
            # Remove dependent molecules
            for j in range(len(Potentials)):
                if Potentials[j][:6] != "amoeba":
                    #u_pklnT[p, :, j * len(Temperatures):(j + 1) * len(Temperatures), :, :] *= float(Independent) / Molecules
                    u_kln[:, j * len(Temperatures):(j + 1) * len(Temperatures), :] *= float(Independent) / Molecules
        
            print("Start4")
            # Now average together the energies and volumes at each state
            for t in range(len(Temperatures)):
                dU[p, t] = np.average(u_kln[t, t, :N_k[t]]) / float(Independent)
                ddU[p, t] = np.std(u_kln[t, t, :N_k[t]]) / N_k[t] ** 0.5 / float(Independent)
        
            print("Start5")
            # convert to nondimensional units from kcal/mol
            for k, beta in enumerate(beta_k):
                u_kln[:, k, :] *= beta
        
            u_kln_save = u_kln.copy()
            N_k_save = N_k.copy()
            print("End!")
        
            print("Number of retained samples")
            print(N_k)
        
            # Now create the full N_k matrix including the roll-backs as well as the free energy container
            # N_k_matrix[i,k] is the total number of snapshots from alchemical state k using in spacing i
            N_k_matrix = np.zeros([spacing + 1, K], np.int32)
            for i in range(spacing + 1):
                N_k_matrix[i, :] = N_k_save.copy()
                N_k_matrix[i, 0: len(Temperatures)] = N_k_matrix[i, 0:len(Temperatures)] * float(i) / float(spacing)
        
            # =============================================================================================
            # COMPUTE FREE ENERGY DIFFERENCE USING MBAR FOR EACH SPACING
            # =============================================================================================
            for i in range(spacing+1):
                if i == 0 and len(Potentials) == 1:
                    continue
                # Initialize MBAR.
                print("Running MBAR...")

                # generate the weights of each of the umbrella set
                mbar = pymbar.MBAR(u_kln, N_k_matrix[i, :], verbose=True)
                print("MBAR Converged...")
           
                hold = mbar.computeEffectiveSampleNumber(verbose=True)
                print(hold)
                 
                # extract self-consistent weights and uncertainties
                (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()

                # extract entropy
                [_, _, Delta_u_ij, _, Delta_s_ij, dDelta_s_ij] = mbar.computeEntropyAndEnthalpy()
                print("Free Energies Optained...")
            
                # Store the dimensionless results in the dA container
                dA[:, p, i, :] = df_i[refk]
                dH_mbar[p, i, :] = Delta_u_ij[0]
                dS_mbar[p, i, :] = Delta_s_ij[0]
                ddS_mbar[p, i, :] = dDelta_s_ij[0]
                print(dA)
            
            # =============================================================================================
            # COMPUTE UNCERTAINTY USING MBAR
            # =============================================================================================
            g_k = np.zeros([K])
            for i in range(spacing + 1):
                if i == 0 and len(Potentials) == 1:
                    continue

                for k in range(K):
                    # subsample correlated data - for now, use energy from current state
                    if N_k_matrix[i, k] > 0:
                        print(N_k_matrix[i, k])
                        g_k[k] = timeseries.statisticalInefficiency(u_kln_save[k, k, 0:100])
                        print("Correlation time for phase (%s), sampled state %d is %10.3f" % (phase, k, g_k[k]))

                        # subsample the data to get statistically uncorrelated data
                        indices = np.array(timeseries.subsampleCorrelatedData(u_kln_save[k, k, 0:N_k_matrix[i, k]],
                                                                              g=g_k[k]))
                        N_k_matrix[i, k] = len(indices)
                        u_kln[k, :, 0:N_k_matrix[i, k]] = u_kln_save[k, :, indices].transpose()  # not sure why we have to transpose
        
                print("Number of retained samples")
                print(N_k)
        
                print("Running MBAR...")
        
                # generate the weights of each state
                mbar = pymbar.MBAR(u_kln, N_k_matrix[i, :], verbose=True)
                print("MBAR Converged...") 
        
                # extract self-consistent weights and uncertainties
                (df_u, ddf_u, theta_u) = mbar.getFreeEnergyDifferences()
        
                # calculate the overlap it necessary
                if len(Temperatures) == 2:
                    O_pij[p, :, :] = mbar.computeOverlap()[2]
        
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
            
                # Write out free energy differences
                print("Free Energy Difference (in units of kcal/mol)")
                for k in range(K):
                    print("%8.3f %8.3f" % (-df_i[k, 0], ddf_u[k, 0]))
        
                # Store the dimensionless results in the ddA container
                ddA[:, p, i, :] = ddf_u[refk]
#            ddA[:, p, i, :] = ddf_i[refk]

        # Saving the files if needed for QHA
        if refdG_type == 'QHA':
            np.save(output_directory + '/dA_raw.npy', dA)
            np.save(output_directory + '/ddA_raw.npy', ddA)

    # =============================================================================================
    # FINALIZE THE RELATIVE FREE ENERGY AND ENTROPY
    # =============================================================================================
    for k in range(len(refT)):
        for i in range(spacing + 1):
            for t, T in enumerate(Temperatures):
                for p in range(len(Polymorphs)):
                    dG[k, p, i, t] = (dA[k, p, i, t] - dA[k, 0, i, t]) / (beta_k[t] * float(Independent)) + float(T) / float(
                        refT[k]) * refdG[k, p]
                    if p == 0:
                        continue
                    dS[k, p, i, t] = (dU[p, t] - dU[0, t] - dG[k, p, i, t]) / float(T)
                    if k == 0:
                        ddG[p, i, t] = ((ddA[k, p, i, t] ** 2 + ddA[k, 0, i, t] ** 2) / (
                                    beta_k[t] * float(Independent)) ** 2) ** 0.5
                        ddS[p, i, t] = (ddU[p, t] ** 2 + ddU[p, t] ** 2 + ddG[p, i, t] ** 2) ** 0.5 / float(T)

    
    # =============================================================================================
    # PLOT THE RELATIVE FREE ENERGY VS TEMPERATURE
    # =============================================================================================

    PlotPress = 1  # Pressure to plot the dGvT curve at
    Temperatures_P = Temperatures[Pressures == PlotPress]

    if not os.path.isdir(output_directory):
        subprocess.call(['mkdir', output_directory])

    np.save(output_directory + '/T_' + molecule + '_' + potential, Temperatures_P)
    for p, Poly in enumerate(Polymorphs):
        np.save(output_directory + '/dGvT_' + molecule + '_' + Poly + '_' + potential, dG[:, p, spacing, Pressures == PlotPress])
        np.save(output_directory + '/ddGvT_' + molecule + '_' + Poly + '_' + potential, ddG[p, spacing, Pressures == PlotPress])
        if len(Potentials) > 1:
            np.save(output_directory + '/dGvT_' + molecule + '_' + Poly + '_' + potential + '_indirect', dG[:, p, 0, :])
            np.save(output_directory + '/ddGvT_' + molecule + '_' + Poly + '_' + potential + '_indirect', ddG[p, 0, :])
            if spacing > 1:
                np.save(output_directory + '/dGvT_' + molecule + '_' + Poly + '_' + potential + '_convergence', dG[:, p, :, :])
                np.save(output_directory + '/ddGvT_' + molecule + '_' + Poly + '_' + potential + '_convergence', ddG[p, :, :])
        np.save(output_directory + '/dS_' + molecule + '_' + Poly + '_' + potential, dS[:, p, spacing, :])
        np.save(output_directory + '/ddS_' + molecule + '_' + Poly + '_' + potential, ddS[p, spacing, :])

    for p, Poly in enumerate(Polymorphs):
        np.save(output_directory + '/UvT_' + molecule + '_' + Poly + '_' + potential, dU[p, :])

    for p, Poly in enumerate(Polymorphs):
        np.save(output_directory + '/VvT_' + molecule + '_' + Poly + '_' + potential, V_avg[p, :])
        np.save(output_directory + '/dVvT_' + molecule + '_' + Poly + '_' + potential, ddV_avg[p, :])

    # =============================================================================================
    # SAVE THE AVERAGE BOX VECTORS AND ANGLES VS TEMPERATURE
    # =============================================================================================

    for p, Poly in enumerate(Polymorphs):
        np.save(output_directory + '/hvT_' + molecule + '_' + Poly + '_' + potential, h_avg[p, :])
        np.save(output_directory + '/dhvT_' + molecule + '_' + Poly + '_' + potential, dh[p, :])
    
    # Save the data for future use.
    for p, Poly in enumerate(Polymorphs):
        np.save(output_directory + '/dUvT_' + molecule + '_' + Poly + '_' + potential, dU[p, :] - dU[0, :])
        np.save(output_directory + '/ddUvT_' + molecule + '_' + Poly + '_' + potential, (ddU[p, :] ** 2 + ddU[0, :] ** 2) ** 0.5)


