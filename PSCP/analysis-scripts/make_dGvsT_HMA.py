# !/bin/python
#
# Create a plot of free energy vs temperature for a polymorph
#
import numpy as np
from pymbar import timeseries  # timeseries analysis
import mdtraj as md
import os
import panedr
from scipy import integrate

def anharmonic_energy(temperature, average_U_anhar_points, T_points):
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)
    return - 1 / (kB * temperature) * np.interp(temperature, T_points, average_U_anhar_points)

def make_dGvsT_HMA(Temperatures_MD=np.array([100, 200, 300]), Molecules=72, molecule='benzene', potential='oplsaa',
                   Polymorphs=['p1', 'p2', 'p3'], output_directory='output_HMA', QHA_raw_files=['', '', ''],
                   HMA_T=[50.]):

    # =============================================================================================
    # READ IN RAW DATA
    # =============================================================================================
    Temperatures_MD = np.sort(np.unique(np.array(Temperatures_MD))).astype(float)

    # Constants.
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)  # Boltzmann constant in kcal/mol/K

    # Parameters
    # How many states?
    K = len(HMA_T)

    #  maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
    N_max = 30000

    # Conversion from kJ to kcal
    kJ_to_kcal = 1 / 4.184

    # N_k[k] is the total number of snapshots from alchemical state k
    N_k = np.zeros(K, np.int32)

    # dA[p,i,k] is the p.interp(Tr, T1, G1)free energy between potential 0 and state k for spacing i in polymorph p
    dA_anharmonic = np.zeros([len(HMA_T), len(Polymorphs)], float)

    # ddA[p,i,k] is the uncertainty in the free energy between potential 0 and state k for spacing i in polymorph p
    #ddA = np.zeros([len(HMA_T), len(Polymorphs)], float)

    # Cycle through all polymorphs
    natoms = len(md.load(Polymorphs[0] + '/temperatures/0/pre_EQ.gro').xyz[0,:,0])
    for p, polymorph in enumerate(Polymorphs):
        # Cycle through all sampled potentials
        raw_QHA = np.laod(QHA_raw_files[p])
        dA_anhar_unint = np.zeros(len(Temperatures_MD))
        for t in range(len(Temperatures_MD)):
            U_MD = np.zeros(N_max)
            F_dot_dr = np.zeros(N_max)
            U_lat = np.interp(Temperatures_MD[t], raw_QHA[:,0], raw_QHA[:,3])

            placement = np.where(t == Temperatures_MD)[0][0]
            dirpath = polymorph + '/temperature/' + str(placement) + '/'
            if os.path.isfile(dirpath + 'FORCES.edr'):
                print("loading " + dirpath + 'FORCES.edr')
                # Loading in the energy file
                all_energy = panedr.edr_to_df(dirpath + 'FORCES.edr')
                [start_production, _, _] = timeseries.detectEquilibration(all_energy['Potential'].values)
                N = len(all_energy['Potential'].values[start_production:])
                U_MD[:N] = all_energy['Potential'].values[start_production:]

                # Loading in the trajectory file
                trajectory = md.load(dirpath + 'FORCES.trr', top=dirpath + 'pre_EQ.gro')[start_production:]
                r = np.average(trajectory.xyz, axis=0)
                F = np.loadtxt(dirpath + 'forces.xvg', comments=['#','@'])[start_production:,1:]
                for k in range(F[:,0]):
                    dr = (trajectory.xyz[k] - r).flatten()
                    F_dot_dr[k] = np.dot(F[k], dr[k])

                dA_anhar_unint[t] = np.average(U_MD / Molecules * kJ_to_kcal - U_lat +
                                               0.5 * F_dot_dr/ Molecules * kJ_to_kcal)
            else:
                dA_anhar_unint[t] = np.nan

        T = Temperatures_MD[~np.isnan(dA_anhar_unint)]
        dA_anhar_unint = dA_anhar_unint[~np.isnan(dA_anhar_unint)]

        for i, t in enumerate(HMA_T):
            dA_anharmonic[i, p] = integrate.quad(anharmonic_energy, 0, t, args=(dA_anhar_unint, T))

    for p, ply in enumerate(Polymorphs):
        G_new = np.zeros(len(HMA_T))
        G_QHA = np.load(QHA_raw_files[p])[:,2]
        T_QHA = np.load(QHA_raw_files[p])[:,0]
        for i in len(HMA_T):
            G_new = np.interp(HMA_T[i], T_QHA, G_QHA) + dA_anharmonic[i, p]

        np.save(output_directory + '/dGvT_QHAnhar_' + molecule + '_' + ply + '_' + potential, G_new)
        np.save(output_directory + '/T_QHAnhar_' + molecule + '_' + ply + '_' + potential, HMA_T)




