import os
import numpy as np
import pymbar
import panedr

def dA_endpoint_MBAR(polymorphs='p1 p2', Molecules=72, Independent=4, Temp=200):
    # Setting constants
    kJ_to_kcal = 1/4.184  # Converting kJ to kcal
    kB = 0.0019872041  # boltzman constant in kcal/(mol*K)

    # Getting the polymorph names
    polymorphs = polymorphs.split()

    # Place to store the free energy differences
    dA = np.zeros(len(polymorphs))
    ddA = np.zeros(len(polymorphs))

    for i, poly in enumerate(polymorphs):
        if os.path.isfile(poly + '/interactions/100/PROD.edr') and os.path.isfile(poly + '/interactions/100/END.edr'):
            # Loading in the difference between the endpoint and production files
            dU = panedr.edr_to_df(poly + '/interactions/100/END.edr')['Potential'].values - panedr.edr_to_df(poly + '/interactions/100/PROD.edr')['Potential'].values

            # Converting the energy differences to go into pymbar
            dW = dU / Molecules * kJ_to_kcal / (kB * Temp)

            # Getting the energy differences with MBAR using Exponential Averaging
            da = np.array(pymbar.EXP(dW)) * kB * Temp

            dA[i] = -da[0]
            ddA[i] = da[1]
        else:
            dA[i] = np.nan
            ddA[i] = np.nan

    # Check to see if there are any nan values
    if np.any(dA == np.nan):
        dA[:] = 0.
        ddA[:] = 0.

    return dA, ddA




