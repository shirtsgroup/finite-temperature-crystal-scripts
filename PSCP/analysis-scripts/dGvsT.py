#!/bin/python
#
# Create a plot of free energy vs temperature for a polymorph
# 
# Copyright Michael R. Shirts, University of Virginia, 2014
#
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis
from optparse import OptionParser # for parsing command-line options
import MBARBootstrap # Bootstrapping algorithm
import random
import os
import usefulFuncs #Useful math functions
import Harvist #Hamiltonian Reweighting Visualization Toolkik
import pdb

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

#=============================================================================================
# READ IN USER INPUTS
#=============================================================================================
parser = OptionParser()
parser.add_option('-p', '--plot', dest = 'plot', help = 'Plot output (default false)', default=True, action = 'store_true')
parser.add_option('-T', dest = 'temp', help = 'Temperature', default = 200)
parser.add_option('-P', dest = 'Pressure', help = 'Pressure', default = 1)
parser.add_option('-n', dest = 'polymorphs', help = 'Polymorphs to analyze', default = 'p1 p2 p3')
parser.add_option('-N', dest = 'molecules', help = 'number of supercell molecules', default = 72)
parser.add_option('-M', dest = 'molecule', help = 'name of the molecule', default = 'benzene')
parser.add_option('-I', dest = 'independent', help = 'number of independent molecules', default = 'Same')
parser.add_option('-u', dest = 'potential', help = 'potential to create the phase diagram in', default = 'oplsaa')
parser.add_option('-i', dest = 'ignoreframes', help = 'Initial frames to ignore', default = 200)
parser.add_option('-j', dest = 'includeframes', help = 'Number of frames to include', default = 100000)
parser.add_option('-z', dest = 'simulation', help = 'The simulation package that was used', default = 'gromacs')
parser.add_option('-d', dest = 'directory', help = 'Parent directory of the reweight directories', default = 'None')
parser.add_option('-E', dest = 'ensemble', help = 'Simulation Ensemble', default = 'NVT')
parser.add_option('-s', dest = 'spacing', help = 'spacing for the plot rolling back the state 0 sampling', default = '1')
parser.add_option('-H', '--hinge', dest = 'hinge', help = 'Optional string at end of jobs', default = 'DefaultHinge')

(options, args) = parser.parse_args()
Temp = options.temp
Pressure = int(options.Pressure)
Molecules = int(options.molecules)
molecule = options.molecule
if options.independent == 'Same':
    Independent = Molecules
else:
    Independent = int(options.independent)
potential=options.potential
ignoreframes = int(options.ignoreframes)
includeframes = int(options.includeframes)
hinge = options.hinge
phase = "solid"
ensemble = options.ensemble
directory = options.directory
if directory == 'None':
    directory = molecule
spacing = int(options.spacing)


"""
Potentials=['gromos', 'designedg']
Charges=['0.1150', '0.1150']
Chargenames=['', 'C01150']
"""

Polys = dict() #Name of the polymorphs
refTs = dict() #Reference temperatures for the PSCP for each system
refdGs = dict() #Reference free energies for the PSCP for each system
refddGs = dict() #Reference uncertainties for the PSCP for each system
refdUs = dict() #Reference lattice minima for each system
absolutedUs = dict() #Absolute lattice energy for each system


#Oplsaa
if potential == "oplsaa":
    Potentials=['oplsaa']
    PotNAME='OPLS'
    Charges=['0.1150', '0.1150', '0.0700', '0.0850', '0.1000', '0.1150', '0.1300', '0.1450', '0.1600', '0.1750', '0.1900', '0.1000', '0.0900', '0.0800', '0.0700','0.0600','0.0500','0.0400','0.0300','0.0200','0.0100']
    Chargenames=['','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','']
    #Chargenames=['C01150', 'C01150', 'C00700', 'C00850', 'C01000', 'C01150', 'C01300', 'C01450', 'C01600', 'C01750', 'C01900', 'C00900', 'C00800', 'C00700', 'C00600', 'C00500', 'C00400', 'C00300', 'C00200', 'C00100']
    PotNAMES=['OPLS']
    SimNAMES=['GRO']

    #Temperatures=numpy.array([45,50]) #Overlap Check
    #Temperatures=numpy.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,70,80,90,100,110,120,130,140])#,80,90,100,110])#Zzzvye
    #Temperatures=numpy.array([10,15,20,25,30,35,40,45,50,60,70,80,90,100,110])#Zzzvye
    Temperatures=numpy.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340])#,360,380,400,420,440,460,480,500])#,320,340,360,380,400]) #benzene
    #Temperatures=numpy.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340])#,360,380,400,420,440,460,480,500])#,320,340,360,380,400]) #kobfud
    #Temperatures=numpy.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340,360,380,400,420,440,460,480,500])#,320,340,360,380,400]) #cbmzpn
    #Temperatures=numpy.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340])#,360,380,400,420,440,460,480,500])#,320,340,360,380,400]) #melfit
    #Temperatures=numpy.array([10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340,360,380,400,420,440,460,480,500])
    #Temperatures=numpy.array([100,130,140,160,180,200,220,240,260,280,300,340,360,380,400])
    #Temperatures=numpy.array([100,120,140,160,180,200]) #MelfitTest
    #ExtraPressures=numpy.array([5000,15000,35000]) #benzene
    #Temperatures=numpy.array([40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300]) #imazol
    #ExtraPressures=numpy.array([5000,10000,16000,20000,26000,31000,36000,40000]) #imazol
    #Temperatures=numpy.array([20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,220,240,260,280,300,300,300,300,300,300,300,300]) #acetac
    #ExtraPressures=numpy.array([5000,5000,5000,5000,5000,10000,16000,20000,26000,31000,36000,40000]) #acetac
    #Temperatures=numpy.array([20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,200,200,200,200,220,240,260,280,300,300,300,300,300,300,300,300]) #formam
    #ExtraPressures=numpy.array([1500,2000,5000,10000,10000,10000,10000,10000,5000,10000,15000,20000,26000,31000,36000,40000]) #formam
    #Temperatures=numpy.array([20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,250,270,290,300,320,340,360,380,400,420,440,460,480,500]) #glycin
    ExtraPressures=[]
    Pressures=numpy.ones(len(Temperatures),int);
    Pressures[len(Pressures)-len(ExtraPressures):len(Pressures)]=ExtraPressures
    #Temperatures=numpy.array([200,200,200,200,200])
    #Pressures=numpy.array([1,500,1000,1500,2000,2500,3000,3500,4000,4500,5000])
    refPot=0
    
    #Benzene
    Polys['benzene']=['p1', 'p2', 'p3']
    refTs['benzene']=200
    refdGs['benzene']=[0.000, 0.185, 0.306]
    refddGs['benzene']=[0.000, 0.019, 0.019]
    refdUs['benzene'] = [0.000, 0.267, 0.240]
    absolutedUs['benzene'] = [-5.624, -5.362, -5.380]
    
    #Formam
    Polys['formam']=['p1', 'p2']
    refTs['formam']=200
    refdGs['formam']=[0.000, 0.118]
    refddGs['formam']=[0.005, 0.005]
    refdUs['formam'] = [0.000, 0.343]
    absolutedUs['formam'] = [-41.04733, -40.7076]   
 
    #Acetac
    Polys['acetac']=['p1', 'p2']
    refTs['acetac']=200
    refdGs['acetac']=[0.000, 0.001]
    refddGs['acetac']=[0.005, 0.005]
    refdUs['acetac'] = [0.000, 0.020]
    absolutedUs['acetac'] = [-36.4858, -36.466]

    
    #Imazol
    Polys['imazol']=['p1', 'p2']
    refTs['imazol']=100
    refdGs['imazol']=[0.000, -0.436]
    refddGs['imazol']=[0.005, 0.005]
    refdUs['imazol'] = [0.000, -0.464]
    absolutedUs['imazol'] = [-15.399, -15.863]

    #Glycin
    Polys['glycin']=['p1', 'p2', 'p4']#, 'p3']
    refTs['glycin']=200
    refdGs['glycin']=[0.000, -0.034, -0.067]#, 1.087]
    refddGs['glycin']=[0.005, 0.005, 0.005]#, 0.005]
    refdUs['glycin'] = [0.000, 0.117, 0.288]#, 1.235]
    absolutedUs['glycin'] = [-131.950, -131.778, -131.609]#, -130.715]

    #Bismev
    Polys['bismev']=['p3', 'p1', 'p2']
    refTs['bismev']=200
    refdGs['bismev']=[0.000, -1.474, -0.186]
    refddGs['bismev']=[0.005, 0.005, 0.005]#, 0.005]
    refdUs['bismev'] = [0.000, -1.459, -0.245]#, 1.235]
    absolutedUs['bismev'] = [-64.3467, -65.8055, -64.5916]#, -130.715]

    #Cbmzpn
    Polys['cbmzpn']=['p3', 'p1']
    refTs['cbmzpn']=200
    refdGs['cbmzpn']=[0.000, 0.195]
    refddGs['cbmzpn']=[0.0, 0.005]
    refdUs['cbmzpn'] = [0.000, 0.565]
    absolutedUs['cbmzpn'] = [-34.3933, -33.8283]

    #Hxacan
    Polys['hxacan']=['p2', 'p1']
    refTs['hxacan']=200
    refdGs['hxacan']=[0.000, -0.619]
    refddGs['hxacan']=[0.005, 0.005]
    refdUs['hxacan'] = [0.000, -0.446]
    absolutedUs['hxacan'] = [-46.965, -47.411]

    #Kobfud
    Polys['kobfud']=['p1', 'p2']
    refTs['kobfud']=200
    refdGs['kobfud']=[0.000, 1.832]
    refddGs['kobfud']=[0.005, 0.005]
    refdUs['kobfud'] = [0.000, 2.255]
    absolutedUs['kobfud'] = [-49.974, -47.719]

    ##Pyrzin
    #Polys['pyrzin']=['p2']#, 'p3', 'p1', 'p5', 'p2']
    #refTs['pyrzin']=200
    #refdGs['pyrzin']=[0.000] #, -0.483, -0.689, -1.055, -1.056]
    #refddGs['pyrzin']=[0.005] #, 0.005, 0.005, 0.005, 0.005]
    #refdUs['pyrzin'] = [0.000] #, -0.176, -0.695, -1.246, 0.951]
    #absolutedUs['pyrzin'] = [-13.348] #, -14.475, -14.993, -15.545, -13.348]


    #Pyrzin
    Polys['pyrzin']=['p4', 'p3', 'p1', 'p5', 'p2']
    refTs['pyrzin']=200
    refdGs['pyrzin']=[0.000, -0.483, -0.689, -1.055, -1.056]
    refddGs['pyrzin']=[0.005, 0.005, 0.005, 0.005, 0.005]
    refdUs['pyrzin'] = [0.000, -0.176, -0.695, -1.246, 0.951]
    absolutedUs['pyrzin'] = [-14.299, -14.475, -14.993, -15.545, -13.348]

    #Qopbed
    Polys['qopbed']=['p1', 'p3']#, 'p2']
    refTs['qopbed']=200
    refdGs['qopbed']=[0.000, 0.878] #1.699]
    refddGs['qopbed']=[0.005, 0.005]
    refdUs['qopbed'] = [0.000, 1.235]# 2.778]
    absolutedUs['qopbed'] = [-41.135, -39.900] # -38.357]

    #Zzzvye
    Polys['zzzvye']=['p3', 'p1']
    refTs['zzzvye']=50
    refdGs['zzzvye']=[0.000, 0.190]
    refddGs['zzzvye']=[0.005, 0.005]
    refdUs['zzzvye'] = [0.000, 0.351]
    absolutedUs['zzzvye'] = [11.103, 11.454]

    #Resora
    Polys['resora']=['p1', 'p2']
    refTs['resora']=200
    refdGs['resora']=[0.000, -0.337]
    refddGs['resora']=[0.005, 0.005]
    refdUs['resora'] = [0.000, -0.278]
    absolutedUs['resora'] = [-22.506, -22.784]

    #Zzzpro
    Polys['zzzpro']=['p1', 'p3']#, 'p2']
    refTs['zzzpro']=200
    refdGs['zzzpro']=[0.000, 0.223]#, 0.409]
    refddGs['zzzpro']=[0.01, 0.01]#, 0.01]
    refdUs['zzzpro'] = [0.000, 0.258]#, 0.480]
    absolutedUs['zzzpro'] = [-19.393, -19.135]#, -18.913]

    #Zzzpus
    Polys['zzzpus']=['p6', 'p2', 'p1']#, 'p3', 'p7']
    refTs['zzzpus']=200
    refdGs['zzzpus']=[0.000, -1.713, -1.713]#, -1.713, -1.713]
    refddGs['zzzpus']=[0.00, 0.028, 0.02]#, 0.02, 0.02]
    refdUs['zzzpus'] = [0.000, -2.226, 1.4102]#, -2.249, -3.067]
    absolutedUs['zzzpus'] = [-84.917, -87.143, -83.507]#, -87.166, -87.985]

    ##Melfit
    #Polys['melfit']=['p7', 'p6']#, 'p8']
    #refTs['melfit']=200
    #refdGs['melfit']=[0.000, 0.52]#, 0.52]
    #refddGs['melfit']=[0.00, 0.015]#, 0.02]
    #refdUs['melfit'] = [0.000, 0.596]#, 1.187]
    #absolutedUs['melfit'] = [-12.456, -11.818]#, -11.268]

    #Melfit
    Polys['melfit']=['p7', 'p6', 'p8', 'p1', 'p5']
    refTs['melfit']=200
    refdGs['melfit']=[0.000, 0.52, 0.52, 0.52, 0.52]
    refddGs['melfit']=[0.02, 0.02, 0.02, 0.02, 0.02]
    refdUs['melfit'] = [0.000, 0.637, 1.187, 1.966, 0.927]
    absolutedUs['melfit'] = [-12.456, -11.818, -11.268, -10.489, -11.529]



    ##Melfit
    #Polys['melfit']=['p7', 'p6', 'p1', 'p5']#, 'p2', 'p3']#, 'p4']
    #refTs['melfit']=200
    #refdGs['melfit']=[0.000, 0.279, 1.581, 0.542]#, -0.962, 1.361]#, -0.679]
    #refddGs['melfit']=[0.02, 0.02, 0.02, 0.02]#, 0.02, 0.02]#, 0.02]
    #refdUs['melfit'] = [0.000, 0.279, 1.581, 0.542]#, -0.962, 1.361]#, -0.679]
    #absolutedUs['melfit'] = [-12.070, -11.792, -10.489, -11.529]#, -13.032, -10.709]#, -12.749]

    ##Bedmig
    #Polys['bedmig']=['p1', 'p2', 'p6', 'p4', 'p7', 'p3', 'p5']
    #refTs['bedmig']=200
    #refdGs['bedmig']=[0.000, 1.678, -0.372, -2.878, 1.955, -0.322, -1.146]
    #refddGs['bedmig']=[0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
    #refdUs['bedmig'] = [0.000, 3.803, -0.649, -1.217, 2.772, 2.777, 5.638]
    #absolutedUs['bedmig'] = [-92.82, -89.01, -93.46, -94.03, -90.03, -90.0378, -87.1763]

    #Bedmig
    Polys['bedmig']=['p1', 'p7']
    refTs['bedmig']=200
    refdGs['bedmig']=[0.000, 1.955]
    refddGs['bedmig']=[0.00, 0.016]
    refdUs['bedmig'] = [0.000, 2.789]
    absolutedUs['bedmig'] = [-92.814, -90.025]

    ##Bedmig
    #Polys['bedmig']=['p5']
    #refTs['bedmig']=200
    #refdGs['bedmig']=[0.000]
    #refddGs['bedmig']=[0.00]
    #refdUs['bedmig'] = [0.000]
    #absolutedUs['bedmig'] = [-87.1763]


    ##Bedmig
    #Polys['bedmig']=['p1', 'p7', 'p3', 'p5', 'p6']
    #refTs['bedmig']=200
    #refdGs['bedmig']=[0.000, 1.955, -0.322, -1.146, -0.372]
    #refddGs['bedmig']=[0.02, 0.02, 0.02, 0.02, 0.02]
    #refdUs['bedmig'] = [0.000, 2.772, 2.777, 5.638, -0.649]
    #absolutedUs['bedmig'] = [-92.82, -90.03, -90.0378, -87.1763, -93.46]

    #Cafine
    Polys['cafine']=['p1', 'p2']
    refTs['cafine']=200
    refdGs['cafine']=[0.000, 1.846]
    refddGs['cafine']=[0.005, 0.005]
    refdUs['cafine'] = [0.000, 2.375]
    absolutedUs['cafine'] = [-28.737, -26.362]

    
    SystemNAME="OPLSAA"

#Amoeba Reweighting
if potential == "amoeba09":
    #Potentials=['amoeba09','designeda']#, 'designeda', 'designeda', 'designeda', 'designeda', 'designeda']
    Potentials=['amoeba09']
    PotNAME='AMO'
    Charges=['0.1150', '0.1150', '0.0700', '0.0850', '0.1000', '0.1150', '0.1300', '0.1450', '0.1600', '0.1750', '0.1900', '0.1000', '0.0900', '0.0800', '0.0700','0.0600','0.0500','0.0400','0.0300','0.0200','0.0100']
    Chargenames=['', '', 'C00700', 'C00850', 'C01000', 'C01150', 'C01300', 'C01450', 'C01600', 'C01750', 'C01900', 'C00900', 'C00800', 'C00700', 'C00600', 'C00500', 'C00400', 'C00300', 'C00200', 'C00100']
    PotNAMES=['AMO', 'DESA']
    SimNAMES=['TIN', 'GRO']
    #PotNAMES=['AMO']
    #SimNAMES=['TIN']
    refPot=0
    Temperatures=numpy.array([30,40,50,60,70,80,90,100,110,130,140,150,160,170,180,190,200,210,220,230,240,250])
    Pressures=numpy.ones(len(Temperatures),int);
    #Temperatures=numpy.array([60,100,140,200])
    
    #Benzene
    Polys['benzene']=['p1', 'p2', 'p3']
    refTs['benzene']=200
    refdGs['benzene']=[0.000, -0.362, -0.291]
    #refdGs['benzene']=[0.000, 0.390, 0.330]
    refddGs['benzene']=[0.000, 0.031, 0.032]
    refdUs['benzene'] = [0.000, 0.781, 0.643]
    absolutedUs['benzene'] = [-3.587, -2.807, -2.944]
    
    #Formam
    Polys['formam']=['p1', 'p2']
    refTs['formam']=100
    refdGs['formam']=[0.000, 0.293]
    refddGs['formam']=[0.000, 0.020]
    refdUs['formam'] = [0.000, 0.487]
    absolutedUs['formam'] = [-26.237, -25.750]
    
    #Acetac
    Polys['acetac']=['p1', 'p2']
    refTs['acetac']=50
    refdGs['acetac']=[0.000, -0.106]
    refddGs['acetac']=[0.02, 0.02]
    refdUs['acetac'] = [0.000, -0.136]
    absolutedUs['acetac'] = [-34.028, -34.163]
    
    #Imazol
    Polys['imazol']=['p1', 'p2']
    refTs['imazol']=200
    refdGs['imazol']=[0.000, 0.125]
    refddGs['imazol']=[0.020, 0.020]
    refdUs['imazol'] = [0.000, 0.201]
    absolutedUs['imazol'] = [-23.186, -22.985]

    #Glycin
    Polys['glycin']=['p1', 'p2', 'p3']
    refTs['glycin']=10
    refdGs['glycin']=[0.000, 0.163, 1.222]
    refddGs['glycin']=[0.000, 0.000, 0.000]
    refdUs['glycin'] = [0.000, 0.163, 1.222]
    absolutedUs['glycin'] = [0.000, 0.000, 0.000]   

 
    SystemName="AMOEBA"

#Gromos
if potential == "gromos":
    Potentials=['gromos54a7']
    PotNAME='GROM'
    Charges=['0.1150', '0.1150', '0.0700', '0.0850', '0.1000', '0.1150', '0.1300', '0.1450', '0.1600', '0.1750', '0.1900', '0.1000', '0.0900', '0.0800', '0.0700','0.0600','0.0500','0.0400','0.0300','0.0200','0.0100']
    Chargenames=['C01150', 'C01150', 'C00700', 'C00850', 'C01000', 'C01150', 'C01300', 'C01450', 'C01600', 'C01750', 'C01900', 'C00900', 'C00800', 'C00700', 'C00600', 'C00500', 'C00400', 'C00300', 'C00200', 'C00100']
    PotNAMES=['GROM']
    SimNAMES=['GRO']
    Temperatures=numpy.array([10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250])
    refPot=0

    #Benzene
    Polys['benzene']=['p1', 'p2', 'p3']
    refTs['benzene']=200
    refdGs['benzene']=[0.000, 0.251, 0.345]
    refddGs['benzene']=[0.000, 0.019, 0.018]
    refdUs['benzene'] = [0.000, 0.314, 0.328]
    absolutedUs['benzene'] = [-10.184, -9.869, -9.853]

SystemName="GROMOS"

Polymorphs=Polys[molecule]
refT=refTs[molecule]
refdG=refdGs[molecule]
refddG=refddGs[molecule]
refdU=refdUs[molecule]
absolutedU=absolutedUs[molecule]
#refdG=[0.000, 0.390, 0.331]


Colors=['b', 'g', 'r','m','c','y','k','g','g','g','g','g','g','g','g','g','g','g','g'];

if (options.plot):
    import matplotlib # for making plots, version 'matplotlib-1.1.0-1'; errors may pop up when using earlier versions
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.font_manager import FontProperties as FP
    font = {'family' : 'normal',
	    'weight' : 'normal',
	    'size'   : 14}
    matplotlib.rc('font', **font)

#=============================================================================================
# ENSURE THAT USER INPUTS ARE SENSIBLE
#=============================================================================================
#Pressure
if Pressure < 0:
    print "Invalid Pressure: " + str(Pressure)
    sys.exit()

#SIMULATION
if (options.simulation == 'gromacs'):
    simulation_short = 'GRO'
elif (options.simulation == 'tinker'):
    simulation_short = 'TIN'
else:
    print "Invalid Simulation Package: " + str(options.simulation)
    sys.exit()

#ENSEMBLE
if ensemble != "NVE" and ensemble != "NVT" and ensemble != "NPT":
    print "Invalid Ensemble"
    print "Supported Ensembles: NVE NVT NPT"
    sys.exit()

#=============================================================================================
# FORMAT INPUTS
#=============================================================================================
#TEMPERATURE
Tnames=[]
refk=-1
for k,temp in enumerate(Temperatures):
    if temp < 10:
        Tnames.append("00" + str(temp)+"K")
    elif temp < 100:
        Tnames.append("0" + str(temp)+"K")
    else:
        Tnames.append(str(temp)+"K")
    if temp == refT and refk==-1:
	refk=k+refPot*len(Temperatures)

#PRESSURE	
Pnames=[]
for k,press in enumerate(Pressures):
    if press < 10:
	    Pnames.append("_00" + str(press) + "P")
    elif press < 100:
	    Pnames.append("_0" + str(press) + "P")
    else:
	    Pnames.append("_" + str(press) + "P")

#OPTIONAL HINGE
# Format the hinge
if hinge == "DefaultHinge":
    hinge = ""
#else:
    #hinge = "_" + hinge



#=============================================================================================
# READ IN RAW DATA
#=============================================================================================
# Constants.
kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184) # Boltzmann constant in kcal/mol/K

# Parameters
K = len(Potentials)*len(Temperatures)  # How many states?
N_max = 3000 # maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
beta_k = 1.0/(kB*Temperatures)   # beta factor for the different temperatures
beta_k = numpy.tile(beta_k,(1,len(Potentials)))[0]
nm_to_M = 10**-9   #Conversion from angstroms into meters
Bar_to_Pa = 100000  #Conversion from bar to pascals
kJ_to_kcal = 0.2390057 #Conversion from kJ to kcal
Na = 6.022*10**23   #Avogadros numbers
convert_units = kJ_to_kcal*numpy.ones([len(Potentials), len(Potentials)],float) #Convert all energies to kcal/mol
convert_units_V = Bar_to_Pa * nm_to_M**3 * 0.001 * kJ_to_kcal * Na * numpy.ones([len(Potentials), len(Potentials)],float) #Convert all pV to kcal/mol
#convert_units[0,0] = J_to_kcal
ignore_symbols = ['#', '@', '@TYPE', 'STEP', '=====================', 'TIME(ps)'];
Efficiency = numpy.zeros(K, float) #This is the sampling efficiency for each potential in each combination of potentials

# Allocate storage for simulation data
u_pklnT = numpy.zeros([len(Polymorphs),K, K, N_max, 20]) # u_pklnT is the complete matrix of all energy data. 'p' is the polymorph, 'k' is the sampled state, 'l' is
								    # the evaluated state, 'n' is the sample number,  and T is the energy term
N_k = numpy.zeros(K,numpy.int32) # N_k[k] is the total number of snapshots from alchemical state k
Terms_l = [] # Terms_l is the list of all energy terms for the 'l' state 
dA = numpy.zeros([len(Polymorphs),spacing+1,K],float) # dA[p,i,k] is the free energy between potential 0 and state k for spacing i in polymorph p
ddA = numpy.zeros([len(Polymorphs),spacing+1,K],float) # ddA[p,i,k] is the uncertainty in the free energy between potential 0 and state k for spacing i in polymorph p
dG = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # dG[p,i,t] is the free energy between polymorph 1 and polymorph p for spacing i and temperature t
ddG = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # ddG[p,i,t] is the uncertanity in the free energy between polymorph 1 and polymorph p for spacing i and temperature t
dS = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # dS[p,i,t] is the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
ddS = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # ddS[p,i,t] is the uncertanity in the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
dS_mbar = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # dS[p,i,t] is the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
ddS_mbar = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # ddS[p,i,t] is the uncertanity in the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
H_mbar = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # dS[p,i,t] is the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
dH_mbar = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # ddS[p,i,t] is the uncertanity in the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
O_pij = numpy.zeros([len(Polymorphs),len(Temperatures),len(Temperatures)]) # O_pij[p,i,j] is the overlap within polymorph p between temperature state i and temperature state j
G_simp = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # dG[p,i,t] is the free energy between polymorph 1 and polymorph p for spacing i and temperature t
ddG_simp = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures)]) # ddG[p,i,t] is the uncertanity in the free energy between polymorph 1 and polymorph p for spacing i and temperature t
dU = numpy.zeros([len(Polymorphs),len(Temperatures)]) # dG[p,t] is the average energy of polymorph p at temperature t
ddU = numpy.zeros([len(Polymorphs),len(Temperatures)]) # dG[p,t] is the standard deviation of the energy of polymorph p at temperature t
u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of configuration n from potential k in potential l
V_pkn = numpy.zeros([len(Polymorphs), len(Temperatures), N_max], float) # V_pkn is the volume of configuration n of polymorph p at temperature k
V_avg = numpy.zeros([len(Polymorphs), len(Temperatures)], float) # V_avg is the average volume of polymorph p at temperature k
ddV_avg = numpy.zeros([len(Polymorphs), len(Temperatures)], float) # V_avg is the standard deviation of the volume of polymorph p at temperature k

#Cycle through all polymorphs
for p,polymorph in enumerate(Polymorphs):
    #Cycle through all sampled potentials
    for i,potential_k in enumerate(Potentials):
    
        #set up the state k hinges for molecules, charges, and alchemical transformations
        if SimNAMES[i] == "GRO":
	    if Chargenames[i] == "":
		charge_hinge_k=""
	    else:
                charge_hinge_k="_"+Chargenames[i]
	    charge_hinge_k=""
            alchemical_hinge="_000L_100G"
    	    if Molecules == Independent:
    	        Molname = str(Independent) + '_'
    	    else:
    	        Molname = str(Molecules) + '_' + str(Independent) + 'ind_'
        elif SimNAMES[i] == "TIN":
            charge_hinge_k=""
            alchemical_hinge=""
    	    Molname = str(Independent) + '_'
    
        for t,Tname in enumerate(Tnames):
    	    k=len(Temperatures)*i+t
	    Pname=Pnames[t]
    	    #Cycle through all evaluated potentials
    	    for j,potential_l in enumerate(Potentials):
    	        l=len(Temperatures)*j
    
    	        #set up the state l hinges for molecules, charges, and alchemical transformations	
    	        if SimNAMES[j] == "GRO":
		    if Chargenames[i] == "":
                	charge_hinge_l=""
            	    else:
                	charge_hinge_l="_"+Chargenames[j]
    	        elif SimNAMES[j] == "TIN":
    	    	    charge_hinge_l=""
		#dirpath = '/home/ecd4bd/crystals/NMA/finishedJobs/feb1GROM/' + molecule + '_' + SimNAMES[i] + '_' + PotNAMES[i] + '_' + polymorph + '_' + Molname + Tname + charge_hinge_k + alchemical_hinge + Pname + hinge + '/'
    	        dirpath = '/oldhome/ecd4bd/finishedJobs_archive/Reweighted/' + directory + '/' + molecule + '_' + SimNAMES[i] + '_' + PotNAMES[i] + '_' + polymorph + '_' + Molname + Tname + charge_hinge_k + alchemical_hinge + Pname + hinge + '/'
    	        #fname_energy = dirpath + 'potenergy.xvg'
		fname_energy = dirpath + potential_l + charge_hinge_l + '_energy.xvg'    
    	        fname_volume = dirpath + 'volume.xvg'
 
    	        #Read in the energy file data for state k evaluated in state l using HARVIST
    	        if k==0:
    	    	    Terms_l.append(Harvist.CreateTerms(fname_energy))
    	        N=Harvist.DetermineLength(fname_energy,ignoreframes=ignoreframes)
    	        N_k[k]=N
    	        u_pklnT[p,k,l,:N,:len(Terms_l[j])] = Harvist.GrabTerms(fname_energy,Terms_l[j],ignoreframes=ignoreframes)[0]
    	        u_kln[k,l,:N] = Harvist.GrabTerms(fname_energy,['ENERGY', 'Potential'], ignoreframes=ignoreframes)[0][:,0]
    	        #Now set these energies over all temperatures
    	        u_pklnT[p,k,l:(l+len(Temperatures)),:N,:len(Terms_l[j])] = u_pklnT[p,k,l,:N,:len(Terms_l[j])]
    	        u_kln[k,l:(l+len(Temperatures)),:N] = u_kln[k,l,:N]
		#Now read in the volumes and add them to the u_kln matrix
		V_pkn[p,t,:N] = Harvist.GrabTerms(fname_volume, ['VOLUME', 'Volume'], ignoreframes=ignoreframes)[0][:N,0]
		#pdb.set_trace()
		#u_kln[k,l:(l+len(Temperatures)),:N] += Pressures.reshape(len(Temperatures),1)*V_pkn[p,t,:N].reshape(1,N)* Bar_to_Pa * nm_to_M**3 * 0.001 * Na
    
    print "Start1"
    #Convert all units to kcal
    u_pklnT[p,:,:,:] *= kJ_to_kcal
    u_kln *= kJ_to_kcal
    
    print "Start2"
    #If this was already in kcal or already fully independent, revert
    for j in range(len(Potentials)):
        if Potentials[j][:6]=="amoeba":
    	    u_pklnT[p,:,j*len(Temperatures):(j+1)*len(Temperatures),:,:] /= kJ_to_kcal
    	    u_kln[:,j*len(Temperatures):(j+1)*len(Temperatures),:] /= kJ_to_kcal
    
    print "Start3"
    #Remove dependent molecules
    for j in range(len(Potentials)):
        if Potentials[j][:6]!="amoeba":
    	    u_pklnT[p,:,j*len(Temperatures):(j+1)*len(Temperatures),:,:] *= float(Independent)/Molecules
    	    u_kln[:,j*len(Temperatures):(j+1)*len(Temperatures),:] *= float(Independent)/Molecules
    print "Start4"

    #Now average together the energies and volumes at each state
    for t in range(len(Temperatures)):
	dU[p,t]=numpy.average(u_kln[t,t,:N_k[t]])/float(Independent)
	ddU[p,t]=numpy.std(u_kln[t,t,:N_k[t]])/N_k[t]**0.5/float(Independent)
	V_avg[p,t]=numpy.average(V_pkn[p,t,:N_k[t]])/float(Independent)
	ddV_avg[p,t]=numpy.std(V_pkn[p,t,:N_k[t]])/N_k[t]**0.5/float(Independent)

    print "Start5"
    # convert to nondimensional units from kcal/mol
    for k,beta in enumerate(beta_k):
        u_kln[:,k,:] *=  beta  
 
    u_kln_save = u_kln.copy()
    N_k_save = N_k.copy()
    print "End!"

    print "Number of retained samples"
    print N_k

    #Now create the full N_k matrix including the roll-backs as well as the free energy container
    N_k_matrix = numpy.zeros([spacing+1,K],numpy.int32) # N_k_matrix[i,k] is the total number of snapshots from alchemical state k using in spacing i
    for i in range(spacing+1):
        N_k_matrix[i,:] = N_k_save.copy()
        N_k_matrix[i,0:len(Temperatures)] = N_k_matrix[i,0:len(Temperatures)]*float(i)/float(spacing)


    
    #=============================================================================================
    # COMPUTE FREE ENERGY DIFFERENCE USING MBAR FOR EACH SPACING
    #=============================================================================================
    
    for i in range(spacing+1):
   	if i==0 and len(Potentials)==1: 
	    continue
        # Initialize MBAR.
        print "Running MBAR..."
        # generate the weights of each of the umbrella set
        #mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive', use_optimized=False)
        mbar = pymbar.MBAR(u_kln, N_k_matrix[i,:], verbose = True)
        print "MBAR Converged..."
    
        # extract self-consistent weights and uncertainties
        (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()
	# extract entropy
	[Delta_f_ij, dDelta_f_ij, Delta_u_ij, dDelta_u_ij, Delta_s_ij, dDelta_s_ij] = mbar.computeEntropyAndEnthalpy()
    
        print "Free Energies Optained..."
    
        #Store the dimensionless results in the dA container
        dA[p,i,:] = df_i[refk]
	dH_mbar[p,i,:] = Delta_u_ij[0]
	dS_mbar[p,i,:] = Delta_s_ij[0]
	ddS_mbar[p,i,:] = dDelta_s_ij[0]
	print dA
    
    #=============================================================================================
    # COMPUTE UNCERTAINTY USING MBAR
    #=============================================================================================
    g_k = numpy.zeros([K])
    N_k_matrix_save = N_k_matrix.copy()
    for i in range(spacing+1):
	if i==0 and len(Potentials)==1:
            continue
        for k in range(K):
    	    #subsample correlated data - for now, use energy from current state
    	    if N_k_matrix[i,k]>0:
    	        print N_k_matrix[i,k]
    	        #print str(i) + ' ' + str(k)
    	        #print u_kln_save[k,k,0:100]
    	        g_k[k] = timeseries.statisticalInefficiency(u_kln_save[k,k,0:100])
    	        print "Correlation time for phase (%s), sampled state %d is %10.3f" % (phase,k,g_k[k])
    	        # subsample the data to get statistically uncorrelated data
    	        indices = numpy.array(timeseries.subsampleCorrelatedData(u_kln_save[k, k, 0:N_k_matrix[i,k]], g=g_k[k]))  # subsample
    	        N_k_matrix[i,k] = len(indices)
    	        u_kln[k,:,0:N_k_matrix[i,k]] = u_kln_save[k,:,indices].transpose()  # not sure why we have to transpose
        
        print "Number of retained samples"
        print N_k
    
        print "Running MBAR..."
    
        # generate the weights of each state
        mbar = pymbar.MBAR(u_kln, N_k_matrix[i,:], verbose = True)
        print "MBAR Converged..." 

        # extract self-consistent weights and uncertainties
        (df_u, ddf_u, theta_u) = mbar.getFreeEnergyDifferences()

	# calculate the overlap it necessary
	if len(Temperatures)==2:
	    O_pij[p,:,:] = mbar.computeOverlap()[2]

    
        # testing
        weights_in_gromos = numpy.zeros(K,float)
        for k in range(K):
    	    w = numpy.exp(mbar.Log_W_nk[:,k])
    	    print "max weight in state %d is %12.7f" % (k,numpy.max(w))
    	    # using Kish (1965) formula.
    	    # effective # of samples =  (\sum_{i=1}^N w_i)^2 / \sum_{i=1}^N w_i^2
    	    #                        =  (\sum_{i=1}^N w_i^2)^-1
    	    neff = 1/numpy.sum(w**2)
    	    #Store the effective number of samples in the target potential with no sampling
    	    if i==0 and k==0:
    	        neff_target=neff
    	    print "Effective number of sample in state %d is %10.3f" % (k,neff)
    	    print "Efficiency for state %d is %d/%d = %10.4f" % (k,neff,len(w),neff/len(w))
    	    Efficiency[k] = neff/len(w) #Store the efficiency
    	    w_0 = numpy.exp(mbar.Log_W_nk[:,0]) #Weights in gromos
    	    initial_configs = numpy.sum(N_k[0:k])
    	    final_configs = numpy.sum(N_k[0:k+1])
    	    print "Total weight in gromos " + str(numpy.sum(w_0[initial_configs:final_configs]))
    	    weights_in_gromos[k] = numpy.sum(w_0[initial_configs:final_configs])
    
        # Write out free energy differences
        print "Free Energy Difference (in units of kcal/mol)"
        for k in range(K):
    	    print "%8.3f %8.3f" % (-df_i[k,0], ddf_u[k,0])   
    
        #Store the dimensionless results in the ddA container
        ddA[p,i,:] = ddf_u[refk]
   
     
    ##Bootstrap the reduced free energies
    #indexVect = numpy.zeros([2,len(N_k)], numpy.float)
    #indexVect[0,:]=refk
    #indexVect[1,:] = numpy.arange(len(N_k)) 
    #datafile = 'BootstrapData/BootstrapData_' + molecule + '_' + polymorph + '_' + potential + '_dfvsT'
    #stdfile = 'BootstrapData/BootstrapStd_' + molecule + '_' + polymorph + '_' + potential + '_dfvsT'

    #if not os.path.isfile(stdfile+'.txt'):
    #    MBARBootstrap.runMBARBootstrap(u_kln, N_k, numpy.ones(len(N_k),float), 1, indexVect, datafile, stdfile, 200)
    
    #ddA[p,spacing,:]=MBARBootstrap.ExtractBootstrap(stdfile+'.txt')
    #if len(Potentials) > 1:
    #	datafile = 'BootstrapData_' + molecule + '_' + polymorph + '_' + potential + '_indirect_dGvsT'
    #    stdfile = 'BootstrapStd_' + molecule + '_' + polymorph + '_' + potential + '_indirect_dGvsT'
    #	if not os.path.isfile(stdfile+'.txt'):
    #	    MBARBootstrap.runMBARBootstrap(u_kln, N_k, numpy.ones(len(N_k),float), 1, indexVect, datafile, stdfile, 200)
    #    ddA[p,0,:]=MBARBootstrap.ExtractBootstrap(stdfile+'.txt')
     

#os.exit()
#Check the overlap it necessary
if len(Temperatures)==2:
    print "Overlap:"
    print O_pij
    pdb.set_trace()

#=============================================================================================
# FINALIZE THE RELATIVE FREE ENERGY AND ENTROPY
#=============================================================================================
dT=Temperatures[1]-Temperatures[0]

#Combine together all bootstrap data
BootstrapNum = 200
i=spacing
dA_Bootstrap0 = numpy.zeros(len(Temperatures), float) #dA_Bootstrap0 is the reduced free energy for the reference polymorph at each temperature t at a particular bootstrap iteration n
dA_Bootstrapp = numpy.zeros(len(Temperatures), float) #dA_Bootstrapp is the reduced free energy for polymorph p  at each temperature t at a particular bootstrap iteration n
dG_Bootstrap = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures), BootstrapNum]) # dG_Bootstrap[p,i,t,n] is the free energy between the reference polymorph and polymorph p for spacing i and temperature t at bootstrap iteration n
dS_Bootstrap = numpy.zeros([len(Polymorphs),spacing+1,len(Temperatures), BootstrapNum]) # dS_Bootstrap[p,i,t,n] is the entropy difference between the reference polymorph and polymorph p for spacing i and temperature t at bootstrap iteration n

if refT < 10:
    refTname="00" + str(refT)+"K"
elif refT < 100:
    refTname="0" + str(refT)+"K"
else:
    refTname=str(refT)+"K"

for p in range(len(Polymorphs)):
    if p==0:
	continue
	

    fnameBootstrapStddGvT = 'BootstrapData/BootstrapStd_' + molecule + '_' + Polymorphs[p] + '_' + potential + '_dGvsT.txt'
    fnameBootstrapStddSvT = 'BootstrapData/BootstrapStd_' + molecule + '_' + Polymorphs[p] + '_' + potential + '_dSvsT.txt'
    
    #if not os.path.isfile(fnameBootstrapStddGvT):
    #
    #    for n in range(BootstrapNum):
    #        dA_Bootstrap0[:] = MBARBootstrap.ExtractBootstrapLine( 'BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[0] + '_' + potential + '_dfvsT.txt', n)
    #        dA_Bootstrapp[:] = MBARBootstrap.ExtractBootstrapLine( 'BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[p] + '_' + potential + '_dfvsT.txt', n)
    #        refdG0_L = MBARBootstrap.ExtractBootstrapLine('BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[0] + '_' + Molname + refTname + '_001P_' + PotNAME + '_L.txt',n)[19]
    #        #refdG0_L = MBARBootstrap.ExtractBootstrapLine('BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[0] + '_' + Molname + refTname + '_001P_' + PotNAME + '_L_50kDihedral.txt',n)[19]
    #        refdG0_G = MBARBootstrap.ExtractBootstrapLine('BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[0] + '_' + Molname + refTname + '_001P_' + PotNAME + '_G.txt',n)[9]
    #        refdGp_L = MBARBootstrap.ExtractBootstrapLine('BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[p] + '_' + Molname + refTname + '_001P_' + PotNAME + '_L.txt',n)[19]
    #        #refdGp_L = MBARBootstrap.ExtractBootstrapLine('BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[p] + '_' + Molname + refTname + '_001P_' + PotNAME + '_L_50kDihedral.txt',n)[19]
    #        refdGp_G = MBARBootstrap.ExtractBootstrapLine('BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[p] + '_' + Molname + refTname + '_001P_' + PotNAME + '_G.txt',n)[9]
    #        #refdG0_G10k = MBARBootstrap.ExtractBootstrapLine('BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[0] + '_' + Molname + refTname + '_001P_' + PotNAME + '_G_50kDihedral2.txt',n)[9]
    #        #refdGp_G10k = MBARBootstrap.ExtractBootstrapLine('BootstrapData/BootstrapData_' + molecule + '_' + Polymorphs[p] + '_' + Molname + refTname + '_001P_' + PotNAME + '_G_50kDihedral2.txt',n)[9]
    #        refPSCP = (refdGp_L + refdGp_G) - (refdG0_L + refdG0_G)
    #        #refPSCP = (refdGp_L - refdGp_G + refdGp_G10k ) - (refdG0_L - refdG0_G + refdG0_G10k) #Flexible Molecules
    #        for t,T in enumerate(Temperatures):
    #            dG_Bootstrap[p,i,t,n] = (dA_Bootstrapp[t]-dA_Bootstrap0[t])/(beta_k[t]*float(Independent)) + float(T)/float(refT)*refPSCP #Compute the free energy difference
    #            dS_Bootstrap[p,i,t,n] = (dU[p,t]-dU[0,t]-dG_Bootstrap[p,i,t,n])/float(T)
    #    
    #    #Now compute the final uncertainties in the entropy and free energy
    #    for t,T in enumerate(Temperatures):
    #        ddG[p,i,t] = numpy.std(dG_Bootstrap[p,i,t,:]) 
    #        ddS[p,i,t] = numpy.std(dS_Bootstrap[p,i,t,:])
    #    #Output the final bootstrapped uncertainties to a MBARBootstrap format file that can be read later
    #    with open(fnameBootstrapStddGvT, "a") as myfile:
    #        myfile.write('Final bootstrapped uncertainties for dGvT')
    #        myfile.write('\n')
    #        myfile.write('Avg:  ' + ' '.join(str(cell) for cell in ddG[p,i,:]))
    #    with open(fnameBootstrapStddSvT, "a") as myfile:
    #        myfile.write('Final bootstrapped uncertainties for dSvT')
    #        myfile.write('\n')
    #        myfile.write('Avg:  ' + ' '.join(str(cell) for cell in ddS[p,i,:]))


    #Load in the free energy differences
    #ddG[p,i,:]=numpy.array(MBARBootstrap.ExtractBootstrap(fnameBootstrapStddGvT))
    #ddS[p,i,:]=numpy.array(MBARBootstrap.ExtractBootstrap(fnameBootstrapStddSvT))





for i in range(spacing+1):
    for t,T in enumerate(Temperatures):
	for p in range(len(Polymorphs)):
	    #if p==0:
	    #	continue
	    #dG[p,i,t] = ((dA[p,i,refPot*len(Temperatures)+t]+(dA[p,i,t] - dA[p,i,refPot*len(Temperatures)+t])          dA[0,i,refPot*len(temperatures)+t]) + (dA[p,i,refPot*len(Temperatures)+t] - dA[0,i,refPot*len(Temperatures)+t]))/(beta_k[t]*float(Independent)) + T/refT*refdG[p]
	    dG[p,i,t] = (dA[p,i,t]-dA[0,i,t])/(beta_k[t]*float(Independent)) + float(T)/float(refT)*refdG[p]
	    ddG[p,i,t] = ((ddA[p,i,t]**2 + ddA[0,i,t]**2)/(beta_k[t]*float(Independent))**2 + float(T)/float(refT)*float(refddG[p])**2)**0.5
	    #ddG[p,i,t] = ((ddA[p,i,t]**2 + ddA[0,i,t]**2 + ddA[p,i,t]**2 + ddA[0,i,t]**2)/(beta_k[t]*float(Independent))**2 + T/refT*refddG[p]**2)**0.5
	    
	    #dS_mbar[p,i,t] /= beta_k[t]
	    
	    if p==0:
		continue
	    dS[p,i,t] = (dU[p,t]-dU[0,t]-dG[p,i,t])/float(T)
	    ddS[p,i,t] = (ddU[p,t]**2 + ddU[p,t]**2 + ddG[p,i,t]**2)**0.5/float(T)

print "Polymorph Free Energy:"
for p in range(len(Polymorphs)):
    print "%8.3f %8.3f" % (dG[p,spacing,len(Temperatures)-1], ddG[p,spacing,len(Temperatures)-1])


##Now print the y-intersept using both the 10K point and the 20K point
#for p in range(len(Polymorphs)):
#    if p==0:
#	continue
#    print "Intercept extrapolated from 20K:"
#    print str(dG[p,spacing,Temperatures==20] + 20.0*dS[p,spacing,Temperatures==20]) + " +/- " + str(20.0*ddS[p,spacing,Temperatures==20]) 
#    print "Intercept extrapolated from 10K:"
#    print str(dG[p,spacing,Temperatures==10] + 10.0*dS[p,spacing,Temperatures==10]) + " +/- " + str(10.0*ddS[p,spacing,Temperatures==10])
#    print "Actual Intercept:"
#    print refdU[p]
#
##Now print the free energies and uncertainties at 100K, 200K, 300K for closure comparison
#print "Polymorph " + Polymorphs[1] + " Free energy at 100K:"
#print "%8.3f %8.3f" % (dG[1,spacing,Temperatures==100], ddG[1,spacing,Temperatures==100])
#print "Polymorph " + Polymorphs[1] + " Free energy at 200K:"
#print "%8.3f %8.3f" % (dG[1,spacing,Temperatures==200], ddG[1,spacing,Temperatures==200])
#print "Polymorph " + Polymorphs[1] + " Free energy at 300K:"
#print "%8.3f %8.3f" % (dG[1,spacing,Temperatures==300], ddG[1,spacing,Temperatures==300]) 

#=============================================================================================
# PLOT THE RELATIVE FREE ENERGY VS TEMPERATURE
#=============================================================================================
f,a = plt.subplots(1,1)
xlabel='Temperature (K)'
ylabel='Relative Free Energy (kcal/mol)'

#dG_AMO_ONLY=numpy.load('dG_AMOOnly.npy')
Temperatures2=numpy.array([0] + [j for j in Temperatures])
Pressures2=numpy.array([1] + [j for j in Pressures])
#Temperatures=numpy.concatenate(0, Temperatures)
dG2=numpy.insert(dG,0,numpy.transpose(numpy.tile(refdU,(1,1))),axis=2)
ddG2=numpy.insert(ddG,0,numpy.zeros([len(Polymorphs),spacing+1]),axis=2)
dS2=numpy.insert(dS,0,numpy.zeros([len(Polymorphs),spacing+1]),axis=2)
ddS2=numpy.insert(ddS,0,numpy.zeros([len(Polymorphs),spacing+1]),axis=2)

##Compute the entropy assuming a linear fit
#dS = numpy.zeros(len(dG),float)
#ddS = numpy.zeros(len(dG), float)
#S=numpy.zeros(200, float)
#polyfit = []
#for p in range(len(Polymorphs)):
#    if p==0:
#        polyfit.append([0.0, 0.0])
#        continue
#    polyfit.append(numpy.polyfit(Temperatures2,dG2[p,spacing,:],1))
#    dS[p] = -1.0*float(polyfit[p][0])
#
#    #Error estimate on the entropy
#    Trandom = numpy.zeros(len(Temperatures2))
#    dGrandom=numpy.zeros(len(Temperatures2))
#    for i in range(200):
#	for j in range(len(Temperatures2)):
#	    randomnum = int(random.random()*len(Temperatures2))
#	    Trandom[j] = Temperatures2[randomnum]
#	    dGrandom[j]=dG2[p,spacing,randomnum]
#	S[i]=numpy.polyfit(Trandom,dGrandom,1)[0]
#    ddS[p]=numpy.std(S)
#	   
##convert to cal from kcal
#dS = dS*1000
#ddS = ddS*1000

#Tempertures=numpy.array(Temperatures)
a.errorbar(Temperatures2,dG2[0,0,:],yerr=numpy.zeros(len(dG2[0,0,:]),float),linestyle='--',marker='.', linewidth=2,alpha=0.6,color='b',label=Polymorphs[0])
PlotPress=1 #Pressure to plot the dGvT curve at
Temperatures_P = Temperatures2[Pressures2==PlotPress]
for p in range(len(Polymorphs)):
    if p==0:
	continue
    if len(Potentials)>1:
        a.errorbar(Temperatures2,dG2[p,0,:],yerr=ddG2[p,0,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color=Colors[p],label=Polymorphs[p])
    a.errorbar(Temperatures_P,dG2[p,1,Pressures2==PlotPress],yerr=ddG2[p,1,Pressures2==PlotPress],linestyle='-',marker='.', linewidth=2,alpha=0.6,color=Colors[p],label=Polymorphs[p])
    #a.plot(Temperatures2, numpy.array(Temperatures2*float(polyfit[p][0])+float(polyfit[p][1])), linestyle='--',marker='None', linewidth=1,alpha=1.0,color='k')

a.set_xlabel(xlabel)
a.set_ylabel(ylabel)
#a.errorbar(Temperatures2,dG2[1,0,:],yerr=ddG2[1,0,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color='g',label="Benzene II")
#a.errorbar(Temperatures2,dG2[2,0,:],yerr=ddG2[2,0,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color='r',label="Benzene III")
#a.errorbar(Temperatures,dU[1,:]-dU[0,:],yerr=ddU[1,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color='g',label="Benzene II")
#a.errorbar(Temperatures,dU[2,:]-dU[0,:],yerr=ddU[2,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color='r',label="Benzene III")
#a.errorbar(Temperatures2,dG2[0,spacing,:],yerr=ddG2[0,spacing,:],linestyle='-',marker='.', linewidth=2,alpha=0.6,color='b',label="Benzene I")
#a.errorbar(Temperatures2,dG2[1,spacing,:],yerr=ddG2[1,spacing,:],linestyle='-',marker='.', linewidth=2,alpha=0.6,color='g',label="Benzene II")
#a.errorbar(Temperatures2,dG2[2,spacing,:],yerr=ddG2[2,spacing,:],linestyle='-',marker='.', linewidth=2,alpha=0.6,color='r',label="Benzene III")
#a.errorbar(Temperatures,dG_AMO_ONLY[1,spacing,:],yerr=ddG[1,spacing,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color='g',label="Benzene II")
#a.errorbar(Temperatures,dG_AMO_ONLY[2,spacing,:],yerr=ddG[2,spacing,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color='r',label="Benzene III")

##Save the data for future use.
#filename=potential
NPdir='NPFiles/'
numpy.save(NPdir+'T_' + molecule + '_' + potential, Temperatures_P)
#for p,Poly in enumerate(Polymorphs):
#    numpy.save(NPdir+'dGvT_' + molecule + '_' + Poly + '_' + potential,dG2[p,spacing,Pressures2==PlotPress])
#    numpy.save(NPdir+'ddGvT_' + molecule + '_' + Poly + '_' + potential, ddG2[p,spacing,Pressures2==PlotPress])
#    if len(Potentials) > 1:
#        numpy.save(NPdir+'dGvT_' + molecule + '_' + Poly + '_' + potential + '_indirect',dG2[p,0,:])
#        numpy.save(NPdir+'ddGvT_' + molecule + '_' + Poly + '_' + potential + '_indirect', ddG2[p,0,:])
#    	#if spacing > 1:
#	#    numpy.save(NPdir+'dGvT_' + molecule + '_' + Poly + '_' + potential + '_convergence',dG2[p,:,:])
#        #    numpy.save(NPdir+'ddGvT_' + molecule + '_' + Poly + '_' + potential + '_convergence', ddG2[p,:,:])
#    numpy.save(NPdir+'dS_' + molecule + '_' + Poly + '_' + potential, dS2[p,spacing,:])
#    numpy.save(NPdir+'ddS_' + molecule + '_' + Poly + '_' + potential, ddS2[p,spacing,:])



#=============================================================================================                                                                                                                                              # PLOT THE RELATIVE FREE ENERGY VS PRESSURE                                                                                                                                                                                                 #=============================================================================================
#f,a = plt.subplots(1,1)
#xlabel='Pressure (GPa)'
#ylabel='Relative Free Energy (kcal/mol)'
#
#PlotTemp=200 #Temperature to plot the dGvsP curve at
#Pressure_T = Pressures[Temperatures==PlotTemp]/float(10000)
#for p in range(len(Polymorphs)):
#    if p==0:
#        continue
#    a.errorbar(Pressure_T,dG2[p,1,Temperatures2==PlotTemp],yerr=ddG2[p,1,Temperatures2==PlotTemp],linestyle='-',marker='.', linewidth=2,alpha=0.6,color=Colors[p],label=Polymorphs[p])
#
#a.set_xlabel(xlabel)
#a.set_ylabel(ylabel)
#
##Save the data for future use.
#filename=potential
#NPdir='NPFiles/'
#numpy.save(NPdir+'P_' + molecule + '_' + potential, Pressure_T)
#for p,Poly in enumerate(Polymorphs):
#    numpy.save(NPdir+'dGvP_' + molecule + '_' + Poly + '_' + potential,dG2[p,spacing,Temperatures2==PlotTemp])
#    numpy.save(NPdir+'ddGvP_' + molecule + '_' + Poly + '_' + potential, ddG2[p,spacing,Temperatures2==PlotTemp])
#    if len(Potentials) > 1:
#        numpy.save(NPdir+'dGvT_' + molecule + '_' + Poly + '_' + potential + '_indirect',dG2[p,0,:])
#        numpy.save(NPdir+'ddGvT_' + molecule + '_' + Poly + '_' + potential + '_indirect', ddG2[p,0,:])
#        if spacing > 1:
#            numpy.save(NPdir+'dGvT_' + molecule + '_' + Poly + '_' + potential + '_convergence',dG2[p,:,:])
#            numpy.save(NPdir+'ddGvT_' + molecule + '_' + Poly + '_' + potential + '_convergence', ddG2[p,:,:])
#numpy.save(NPdir+'dS_' + molecule + '_' + potential, dS)
#numpy.save(NPdir+'ddS_' + molecule + '_' + potential, ddS)


#=============================================================================================
# PLOT THE RELATIVE ENTROPY VS TEMPERATURE
#=============================================================================================
f,a = plt.subplots(1,1)
plt.xlabel('Temperature (K)',fontsize=20)
plt.ylabel('Relative Entropy (kcal/molK)',fontsize=20)
for p in range(len(Polymorphs)):
    a.errorbar(Temperatures2,dS2[p,spacing,:],yerr=ddS2[p,spacing,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color=Colors[p],label=Polymorphs[p])

#=============================================================================================
# PLOT THE AVERAGE ENERGY VS TEMPERATURE
#=============================================================================================

f,a = plt.subplots(1,1)
plt.xlabel('Temperature (K)',fontsize=20)
plt.ylabel('Average Energy (kcal/mol)',fontsize=20)
dU2=numpy.insert(dU,0,numpy.transpose(absolutedU),axis=1)
ddU2=numpy.insert(ddU,0,0,axis=1)
for p in range(len(Polymorphs)):
    a.errorbar(Temperatures2,dU2[p,:],yerr=ddU2[p,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color=Colors[p],label=Polymorphs[p])
plt.legend(loc='upper left')
NPdir='NPFiles/'
for p,Poly in enumerate(Polymorphs):
    numpy.save(NPdir+'UvT_' + molecule + '_' + Poly + '_' + potential,dU2[p,:])

#=============================================================================================
# PLOT THE AVERAGE BOX VOLUME VS TEMPERATURE
#=============================================================================================

f,a = plt.subplots(1,1)
plt.xlabel('Temperature (K)',fontsize=20)
plt.ylabel('Average Volume (nm^3/mol)',fontsize=20)
#dU2=numpy.insert(dU,0,numpy.transpose(absolutedU),axis=1)
#ddU2=numpy.insert(ddU,0,0,axis=1)
for p in range(len(Polymorphs)):
    a.errorbar(Temperatures,V_avg[p,:],yerr=ddV_avg[p,:],linestyle='--',marker='.', linewidth=2,alpha=0.6,color=Colors[p],label=Polymorphs[p])

NPdir='NPFiles/'
for p,Poly in enumerate(Polymorphs):
    numpy.save(NPdir+'VvT_' + molecule + '_' + Poly + '_' + potential,V_avg[p,:])

#=============================================================================================
# PLOT THE DIFFERENCE IN AVERAGE ENERGY VS TEMPERATURE
#=============================================================================================

f,a = plt.subplots(1,1)
plt.xlabel('Temperature (K)',fontsize=20)
plt.ylabel('Average Energy Difference (kcal/mol)',fontsize=20)
for p in range(len(Polymorphs)):
    if p==0:
	yerror=numpy.zeros(len(ddU2[0,:]))
    else:
	yerror=[]
	for t in range(len(ddU2[0,:])):
	    yerror.append(float((ddU2[0,t]**2+ddU2[p,t]**2)**0.5))
    a.errorbar(Temperatures2,dU2[p,:]-dU2[0,:],yerr=yerror,linestyle='--',marker='.', linewidth=2,alpha=0.6,color=Colors[p],label=Polymorphs[p])

#Save the data for future use.
NPdir='NPFiles/'
for p,Poly in enumerate(Polymorphs):
    numpy.save(NPdir+'dUvT_' + molecule + '_' + Poly + '_' + potential,dU2[p,:]-dU2[0,:])
    numpy.save(NPdir+'ddUvT_' + molecule + '_' + Poly + '_' + potential, (ddU2[p,:]**2 + ddU2[0,:]**2)**0.5)

plt.show()