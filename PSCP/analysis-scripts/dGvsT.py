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
from optparse import OptionParser # for parsing command-line options
import mdtraj as md
import MBARBootstrap # Bootstrapping algorithm
import random
import os
import usefulFuncs #Useful math functions
import Harvist #Hamiltonian Reweighting Visualization Toolkik
import pdb
import sys
import panedr
import matplotlib
import matplotlib.pyplot as plt

font = {'family': 'normal',
        'weight': 'normal',
        'size': 16}



def old_systems_dictionary(potential, molecule):
    Polys = dict()  # Name of the polymorphs
    refTs = dict()  # Reference temperatures for the PSCP for each system
    refdGs = dict()  # Reference free energies for the PSCP for each system
    refddGs = dict()  # Reference uncertainties for the PSCP for each system
    refdUs = dict()  # Reference lattice minima for each system
    absolutedUs = dict()  # Absolute lattice energy for each system
    
    
    # Oplsaa
    if potential == "oplsaa":
        Potentials = ['oplsaa']
        PotNAME = 'OPLS'
        Charges = ['0.1150', '0.1150', '0.0700', '0.0850', '0.1000', '0.1150', '0.1300', '0.1450', '0.1600', '0.1750',
                   '0.1900', '0.1000', '0.0900', '0.0800', '0.0700', '0.0600', '0.0500', '0.0400', '0.0300', '0.0200',
                   '0.0100']
        Chargenames = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
                       '', '', '', '', '', '', '', '']
        #Chargenames=['C01150', 'C01150', 'C00700', 'C00850', 'C01000', 'C01150', 'C01300', 'C01450', 'C01600', 'C01750', 'C01900', 'C00900', 'C00800', 'C00700', 'C00600', 'C00500', 'C00400', 'C00300', 'C00200', 'C00100']
        PotNAMES=['OPLS']
        SimNAMES=['GRO']
        Temperatures = np.array([50, 100, 200, 300])
    
        #Temperatures=np.array([45,50]) #Overlap Check
        #Temperatures=np.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,70,80,90,100,110,120,130,140])#,80,90,100,110])#Zzzvye
        #Temperatures=np.array([10,15,20,25,30,35,40,45,50,60,70,80,90,100,110])#Zzzvye
        #Temperatures = np.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340])#,360,380,400,420,440,460,480,500])#,320,340,360,380,400]) #benzene
        #Temperatures=np.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340])#,360,380,400,420,440,460,480,500])#,320,340,360,380,400]) #kobfud
        #Temperatures=np.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340,360,380,400,420,440,460,480,500])#,320,340,360,380,400]) #cbmzpn
        #Temperatures=np.array([10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340])#,360,380,400,420,440,460,480,500])#,320,340,360,380,400]) #melfit
        #Temperatures=np.array([10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300])#,320,340,360,380,400,420,440,460,480,500])
        #Temperatures=np.array([100,130,140,160,180,200,220,240,260,280,300,340,360,380,400])
        #Temperatures=np.array([100,120,140,160,180,200]) #MelfitTest
        #ExtraPressures=np.array([5000,15000,35000]) #benzene
        #Temperatures=np.array([40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300]) #imazol
        #ExtraPressures=np.array([5000,10000,16000,20000,26000,31000,36000,40000]) #imazol
        #Temperatures=np.array([20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,220,240,260,280,300,300,300,300,300,300,300,300]) #acetac
        #ExtraPressures=np.array([5000,5000,5000,5000,5000,10000,16000,20000,26000,31000,36000,40000]) #acetac
        #Temperatures=np.array([20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,200,200,200,200,220,240,260,280,300,300,300,300,300,300,300,300]) #formam
        #ExtraPressures=np.array([1500,2000,5000,10000,10000,10000,10000,10000,5000,10000,15000,20000,26000,31000,36000,40000]) #formam
        #Temperatures=np.array([20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,250,270,290,300,320,340,360,380,400,420,440,460,480,500]) #glycin
        ExtraPressures = []
        Pressures = np.ones(len(Temperatures), int)
        Pressures[len(Pressures) - len(ExtraPressures): len(Pressures)] = ExtraPressures
        #Temperatures=np.array([200,200,200,200,200])
        #Pressures=np.array([1,500,1000,1500,2000,2500,3000,3500,4000,4500,5000])
        refPot = 0
        
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
        Polys['zzzpus']=['p6', 'p2'] #, 'p1']#, 'p3', 'p7']
        refTs['zzzpus']=200
        refdGs['zzzpus']=[0.000, -1.713] #, -1.713]#, -1.713, -1.713]
        refddGs['zzzpus']=[0.00, 0.028] #, 0.02]#, 0.02, 0.02]
        refdUs['zzzpus'] = [0.000, -2.226] #, 1.4102]#, -2.249, -3.067]
        absolutedUs['zzzpus'] = [-84.917, -87.143] #, -83.507]#, -87.166, -87.985]
    
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
    
        SystemNAME = "OPLSAA"
    
    # Amoeba Reweighting
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
        Temperatures=np.array([30,40,50,60,70,80,90,100,110,130,140,150,160,170,180,190,200,210,220,230,240,250])
        Pressures=np.ones(len(Temperatures),int);
        #Temperatures=np.array([60,100,140,200])
        
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
    
        SystemName = "AMOEBA"
    
    #Gromos
    if potential == "gromos":
        Potentials=['gromos54a7']
        PotNAME='GROM'
        Charges=['0.1150', '0.1150', '0.0700', '0.0850', '0.1000', '0.1150', '0.1300', '0.1450', '0.1600', '0.1750', '0.1900', '0.1000', '0.0900', '0.0800', '0.0700','0.0600','0.0500','0.0400','0.0300','0.0200','0.0100']
        Chargenames=['C01150', 'C01150', 'C00700', 'C00850', 'C01000', 'C01150', 'C01300', 'C01450', 'C01600', 'C01750', 'C01900', 'C00900', 'C00800', 'C00700', 'C00600', 'C00500', 'C00400', 'C00300', 'C00200', 'C00100']
        PotNAMES=['GROM']
        SimNAMES=['GRO']
        Temperatures=np.array([10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250])
        refPot=0
    
        #Benzene
        Polys['benzene']=['p1', 'p2', 'p3']
        refTs['benzene']=200
        refdGs['benzene']=[0.000, 0.251, 0.345]
        refddGs['benzene']=[0.000, 0.019, 0.018]
        refdUs['benzene'] = [0.000, 0.314, 0.328]
        absolutedUs['benzene'] = [-10.184, -9.869, -9.853]
    
        SystemName = "GROMOS"
    
    Polymorphs = Polys[molecule]
    refT = refTs[molecule]
    refdG = refdGs[molecule]
    refddG = refddGs[molecule]
    refdU = refdUs[molecule]
    absolutedU = absolutedUs[molecule]
    return Polymorphs, refT, refdG, refddG, refdU, absolutedU

def get_potential_info(potential):
    if potential == 'oplsaa':
        SimNAMES=['GRO']
        Chargenames = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
                       '', '', '', '', '', '', '', '', '']
        PotNAMES = ['OPLS']
    elif potential == 'amoeba09':
        SimNAMES=['TIN', 'GRO']
        Chargenames = ['', '', 'C00700', 'C00850', 'C01000', 'C01150', 'C01300', 'C01450', 'C01600', 'C01750', 'C01900',
                       'C00900', 'C00800', 'C00700', 'C00600', 'C00500', 'C00400', 'C00300', 'C00200', 'C00100']
        PotNAMES = ['AMO', 'DESA']
    elif potential == "gromos":
        SimNAMES = ['GRO']
        Chargenames = ['C01150', 'C01150', 'C00700', 'C00850', 'C01000', 'C01150', 'C01300', 'C01450', 'C01600',
                       'C01750', 'C01900', 'C00900', 'C00800', 'C00700', 'C00600', 'C00500', 'C00400', 'C00300',
                       'C00200', 'C00100']
        PotNAMES = ['GROM']
    return SimNAMES, Chargenames, PotNAMES

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

def dGvsT(plot_out=False, Temperatures=np.array([100,200,300]), Pressure=1, Molecules=72, molecule='benzene', 
          Independent=0, potential='oplsaa', ignoreframes=200, includeframes=100000,
          simulation='gromacs', directory='', ensemble='NVT', spacing=1, hinge='DefaultHinge', phase='solid',
          Polymorphs=['p1', 'p2', 'p3'], refT=200, refdG=[0.000, 0.185, 0.306], refddG=[0.000, 0.019, 0.019],
          refdU=[0.000, 0.267, 0.240], absolutedU=[-5.624, -5.362, -5.380]):
#NSA: Is this needed?
    Colors = ['b', 'g', 'r', 'm', 'c', 'y', 'k', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g']

    if Independent == 0:
        Independent = Molecules

    # Hard set from old dictionary funciton
    refPot = 0
    ExtraPressures = []
    Pressures = np.ones(len(Temperatures), int)
    Pressures[len(Pressures) - len(ExtraPressures): len(Pressures)] = ExtraPressures
    Potentials = [potential]
    SimNAMES, Chargenames, PotNAMES = get_potential_info(potential)

    if (plot_out):
#        import matplotlib.cm as cm
#        from matplotlib.font_manager import FontProperties as FP
        font = {'family': 'normal',
                'weight': 'normal',
                'size': 14}
        matplotlib.rc('font', **font)
    
    # =============================================================================================
    # ENSURE THAT USER INPUTS ARE SENSIBLE
    # =============================================================================================
    # Pressure
    if Pressure < 0:
        print("Invalid Pressure: " + str(Pressure))
        sys.exit()
    
    # ENSEMBLE
    if ensemble != "NVE" and ensemble != "NVT" and ensemble != "NPT":
        print("Invalid Ensemble")
        print("Supported Ensembles: NVE NVT NPT")
        sys.exit()
    
    # =============================================================================================
    # FORMAT INPUTS
    # =============================================================================================
    # TEMPERATURE
    refk = -1
    for k, temp in enumerate(Temperatures):
        if temp == refT and refk == -1:
            refk = k + refPot * len(Temperatures)

    # =============================================================================================
    # READ IN RAW DATA
    # =============================================================================================
    # Constants.
    kB = 1.3806488e-23 * 6.0221413e23 / (1000.0 * 4.184)  # Boltzmann constant in kcal/mol/K
    
    # Parameters
    # How many states?
    K = len(Potentials) * len(Temperatures)
    
    #  maximum number of snapshots/simulation (could make this automated) - doesn't matter, as long as it's long enough.
    N_max = 3000
    
    # beta factor for the different temperatures
    beta_k = 1.0 / (kB * Temperatures)
    beta_k = np.tile(beta_k, (1, len(Potentials)))[0]
    
    # Conversion from kJ to kcal
    kJ_to_kcal = 0.2390057

    # This is the sampling efficiency for each potential in each combination of potentials
    Efficiency = np.zeros(K, float)
    
    # Allocate storage for simulation data
    # u_pklnT is the complete matrix of all energy data. 'p' is the polymorph, 'k' is the sampled state, 'l' is the
    #    evaluated state, 'n' is the sample number,  and T is the energy term
    #u_pklnT = np.zeros([len(Polymorphs), K, K, N_max, 20])
    
    # N_k[k] is the total number of snapshots from alchemical state k
    N_k = np.zeros(K, np.int32)
    
    # Terms_l is the list of all energy terms for the 'l' state
    Terms_l = []
    
    # dA[p,i,k] is the free energy between potential 0 and state k for spacing i in polymorph p
    dA = np.zeros([len(Polymorphs), spacing + 1, K], float)
    
    # ddA[p,i,k] is the uncertainty in the free energy between potential 0 and state k for spacing i in polymorph p
    ddA = np.zeros([len(Polymorphs), spacing + 1, K], float)
    
    # dG[p,i,t] is the free energy between polymorph 1 and polymorph p for spacing i and temperature t
    dG = np.zeros([len(Polymorphs), spacing + 1, len(Temperatures)])
    
    # ddG[p,i,t] is the uncertanity in the free energy between polymorph 1 and polymorph p for spacing i and temperature t
    ddG = np.zeros([len(Polymorphs), spacing + 1, len(Temperatures)])
    
    # dS[p,i,t] is the relative entropy between polymorph 1 and polymorph p for spacing i and temperature t
    dS = np.zeros([len(Polymorphs), spacing + 1, len(Temperatures)])
    
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
    for p, polymorph in enumerate(Polymorphs):
        # Cycle through all sampled potentials
        for i, potential_k in enumerate(Potentials):
            for t in range(len(Temperatures)):
                k = len(Temperatures) * i + t
                # Cycle through all evaluated potentials
                for j, potential_l in enumerate(Potentials):
                    l = len(Temperatures) * j
    
                    # set up the state l hinges for molecules, charges, and alchemical transformations
                    if SimNAMES[j] == "GRO":
                        if Chargenames[i] == "":
                            charge_hinge_l = ""
                        else:
                            charge_hinge_l = "_" + Chargenames[j]
                    elif SimNAMES[j] == "TIN":
                        charge_hinge_l = ""
                    dirpath = polymorph + '/temperature/' + str(t) + '/'
                    fname_energy = dirpath + potential_l + charge_hinge_l + '_energy.xvg'

                    all_energy = panedr.edr_to_df(dirpath + 'PROD.edr')
                    [start_production, _, _] = timeseries.detectEquilibration(np.array(all_energy['Total Energy']))

                    # Read in the energy file data for state k evaluated in state l using HARVIST
                    #if k == 0:
                    #    Terms_l.append(Harvist.CreateTerms(fname_energy))

                    # Setting the end point of the simulation
                    N = len(np.array(all_energy['Total Energy'])[start_production:])
                    N_k[k] = N

                    #u_pklnT[p, k, l, :N, :len(Terms_l[j])] = Harvist.GrabTerms(fname_energy, Terms_l[j],
                    #                                                           ignoreframes=start_production)[0]
                    u_kln[k, l, :N] = np.array(all_energy['Potential'])[start_production:]

                    # Now set these energies over all temperatures
                    #u_pklnT[p, k, l: (l + len(Temperatures)), :N, :len(Terms_l[j])] = u_pklnT[p, k, l, :N, :len(Terms_l[j])]
                    u_kln[k, l:(l + len(Temperatures)), :N] = u_kln[k, l, :N]

                    # Now read in the volumes and average them
                    V_pkn[p, t, :N] = np.array(all_energy['Volume'])[start_production:]
                    V_avg[p, t] = np.average(V_pkn[p, t, :N]) / float(Independent)
                    ddV_avg[p, t] = np.std(V_pkn[p, t, :N]) / N ** 0.5 / float(Independent)

                    # Now read in the lattice tensor and average them
                    if 'Box-XX' in list(all_energy):
                        box_letters = ['XX', 'YY', 'ZZ', 'YX', 'ZX', 'ZY']
                    else:
                        box_letters = ['X', 'Y', 'Z']

                    sign = np.sign(md.load(dirpath + 'pre_EQ.gro').unitcell_vectors[0].T)
                    for s in range(3):
                        for j in range(3):
                            if sign[s, j] == 0.:
                                # Correcting for the sign of the lattice parameters
                                sign[s, j] = 1.
 
                    for b in range(len(box_letters)):
                        C_pkn[p, t, :N, box_place[b, 0], box_place[b, 1]] = np.array(all_energy['Box-' + box_letters[b]])[start_production:] * \
                                sign[box_place[b, 0], box_place[b, 1]] * 10
                    C_avg = np.average(C_pkn[p, t, :N], axis=0)
                    dC = np.std(C_pkn[p, t, :N], axis=0)
                    h_avg[p, t] = crystal_matrix_to_lattice_parameters(C_avg) 
                    dh[p, t] = np.absolute(crystal_matrix_to_lattice_parameters(C_avg + dC) - h_avg[p, t])

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
        
            # extract self-consistent weights and uncertainties
            (df_i, ddf_i, theta_i) = mbar.getFreeEnergyDifferences()

            # extract entropy
            [_, _, Delta_u_ij, _, Delta_s_ij, dDelta_s_ij] = mbar.computeEntropyAndEnthalpy()
            print("Free Energies Optained...")
        
            # Store the dimensionless results in the dA container
            dA[p, i, :] = df_i[refk]
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
#NSA: see next comment
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
            ddA[p, i, :] = ddf_u[refk]

    # Check the overlap it necessary
    if len(Temperatures) == 2:
        print("Overlap:")
        print(O_pij)
        pdb.set_trace()
    
    # =============================================================================================
    # FINALIZE THE RELATIVE FREE ENERGY AND ENTROPY
    # =============================================================================================
    for i in range(spacing + 1):
        for t, T in enumerate(Temperatures):
            for p in range(len(Polymorphs)):
                #print('HERE!!!!!', dA[p, i, t], dA[0, i, t], (dA[p, i, t] - dA[0, i, t]), (beta_k[t] * float(Independent)), float(T),  float(refT), refdG[p])
                dG[p, i, t] = (dA[p, i, t] - dA[0, i, t]) / (beta_k[t] * float(Independent)) + float(T) / float(refT) * \
                                                                                               refdG[p]
                ddG[p, i, t] = ((ddA[p, i, t] ** 2 + ddA[0, i, t] ** 2) / (beta_k[t] * float(Independent)) ** 2 +
                                float(T) / float(refT) * float(refddG[p]) ** 2) ** 0.5
                if p == 0:
                    continue
                dS[p, i, t] = (dU[p, t] - dU[0, t] - dG[p, i, t]) / float(T)
                ddS[p, i, t] = (ddU[p, t] ** 2 + ddU[p, t] ** 2 + ddG[p, i, t] ** 2) ** 0.5 / float(T)
    
    print("Polymorph Free Energy:")
    for p in range(len(Polymorphs)):
        print("%8.3f %8.3f" % (dG[p, spacing, len(Temperatures) - 1], ddG[p, spacing, len(Temperatures) - 1]))
    
    # =============================================================================================
    # PLOT THE RELATIVE FREE ENERGY VS TEMPERATURE
    # =============================================================================================

#    Temperatures2 = np.array([0] + [j for j in Temperatures])
#    Pressures2 = np.array([1] + [j for j in Pressures])
#    dG2 = np.insert(dG, 0, np.transpose(np.tile(refdU, (1, 1))), axis=2)
#    ddG2 = np.insert(ddG, 0, np.zeros([len(Polymorphs), spacing + 1]), axis=2)
#    dS2 = np.insert(dS, 0, np.zeros([len(Polymorphs), spacing + 1]), axis=2)
#    ddS2 = np.insert(ddS, 0, np.zeros([len(Polymorphs), spacing + 1]), axis=2)
    PlotPress = 1  # Pressure to plot the dGvT curve at
    Temperatures_P = Temperatures[Pressures == PlotPress]

    if plot_out == True:
        f, a = plt.subplots(1, 1)
        xlabel = 'Temperature (K)'
        ylabel = 'Relative Free Energy (kcal/mol)'

        a.errorbar(Temperatures, dG[0, 0, :], yerr=np.zeros(len(dG[0, 0, :]), float), linestyle='--', marker='.',
                   linewidth=2, alpha=0.6, color='b', label=Polymorphs[0])
        for p in range(len(Polymorphs)):
            if p == 0:
                continue
            if len(Potentials) > 1:
                a.errorbar(Temperatures, dG[p, 0, :], yerr=ddG[p, 0, :], linestyle='--', marker='.', linewidth=2,
                           alpha=0.6, color=Colors[p], label=Polymorphs[p])
            a.errorbar(Temperatures_P, dG[p, 1, Pressures == PlotPress], yerr=ddG[p, 1, Pressures == PlotPress],
                       linestyle='-', marker='.', linewidth=2, alpha=0.6, color=Colors[p], label=Polymorphs[p])
    
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)

    np.save('T_' + molecule + '_' + potential, Temperatures_P)
    for p, Poly in enumerate(Polymorphs):
        np.save('dGvT_' + molecule + '_' + Poly + '_' + potential, dG[p, spacing, Pressures == PlotPress])
        np.save('ddGvT_' + molecule + '_' + Poly + '_' + potential, ddG[p, spacing, Pressures == PlotPress])
        if len(Potentials) > 1:
            np.save('dGvT_' + molecule + '_' + Poly + '_' + potential + '_indirect', dG[p, 0, :])
            np.save('ddGvT_' + molecule + '_' + Poly + '_' + potential + '_indirect', ddG[p, 0, :])
            if spacing > 1:
                np.save('dGvT_' + molecule + '_' + Poly + '_' + potential + '_convergence', dG[p, :, :])
                np.save('ddGvT_' + molecule + '_' + Poly + '_' + potential + '_convergence', ddG[p, :, :])
        np.save('dS_' + molecule + '_' + Poly + '_' + potential, dS[p, spacing, :])
        np.save('ddS_' + molecule + '_' + Poly + '_' + potential, ddS[p, spacing, :])
    
    # =============================================================================================
    # PLOT THE RELATIVE ENTROPY VS TEMPERATURE
    # =============================================================================================

    if plot_out == True:
        f, a = plt.subplots(1, 1)
        plt.xlabel('Temperature (K)', fontsize=20)
        plt.ylabel('Relative Entropy (kcal/molK)', fontsize=20)
        for p in range(len(Polymorphs)):
            a.errorbar(Temperatures2, dS[p, spacing, :], yerr=ddS[p, spacing, :], linestyle='--', marker='.', linewidth=2,
                       alpha=0.6, color=Colors[p], label=Polymorphs[p])
    
    # =============================================================================================
    # PLOT THE AVERAGE ENERGY VS TEMPERATURE
    # =============================================================================================

#    dU2 = np.insert(dU, 0, np.transpose(absolutedU), axis=1)
#    ddU2 = np.insert(ddU, 0, 0, axis=1)

    if plot_out == True:
        f, a = plt.subplots(1, 1)
        plt.xlabel('Temperature (K)', fontsize=20)
        plt.ylabel('Average Energy (kcal/mol)', fontsize=20)
        for p in range(len(Polymorphs)):
            a.errorbar(Temperatures2, dU[p, :], yerr=ddU[p, :], linestyle='--', marker='.', linewidth=2, alpha=0.6,
                       color=Colors[p], label=Polymorphs[p])
        plt.legend(loc='upper left')

    for p, Poly in enumerate(Polymorphs):
        np.save('UvT_' + molecule + '_' + Poly + '_' + potential, dU[p, :])
    
    # =============================================================================================
    # PLOT THE AVERAGE BOX VOLUME VS TEMPERATURE
    # =============================================================================================
    
    if plot_out == True:
        f, a = plt.subplots(1, 1)
        plt.xlabel('Temperature (K)', fontsize=20)
        for p in range(len(Polymorphs)):
            a.errorbar(Temperatures, V_avg[p, :], yerr=ddV_avg[p, :], linestyle='--', marker='.', linewidth=2,
                       alpha=0.6, color=Colors[p], label=Polymorphs[p])

    for p, Poly in enumerate(Polymorphs):
        np.save('VvT_' + molecule + '_' + Poly + '_' + potential, V_avg[p, :])

    # =============================================================================================
    # SAVE THE AVERAGE BOX VECTORS AND ANGLES VS TEMPERATURE
    # =============================================================================================

    for p, Poly in enumerate(Polymorphs):
        np.save('hvT_' + molecule + '_' + Poly + '_' + potential, h_avg[p, :])
        np.save('dhvT_' + molecule + '_' + Poly + '_' + potential, dh[p, :])

    # =============================================================================================
    # PLOT THE DIFFERENCE IN AVERAGE ENERGY VS TEMPERATURE
    # =============================================================================================
    
    if plot_out == True:
        f, a = plt.subplots(1, 1)
        plt.xlabel('Temperature (K)', fontsize=20)
        plt.ylabel('Average Energy Difference (kcal/mol)', fontsize=20)
        for p in range(len(Polymorphs)):
            if p == 0:
                yerror = np.zeros(len(ddU[0, :]))
            else:
                yerror = []
                for t in range(len(ddU[0, :])):
                    yerror.append(float((ddU[0, t] ** 2 + ddU[p, t] ** 2) ** 0.5))
            a.errorbar(Temperatures2, dU[p, :] - dU[0, :], yerr=yerror, linestyle='--', marker='.', linewidth=2,
                       alpha=0.6, color=Colors[p], label=Polymorphs[p])
    
    # Save the data for future use.
    for p, Poly in enumerate(Polymorphs):
        np.save('dUvT_' + molecule + '_' + Poly + '_' + potential, dU[p, :] - dU[0, :])
        np.save('ddUvT_' + molecule + '_' + Poly + '_' + potential, (ddU[p, :] ** 2 + ddU[0, :] ** 2) ** 0.5)

    if plot_out == True:
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    # =============================================================================================
    # READ IN USER INPUTS
    # =============================================================================================
    parser = OptionParser()
    parser.add_option('-p', '--plot', dest='plot', help='Plot output (default false)', default=True, action='store_true')
#    parser.add_option('-T', dest='temp', help='Temperature', default=200)
    parser.add_option('-P', dest='Pressure', help='Pressure', default=1)
#    parser.add_option('-n', dest='polymorphs', help='Polymorphs to analyze', default='p1 p2 p3')
    parser.add_option('-N', dest='molecules', help='number of supercell molecules', default=72)
    parser.add_option('-M', dest='molecule', help='name of the molecule', default='benzene')
    parser.add_option('-I', dest='independent', help='number of independent molecules', default='Same')
    parser.add_option('-u', dest='potential', help='potential to create the phase diagram in', default='oplsaa')
    parser.add_option('-i', dest='ignoreframes', help='Initial frames to ignore', default=200)
    parser.add_option('-j', dest='includeframes', help='Number of frames to include', default=100000)
    parser.add_option('-z', dest='simulation', help='The simulation package that was used', default='gromacs')
    parser.add_option('-d', dest='directory', help='Parent directory of the reweight directories', default='None')
    parser.add_option('-E', dest='ensemble', help='Simulation Ensemble', default='NVT')
    parser.add_option('-s', dest='spacing', help='spacing for the plot rolling back the state 0 sampling', default=1)
    parser.add_option('-H', '--hinge', dest='hinge', help='Optional string at end of jobs', default='DefaultHinge')
    
    (options, args) = parser.parse_args()
    plot_out = options.plot
    Temp = options.temp
    Pressure = int(options.Pressure)
    Molecules = int(options.molecules)
    molecule = options.molecule
    if options.independent == 'Same':
        Independent = Molecules
    else:
        Independent = int(options.independent)
    potential = options.potential
    ignoreframes = int(options.ignoreframes)
    includeframes = int(options.includeframes)
    simulation = options.simulation
    directory = options.directory
    if directory == 'None':
        directory = ''
    ensemble = options.ensemble
    spacing = int(options.spacing)
    hinge = options.hinge
    phase = "solid"
    
    Polymorphs, refT, refdG, refddG, refdU, absolutedU = old_systems_dictionary(potential, molecule)
    dGvsT(plot_out=plot_out, Temperatures=Temp, Pressure=Pressure, Molecules=Molecules, molecule=molecule,
          Independent=Independent, potential=potential, ignoreframes=ignoreframes, includeframes=includeframes, 
          simulation=simulation, directory=directory, ensemble=ensemble, spacing=spacing, hinge=hinge, phase=phase, 
          Polymorphs=Polymorphs, refT=refT, refdG=refdG, refddG=refddG, refdU=refdU, absolutedU=absolutedU)

