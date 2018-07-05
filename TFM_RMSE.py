# -*- coding: cp1252 -*-
# ----------------------------------------------------------
# Copyright (C) 2018 PHARAMACELERA S.L.
# All rights reserved.
# 
# File: TFM_RMSE.py
#
# Created on 08/06/2018
# ----------------------------------------------------------

##ELISA G. DE LOPE
##ATOM TYPES - MOPAC PROJECT

import os
import mol2_sdf
import sys
from math import sqrt

atdic = {
"Br":[-0.2726,-0.3234,1.3069],
"C1":[-0.0122,-0.0453,0.1518],
"C10":[-0.0175,-0.0404,0.1401],
"C11":[0.0080,-0.0330,0.1075],
"C12":[0.0100,-0.0284,0.1134],
"C13":[-0.0502,-0.0695,0.2990],
"C14":[-0.0354,-0.0889,0.3831],
"C15":[-0.0215,-0.0668,0.2874],
"C16":[-0.0240,-0.0640,0.2754],
"C17":[-0.0229,-0.0597,0.2568],
"C18":[-0.0413,-0.1009,0.3537],
"C19":[-0.0457,-0.0797,0.3428],
"C2":[0.0093,-0.0293,0.1035],
"C20":[-0.0448,-0.0800,0.3432],
"C21":[-0.0463,-0.0802,0.3440],
"C22":[-0.0360,-0.0786,0.3306],
"C23":[-0.0292,-0.0865,0.3573],
"C24":[-0.0502,-0.0736,0.3125],
"C25":[-0.0443,-0.0830,0.3557],
"C26":[-0.0447,-0.0928,0.3468],
"C3":[-0.0211,-0.0415,0.1416],
"C4":[0.0087,-0.0315,0.1047],
"C5":[-0.0467,-0.0862,0.3572],
"C6":[-0.0459,-0.0953,0.3476],
"C7":[-0.0690,-0.1625,0.6886],
"C8":[-0.0174,-0.0572,0.1911],
"C9":[-0.0557,-0.0568,0.1897],
"CS":[-0.0490,-0.0696,0.2974],
"Cl":[-0.1327,-0.3084,1.6361],
"F":[-0.2817,-0.2491,0.8556],
"H1":[-0.0658,-0.2059,0.4760],
"H2":[-0.7093,-0.1688,0.2321],
"H3":[-0.4996,-0.1147,0.1794],
"H4":[-1.1094,-0.1715,0.2369],
"HS":[-1.2935,-0.1769,0.2345],
"I":[-0.1375,-0.3402,1.3470],
"N1":[-0.8965,-0.0757,1.0249],
"N11":[-0.0442,-0.1489,0.3481],
"N12":[-0.1035,-0.0754,0.1520],
"N13":[-0.0946,-0.0849,0.1718],
"N14":[-3.4998,-0.1544,0.3058],
"N2":[-0.2892,-0.1420,0.7459],
"N3":[-0.9094,-0.0677,1.0151],
"N4":[-0.2883,-0.1403,0.7462],
"N5":[-1.4969,-0.2375,1.2650],
"N6":[-0.7002,-0.1667,0.3362],
"N7":[-0.0553,-0.0810,0.1634],
"N8":[-0.0683,-0.0786,0.1586],
"N9":[-1.5045,-0.2630,0.5295],
"NS":[-0.0967,-0.1266,0.1583],
"O1":[-0.2992,-0.1642,0.7102],
"O10":[-1.7004,-0.2423,1.0566],
"O11":[-2.4989,-0.2355,1.0449],
"O2":[-0.7069,-0.1164,0.7431],
"O3":[-0.2988,-0.1614,0.6999],
"O4":[-0.7008,-0.1598,0.6921],
"O5":[-0.9066,-0.2467,1.0694],
"O6":[-2.9009,-0.2301,1.0076],
"O7":[-2.0599,-0.2453,1.0744],
"O8":[-1.5026,-0.2435,1.0634],
"O9":[-2.0986,-0.2422,1.0598],
"OS":[-4.6982,-0.2284,0.9820],
"P":[-0.0946,-0.1079,0.8818],
"S1":[-0.1079,-0.2495,1.5100],
"S2":[-0.4981,-0.0995,0.5436],
"S3":[-0.0643,-0.2439,1.5347]}

def get_data(myfile):
    """Transform each molecule in the file into a molecule object with attributes (check mol2_sdf.py file)"""
    my_mol=mol2_sdf.read_molecule(myfile)
    return my_mol

def countats(myat,myats_dict):
    """Adds 1 for the atom-type counter when it is found"""
    if myat not in myats_dict:
        myats_dict[myat] = int(1)
    else:
        myats_dict[myat]+= 1
    return (myats_dict)

def comparevalues(atdic, valueslist_at, valueslist_mopac, myats_dict, errordic):
    """Generates an errordic dictionary with the squared error addition for each atom-type.
        Generates a myats_dict dictionary with the number of atoms for each atom-type.""" 
    for i in range(0, len(valueslist_at)):  
        if valueslist_at[i] in atdic.values(): #check which atom-type
            myat = atdic.keys()[atdic.values().index(valueslist_at[i])] #get atom-type
            myats_dict = countats(myat, myats_dict)
            if myat not in errordic.keys():              
                hele = (valueslist_mopac[i][0] - valueslist_at[i][0])**2
                hcav = (valueslist_mopac[i][1] - valueslist_at[i][1])**2
                hvwa = (valueslist_mopac[i][2] - valueslist_at[i][2])**2
                errordic[myat]=[hele,hcav,hvwa]
            else:
                errordic[myat][0] += (valueslist_mopac[i][0] - valueslist_at[i][0])**2
                errordic[myat][1] += (valueslist_mopac[i][1] - valueslist_at[i][1])**2
                errordic[myat][2] += (valueslist_mopac[i][2] - valueslist_at[i][2])**2
        else:
            print "Atom type not found!"
    return (errordic, myats_dict)



def mse(errordic, myats_dict, msestats):
    for element in errordic.keys():
        mse_ele = float(errordic[element][0]/myats_dict[element])
        mse_cav = float(errordic[element][1]/myats_dict[element])
        mse_vwa = float(errordic[element][2]/myats_dict[element])
        msestats[element] = [mse_ele,mse_cav,mse_vwa]
    return(msestats)

def rmse(msestats, rmsestats):
    for element in msestats.keys():
        rmse_ele = sqrt(msestats[element][0])
        rmse_cav = sqrt(msestats[element][1])
        rmse_vwa = sqrt(msestats[element][2])
        rmsestats[element] = [rmse_ele,rmse_cav,rmse_vwa]
    return(rmsestats)

#-------------------------------------------------------------------------------------------------------------
##MAIN

atfile = open("ace_pruebaats.mol2", 'r')
mopacfile = open("ace_pruebastandard.mol2", 'r')

my_mol_at = get_data(atfile)
my_mol_mopac = get_data(mopacfile)
errordic={}
myats_dict={}
while my_mol_at != None and my_mol_mopac != None:
    atslist_at = my_mol_at.lAtomTypes
    atslist_mopac =  my_mol_mopac.lAtomTypes
    valueslist_at = my_mol_at.lMopacValues
    valueslist_mopac = my_mol_mopac.lMopacValues
    (errordic, myats_dict) = comparevalues(atdic, valueslist_at, valueslist_mopac, myats_dict, errordic)
        
    my_mol_at = get_data(atfile)
    my_mol_mopac = get_data(mopacfile)
    
msestats ={}
rmsestats ={}
msestats = mse(errordic, myats_dict, msestats)
rmsestats = rmse(msestats, rmsestats)

MSE_output = open("MSE_error.txt", "a")
RMSE_output = open("RMSE_error.txt", "a")

print >> MSE_output, "AT", "MSE_hele", "MSE_hcav","MSE_hvwa"  
print >> RMSE_output, "AT", "RMSE_hele", "RMSE_hcav","RMSE_hvwa"  

for at in sorted(msestats.keys()):
    print at, msestats[at][0], msestats[at][1], msestats[at][2]
    print at, rmsestats[at][0], rmsestats[at][1], rmsestats[at][2]
    print >> MSE_output, at, msestats[at][0], msestats[at][1], msestats[at][2]  
    print >> RMSE_output, at, rmsestats[at][0], rmsestats[at][1], rmsestats[at][2] 

MSE_output.close()
RMSE_output.close()

