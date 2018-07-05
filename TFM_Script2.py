# -*- coding: cp1252 -*-
# ----------------------------------------------------------
# Copyright (C) 2018 PHARAMACELERA S.L.
# All rights reserved.
# 
# File: TFM_Script2.py
#
# Created on 05/05/2018
# ----------------------------------------------------------

##ELISA G. DE LOPE
##ATOM TYPES - MOPAC PROJECT


import os
import matplotlib.pyplot as plt
import mol2
import sys
import numpy as np

def plotting_contributions(filename, datafile_list, mode_list, c_type, mydir, outliers):
    """Generates plots for LogP contributions distributions of values"""
    if outliers == True:
        if c_type == "logP_":
            datafile_list = [k for k in datafile_list if k < 10 and k > -10]
        elif c_type == "logPele_":
            datafile_list = [k for k in datafile_list if k < 2 and k > -8]            
        elif c_type == "logPcav_":
            datafile_list = [k for k in datafile_list if k < 1 and k > -1]
        elif c_type == "logPvwa_":
            datafile_list = [k for k in datafile_list if k < 3 and k > -1]
    else:
        pass

    n, b, patches = plt.hist(x=datafile_list, bins = 100, range=(-10,10), histtype='bar', align='mid', orientation='vertical',facecolor='#FC6600')    
    plt.xlim(-10,10)
    plt.xlabel('Value')
    plt.ylabel('Count')
    plt.grid(True)

    #construct mode interval and compute the mean:
    bin_max = list(n).index(n.max())
    mode_int = [m for m in datafile_list if m > b[bin_max] and m < b[bin_max + 1]]
    mymode = np.mean(mode_int)
        
    if c_type.startswith("logP") and c_type != "logP":
        plt.title(filename[8:-4]+'_'+c_type)
        plt.savefig(os.path.join(mydir, filename[8:-4]))
        s = filename[8:-4]+"\t"+str(mymode)
    else:
        plt.title(filename[5:-4]+'_'+c_type)
        plt.savefig(os.path.join(mydir, filename[5:-4]))
        s = filename[5:-4]+"\t"+str(mymode)

    plt.close()  
    mode_list.append(s)
    
    return (mode_list)
    
#----------------MAIN-------------------------------------------
#Plot data
outliers = eval(sys.argv[1])
    
mydir_logP = os.path.join(os.getcwd(),"logP_plots")
mydir_logPele = os.path.join(os.getcwd(),"logPele_plots")
mydir_logPcav = os.path.join(os.getcwd(),"logPcav_plots")
mydir_logPvwa = os.path.join(os.getcwd(),"logPvwa_plots")
mydir_modes = os.path.join(os.getcwd(),"mode_values")
if not os.path.exists(mydir_logP):
    os.makedirs(mydir_logP)
    os.makedirs(mydir_logPele)
    os.makedirs(mydir_logPcav)
    os.makedirs(mydir_logPvwa)
    os.makedirs(mydir_modes)

mode_list_logP = []
mode_list_logPele = []
mode_list_logPcav = []
mode_list_logPvwa = []

for filename in sorted(os.listdir(os.path.join(os.getcwd(),"temp"))):
    print filename
    datafile_list = []
    with open("temp/"+filename) as infile:
        for line in infile:
            datafile_list.append(float(line))
    if filename.endswith(".txt") and filename.startswith("logP_"):
        mode_list_logP = plotting_contributions(filename, datafile_list, mode_list_logP, 'logP', mydir_logP, outliers)
    elif filename.endswith(".txt") and filename.startswith("logPele_"):
        mode_list_logPele = plotting_contributions(filename, datafile_list, mode_list_logPele, 'logPele', mydir_logPele, outliers)
    elif filename.endswith(".txt") and filename.startswith("logPcav_"):
        mode_list_logPcav = plotting_contributions(filename, datafile_list, mode_list_logPcav, 'logPcav', mydir_logPcav, outliers)
    elif filename.endswith(".txt") and filename.startswith("logPvwa_"):
        mode_list_logPvwa = plotting_contributions(filename, datafile_list, mode_list_logPvwa, 'logPvwa', mydir_logPvwa, outliers)        
#    os.remove("temp/"+filename)
#os.rmdir("temp")


#Print mode data
modefile_logP = open("mode_values/mode_logP.txt", "w")
modefile_logPele = open("mode_values/mode_logPele.txt", "w")
modefile_logPcav = open("mode_values/mode_logPcav.txt", "w")
modefile_logPvwa = open("mode_values/mode_logPvwa.txt", "w")

for i in range(0, len(mode_list_logP)):
    modefile_logP.write(mode_list_logP[i]+"\n")
    modefile_logPele.write(mode_list_logPele[i]+"\n")
    modefile_logPcav.write(mode_list_logPcav[i]+"\n")
    modefile_logPvwa.write(mode_list_logPvwa[i]+"\n")
modefile_logP.close()
modefile_logPele.close()
modefile_logPcav.close()
modefile_logPvwa.close()

