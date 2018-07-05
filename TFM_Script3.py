# -*- coding: cp1252 -*-
# ----------------------------------------------------------
# Copyright (C) 2018 PHARAMACELERA S.L.
# All rights reserved.
# 
# File: TFM_Script3.py
#
# Created on 28/04/2018
# ----------------------------------------------------------

##ELISA G. DE LOPE
##ATOM TYPES - MOPAC PROJECT

import os
import mol2
import sys
from math import sqrt


def outliers_counter(filename, c_type, outliers):
    """Collects data from outliers"""
    counter = 0
    with open(filename) as infile:
        for line in infile:
            if c_type == "logP_" and (float(line) <= -10 or float(line) >= 10):
                counter = counter + 1
            elif c_type == "logPele_" and (float(line) <= -8 or float(line) >= 2):
                counter = counter + 1            
            elif c_type == "logPcav_" and (float(line) <= -1 or float(line) >= 1):
                counter = counter + 1
            elif c_type == "logPvwa_" and (float(line) <= -1 or float(line) >= 3):
                counter = counter + 1
    return (counter)

def outliers_check (val, c_type):
    """Checks for outliers"""
    if c_type == "logP_" and (val <= -10 or val >= 10):
        out = True
    elif c_type == "logPele_" and (val <= -8 or val >= 2):
        out = True            
    elif c_type == "logPcav_" and (val <= -1 or val >= 1):
        out = True
    elif c_type == "logPvwa_" and (val <= -1 or val >= 3):
        out = True
    else:
        out = False
    return (out)
    
def stdev (mynumbers, mymean, n_lines):
    """Computes standard deviation"""
    mysd = 0
    for element in mynumbers:
        mysd += (element - mymean)**2
    if n_lines == 1:
        mysd = 0.0
    else:
        mysd = sqrt(mysd / float(n_lines-1))
    return (mysd)
        
def stats (filename, dict_stats, c_type, outliers):
    """Compute stats"""
    mynumbers = []
    mysum = 0
    n_lines = 0
    total = 0
    mymax = None
    mymin = None
    with open(filename) as infile:
        for line in infile:
            total = total + 1
            if outliers == True and outliers_check(float(line), c_type) == True:
                linecount = False
            elif float(line) < -10 or float(line) > 10:
                linecount = False
            else:
                linecount = True
            if linecount == True: #if line not in range it will not be taken into account for statistics.
                mynumbers.append(float(line))
                mysum = mysum + float(line)
                n_lines = n_lines + 1            
                if float(line) > mymax or mymax == None:
                    mymax = float (line)
                if float(line) < mymin or mymin == None:
                    mymin = float (line)
    if n_lines == 0:
        mymean = None
        mysd = None
    else:
        mymean = round(mysum/n_lines,4)
        mysd = round(stdev(mynumbers, mymean, n_lines),4)
    outliers_count = outliers_counter(filename, c_type, outliers)

    if c_type == "logP_":
        dict_stats[filename[5:-4]] = [str(mymean), str(mysd), str(mymax), str(mymin), str(outliers_count), str(total), str(round(float(outliers_count)/float(total),6))] 
    else:
        dict_stats[filename[8:-4]] = [str(mymean), str(mysd), str(mymax), str(mymin), str(outliers_count), str(total), str(round(float(outliers_count)/float(total),6))] 
    return (dict_stats)


def statstofile (filename, dictstats):
    """Prints stats data to file"""
    fileout = open("stats/"+filename, "w")
    fileout.write("AtomType\tMean\tSDev\tMax\tMin\t#Outliers\t#Total\tOutfreq\n")
    s = "\t"
    for key in sorted(dictstats.keys()):
	print key
        print >> fileout, key[5:], "\t", s.join(dictstats[key])
    fileout.close()


 
#----------------MAIN-------------------------------------------
#Compute statistics 
outliers = eval(sys.argv[1])
    
TlogP_stats = {}
logPele_stats = {}
logPcav_stats = {}
logPvwa_stats = {}
for filename in sorted(os.listdir(os.path.join(os.getcwd(),"temp"))):
    print filename
    if filename.endswith(".txt") and filename.startswith("logP_"):
        (TlogP_stats)=stats("temp/"+filename, TlogP_stats, "logP_", outliers)
    elif filename.endswith(".txt") and filename.startswith("logPele_"):
        (logPele_stats)=stats("temp/"+filename, logPele_stats, "logPele_", outliers)
    elif filename.endswith(".txt") and filename.startswith("logPcav_"):
        (logPcav_stats)=stats("temp/"+filename, logPcav_stats, "logPcav_", outliers)
    elif filename.endswith(".txt") and filename.startswith("logPvwa_"):
        (logPvwa_stats)=stats("temp/"+filename, logPvwa_stats, "logPvwa_", outliers)
#    os.remove("temp/"+filename)
#os.rmdir("temp")



#----------------------------------------------------------------
#Write statistics to file
if not os.path.exists("stats"): 
    os.makedirs("stats")

statstofile("logP_statistics.txt", TlogP_stats)
statstofile("logPele_statistics.txt", logPele_stats)
statstofile("logPcav_statistics.txt", logPcav_stats)
statstofile("logPvwa_statistics.txt", logPvwa_stats)





