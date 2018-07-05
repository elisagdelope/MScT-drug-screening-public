# -*- coding: cp1252 -*-
# ----------------------------------------------------------
# Copyright (C) 2018 PHARAMACELERA S.L.
# All rights reserved.
# 
# File: TFM_Script1.py
#
# Created on 17/04/2018
# ----------------------------------------------------------

##ELISA G. DE LOPE
##ATOM TYPES - MOPAC PROJECT


import os
import mol2_sdf
import sys

def get_data(myfile):
    """Transform each molecule in the file into a molecule object with attributes (check mol2.py file)"""
    my_mol=mol2_sdf.read_molecule(myfile)
    return my_mol

def contributions_mopac(at_list,mopac_list, countNA, countats):
    """Collects data in four dictionaries, one per contribution type, so that {C1:[Pele1,Pele2,Pele3,Pele4,...],C2:[Pele1,Pele2,Pele3,...]},
        {C1:[Pcav1,Pcav2,Pcav3,Pcav4,...],C2:[Pcav1,Pcav2,Pcav3,...]}, same for Pvwa} and another one for total logP values"""
    
    for i in range(0,len(at_list)):
        countats = countats + 1
        if at_list[i] != "N/A":
            mopac_values = mopac_list[i]
            logP=sum(mopac_values)
            filename_logP = "temp/logP_"+at_list[i]+".txt"
            filename_logPele = "temp/logPele_"+at_list[i]+".txt"
            filename_logPcav = "temp/logPcav_"+at_list[i]+".txt"
            filename_logPvwa = "temp/logPvwa_"+at_list[i]+".txt"
            if os.path.isfile(filename_logP):
                myfile_logP=open (filename_logP, "a")
                myfile_logPele=open (filename_logPele, "a")
                myfile_logPcav=open (filename_logPcav, "a")
                myfile_logPvwa=open (filename_logPvwa, "a")
                print >> myfile_logP, logP
                print >> myfile_logPele, mopac_values[0]
                print >> myfile_logPcav, mopac_values[1]
                print >> myfile_logPvwa, mopac_values[2]                 
            else:
                myfile_logP=open (filename_logP, "w")
                myfile_logPele=open (filename_logPele, "w")
                myfile_logPcav=open (filename_logPcav, "w")
                myfile_logPvwa=open (filename_logPvwa, "w")
                print >> myfile_logP, logP
                print >> myfile_logPele, mopac_values[0]
                print >> myfile_logPcav, mopac_values[1]
                print >> myfile_logPvwa, mopac_values[2]            
            myfile_logP.close()
            myfile_logPele.close()
            myfile_logPcav.close()
            myfile_logPvwa.close()
        else:
            countNA = countNA + 1
    return (countNA,countats)
        
def build_atlist (at_list, complete_list):
    at_list = list(set(at_list))
    for atom in at_list:
        if atom not in complete_list:
            complete_list.append(atom)
    return(complete_list)


def check_outlier(logP, mopac_values, c_type):
    """Checks if value is an outlier"""
    if c_type == "logP" and (logP <= -10 or logP >= 10):
        out = True
    elif c_type == "logPele" and (mopac_values[0] <= -8 or mopac_values[0] >= 2):
        out = True
    elif c_type == "logPcav" and (mopac_values[1] <= -1 or mopac_values[1] >= 1):
        out = True
    elif c_type == "logPvwa" and (mopac_values[2] <= -1 or mopac_values[2] >= 3):
        out = True
    else:
        out = False
    return (out)

def outliers_molcounter(dict_outliers, c_type, mol_name, at_list, mopac_list, mol_outlist):
    """Collects data from outliers"""
    mol_outliers = False
    at_uniques = list(set(at_list))
    for at in at_uniques:
        if at not in dict_outliers:
            dict_outliers[at] = int(0)  
    for i in range(0,len(at_list)):
        if not at_list[i] in at_uniques:
            pass
        else:            
            mopac_values = mopac_list[i]
            logP=sum(mopac_values)
            (at_outlier) = check_outlier (logP, mopac_values, c_type)
            if at_outlier == True:
                dict_outliers[at_list[i]] += 1
                at_uniques.remove(at_list[i])
                mol_outliers = True
    if mol_outliers == True:
        mol_outlist.append(mol_name.rstrip())
    return (dict_outliers, mol_outlist)


def allouts_count(dictionary, c_type, at_list, mopac_list):
    """Counts outliers"""
    for i in range(0,len(at_list)):
        if at_list[i] not in dictionary:
            dictionary[at_list[i]] = int(0)        
        mopac_values = mopac_list[i]
        logP=sum(mopac_values)
        (at_outlier) = check_outlier (logP, mopac_values, c_type)
        if at_outlier == True:
            dictionary[at_list[i]] += 1
    return(dictionary)


def print_outliers (atlist, dictlogP, dictlogPele, dictlogPcav, dictlogPvwa, fileout):
    """Prints outliers matrix"""
    print >> fileout, "\tlogP\tlogPele\tlogPcav\tlogPvwa"
    for element in sorted(my_atlist):
        vals = [0,0,0,0]
        if element in dictlogP.keys():
            vals[0] = dictlogP[element]
        if element in dictlogPele.keys():
            vals[1] = dictlogPele[element]
        if element in dictlogPcav.keys():
            vals[2] = dictlogPcav[element]
        if element in dictlogPvwa.keys():
            vals[3] = dictlogPvwa[element]
        vals = [str(v) for v in vals]
        print >> fileout, element, "\t", "\t".join(vals)


        
#----------------MAIN-------------------------------------------
#Open file with atom types and file with contributions

atfile = open(sys.argv[1], 'r')
mopacfile = open(sys.argv[2], 'r')


#---------------------------------------------------------------
#Get the data from each molecule and atom type and store it in files
my_mol_at = get_data(atfile)
my_mol_mopac = get_data(mopacfile)

#store data in files
if not os.path.exists("temp"):
    os.makedirs("temp")

#store outliers data
countNA = 0
countats = 0
dict_outlogP_mol = {}
dict_outlogPele_mol = {}
dict_outlogPcav_mol = {}
dict_outlogPvwa_mol = {}
dict_outlogP = {}
dict_outlogPele = {}
dict_outlogPcav = {}
dict_outlogPvwa = {}
mol_outlist_logP = []
mol_outlist_logPele = []
mol_outlist_logPcav = []
mol_outlist_logPvwa = []
my_atlist = []
num_mols = 0
while my_mol_at != None and my_mol_mopac != None:
    (countNA,countats) = contributions_mopac(my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues, countNA, countats) #get and store data to files
    #keep outliers info per molecule
    (dict_outlogP_mol, mol_outlist_logP)=outliers_molcounter(dict_outlogP_mol, "logP", my_mol_at.sName, my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues, mol_outlist_logP)
    (dict_outlogPele_mol, mol_outlist_logPele)=outliers_molcounter(dict_outlogPele_mol, "logPele", my_mol_at.sName, my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues, mol_outlist_logPele)
    (dict_outlogPcav_mol, mol_outlist_logPcav)=outliers_molcounter(dict_outlogPcav_mol, "logPcav", my_mol_at.sName, my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues, mol_outlist_logPcav)
    (dict_outlogPvwa_mol, mol_outlist_logPvwa)=outliers_molcounter(dict_outlogPvwa_mol, "logPvwa", my_mol_at.sName, my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues, mol_outlist_logPvwa)
    #keep outliers total info
    (dict_outlogP)=allouts_count(dict_outlogP, "logP", my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues)
    (dict_outlogPele)=allouts_count(dict_outlogPele, "logPele", my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues)
    (dict_outlogPcav)=allouts_count(dict_outlogPcav, "logPcav", my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues)
    (dict_outlogPvwa)=allouts_count(dict_outlogPvwa, "logPvwa", my_mol_at.lAtomTypes, my_mol_mopac.lMopacValues)
    #build atom types list
    my_atlist = build_atlist(my_mol_at.lAtomTypes, my_atlist)
    num_mols += 1
    my_mol_at = get_data(atfile)
    my_mol_mopac = get_data(mopacfile)


#---------------------------------------------------------------
#Print N/A info

myfile_NA = open ("NA_ats.txt", "w")
print >> myfile_NA, "Number of N/A atoms:", str(countNA)
print >> myfile_NA, "Number of total atoms", str(countats)
print >> myfile_NA, "% of N/As:", str(float(countNA/countats)*100)
myfile_NA.close()

#---------------------------------------------------------------
#Print outliers info to matrix-like file
myfile_outliers_atmatrix=open ("at_outliers_matrix.txt", "w")
myfile_outliers_molmatrix=open ("mol_outliers_matrix.txt", "w")    

print_outliers(my_atlist, dict_outlogP_mol, dict_outlogPele_mol, dict_outlogPcav_mol, dict_outlogPvwa_mol, myfile_outliers_molmatrix)
print_outliers(my_atlist, dict_outlogP, dict_outlogPele, dict_outlogPcav, dict_outlogPvwa, myfile_outliers_atmatrix)

myfile_outliers_atmatrix.close()
myfile_outliers_molmatrix.close()


#---------------------------------------------------------------
#Print outliers info to mol summary file
       
myfile_outliers_molsummary=open ("mol_outliers_summary.txt", "w")
print >> myfile_outliers_molsummary, "Number of molecules containing at least one atom type with outlier logP:", len(mol_outlist_logP)
print >> myfile_outliers_molsummary, "Number of molecules containing at least one atom type with outlier logPele:", len(mol_outlist_logPele)
print >> myfile_outliers_molsummary, "Number of molecules containing at least one atom type with outlier logPcav:", len(mol_outlist_logPcav)
print >> myfile_outliers_molsummary, "Number of molecules containing at least one atom type with outlier logPvwa:", len(mol_outlist_logPvwa)
print >> myfile_outliers_molsummary, "Total number of analyzed molecules:", num_mols, "\n"
print >> myfile_outliers_molsummary, "LogP outlier-containing molecules are:", "\t".join(mol_outlist_logP)
print >> myfile_outliers_molsummary, "LogPele outlier-containing molecules are:", "\t".join(mol_outlist_logPele)
print >> myfile_outliers_molsummary, "LogPcav outlier-containing molecules are:", "\t".join(mol_outlist_logPcav)
print >> myfile_outliers_molsummary, "LogPvwa outlier-containing molecules are:", "\t".join(mol_outlist_logPvwa)
myfile_outliers_molsummary.close()
