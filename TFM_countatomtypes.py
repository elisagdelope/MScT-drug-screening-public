# -*- coding: cp1252 -*-
# ----------------------------------------------------------
# Copyright (C) 2018 PHARAMACELERA S.L.
# All rights reserved.
# 
# File: TFM_countatomtypes.py
#
# Created on 12/05/2018
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

def contributions_mopac(at_list,countNA, countats,myats_dict):
    """Collects data in four dictionaries, one per contribution type, so that {C1:[Pele1,Pele2,Pele3,Pele4,...],C2:[Pele1,Pele2,Pele3,...]}, {C1:[Pcav1,Pcav2,Pcav3,Pcav4,...],C2:[Pcav1,Pcav2,Pcav3,...]}, same for Pvwa} and another one for total logP values"""
    
    for i in range(0,len(at_list)):
        countats = countats + 1
        if at_list[i] != "N/A":
            if at_list[i] not in myats_dict:
                myats_dict[at_list[i]] = int(1)
            else:
                myats_dict[at_list[i]]+= 1
        else:
            countNA = countNA + 1
    return (countNA,countats, myats_dict)


atfile = open(sys.argv[1], 'r')


myats_dict= {}
my_mol_at = get_data(atfile)
countNA = 0
countats = 0
while my_mol_at != None:
    (countNA,countats, myats_dict) = contributions_mopac(my_mol_at.lAtomTypes, countNA, countats, myats_dict) #get and store data to files
    my_mol_at = get_data(atfile)
print "# NA:", countNA
print "# Total ats:", countats
for element in sorted(myats_dict.keys()):
    print "#", element, ":", myats_dict[element]
