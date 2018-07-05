# -*- coding: cp1252 -*-
# ----------------------------------------------------------
# Copyright (C) 2018 PHARAMACELERA S.L.
# All rights reserved.
# 
# File: TFM_extract_molsbyat.py
#
# Created on 20/05/2018
# ----------------------------------------------------------

##ELISA G. DE LOPE
##ATOM TYPES - MOPAC PROJECT

#Script to extract the names of the molecules containing a certain atom type.
#USAGE: python extract_molsbyat.py file_with_ats atomtype


import mol2
import sys

def get_data(myfile):
    """Transform each molecule in the file into a molecule object with attributes (check mol2.py file)"""
    my_mol=mol2.read_molecule(myfile)
    return my_mol

###--------MAIN----------------------------------------------------------
atfile = open(sys.argv[1], 'r')
target_at = sys.argv[2]
target_at = target_at.upper()

myoutput_file = open("mols_"+target_at+".txt", "w")

my_mol_at = get_data(atfile)
while my_mol_at != None:
    my_mol_at.lAtomTypes = [element.upper() for element in my_mol_at.lAtomTypes] 
    if target_at in my_mol_at.lAtomTypes:
        print target_at
        s = my_mol_at.sName
        myoutput_file.write(s)
    my_mol_at = get_data(atfile)

myoutput_file.close()
