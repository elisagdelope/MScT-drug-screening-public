# -*- coding: cp1252 -*-
# ----------------------------------------------------------
# Copyright (C) 2018 PHARAMACELERA S.L.
# All rights reserved.
# 
# File: TFM_getH3mols.py
#
# Created on 14/06/2018
# ----------------------------------------------------------

##ELISA G. DE LOPE
##ATOM TYPES - MOPAC PROJECT


import mol2_sdf
import sys
def get_data(myfile):
    """Transform each molecule in the file into a molecule object with attributes (check mol2.py file)"""
    my_mol=mol2_sdf.read_molecule(myfile)
    return my_mol

myoutput_zero = open("H3_analysis.txt", "w")
atfile = open(sys.argv[1], 'r')
mopacfile = open(sys.argv[2], 'r')

my_mol_at = get_data(atfile)
my_mol_mopac = get_data(mopacfile)
while my_mol_at != None and my_mol_mopac != None:
    if "H3" in my_mol_at.lAtomTypes:
        myindex = my_mol_at.lAtomTypes.index("H3")
        myhelevalue = my_mol_mopac.lMopacValues[myindex][0]
        difsq = float((-0.7093-myhelevalue)**2)
        if difsq > 167.049 :
            s = my_mol_at.sName.rstrip("\n")
            myoutput_zero.write(s)
    my_mol_at = get_data(atfile)
    my_mol_mopac = get_data(mopacfile)

myoutput_zero.close()
