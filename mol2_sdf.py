# -*- coding: cp1252 -*-
# ----------------------------------------------------------
# Copyright (C) 2015 PHARAMACELERA S.L.
# All rights reserved.
# 
# File: mol2_sdf.py
#
# Created on 10/04/2018
# ----------------------------------------------------------

##ELISA G. DE LOPE
##ATOM TYPES - MOPAC PROJECT

import re
import math
# ----------------------------------------------------------
# Wrapper of a MOL2 molecule.
# Contains:
#  - number of atoms
#  - number of bonds
#  - list of strings (lines from the MOL2 file that belongs
#    to the molecule
#  - list of atom types
#  - list of mopac contributions
# ----------------------------------------------------------
class Mol2Object(object):
	def __init__(self, s_name, i_natoms, i_nbonds, l_atoms_data, l_at_list, l_mopac_list, l_lines):
		self.sName = s_name
		self.iNumAtoms = i_natoms
		self.iNumBonds = i_nbonds
		self.lAtomsData = l_atoms_data
		self.lAtomTypes = l_at_list
		self.lMopacValues = l_mopac_list
		self.lLines = l_lines

# ----------------------------------------------------------
# Reads the next molecule from MOL2 file and returns the
# associated Mol2Object instance
# ----------------------------------------------------------
def read_molecule(f):
	line_list = []
	atoms_data = []
	at_list = []
	mopac_list = []
	num = 0

	# Molecule keyword
	while True:
		s_line = f.readline()
		if not s_line:
			return None
			break
		if (len(re.findall('MOLECULE', s_line))):
			break
	num = 0
	line_list.append(s_line)

	# Name of molecule
	s_line = f.readline()
	s_name = s_line
	line_list.append(s_line)
	num = num + 1

	# Num atoms and bonds
	s_line = f.readline()
	line_list.append(s_line)
	num = num + 1
	num_atoms = int(s_line.split()[0])
	num_bonds = int(s_line.split()[1])

	#Atoms
	while True:
		s_line = f.readline()
		if not s_line:
			return None
			break
		line_list.append(s_line)
		num = num + 1
		if (len(re.findall('ATOM', s_line))):
			break
	i = 0
	i_first_atom_line = num + 1
	while (i < num_atoms):
		s_line = f.readline()
		line_list.append(s_line)
                split_line = s_line.split()
                atoms_data.append(split_line) #THE REAL INTEREST
		
                if len(split_line) == 11: #mol2 with atom types
                        atomtype = split_line[9]
                        at_list.append(atomtype)
                elif len(split_line) == 13: #mol2 with contributions
                        mopac_values = [float(k) for k in split_line[10:13]]
		        mopac_list.append(mopac_values)
                else:
                        print "not contributions nor atomtypes file"
                        return None
                        break
                
                
		
		num = num + 1		
		i = i + 1
	s_line = f.readline()
	line_list.append(s_line)
	num = num + 1

	# Atom keyword
	i = 0
	i_first_bond_line = num + 1
	while (i < num_bonds):
		s_line = f.readline()
		if not s_line:
			return None
			break
		line_list.append(s_line)
		num = num + 1
		i = i + 1

	mol2 = Mol2Object(s_name, num_atoms, num_bonds, atoms_data, at_list, mopac_list, line_list)
	return mol2

# ----------------------------------------------------------
# Writes a MOL2 molecule to a file
# ----------------------------------------------------------
def write_molecule(f, obj_mol2):
	for line in obj_mol2.lLines:
		f.write(line)

# ----------------------------------------------------------
# Writes a MOL2 molecule to a file
# ----------------------------------------------------------
def toString(obj_mol2):
	out = ""
	for line in obj_mol2.lLines:
		out = out + line
	return out

