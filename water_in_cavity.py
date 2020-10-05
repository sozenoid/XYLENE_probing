#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 14:09:06 2020

@author: macenrola
"""

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
from subprocess import call
import os
from mol_ops_amber import make_pdb_with_named_residues
# sys.path.append('/home/macenrola/Thesis/Molecule_operations/')
working_directory = '/home/macenrola/Documents/water_in_cavity_CB8'
antepath = '/home/macenrola/anaconda3/envs/chemts/bin/antechamber'



if __name__ =="__main__":
	os.chdir(working_directory)