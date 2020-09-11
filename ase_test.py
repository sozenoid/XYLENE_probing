#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 00:45:20 2020

@author: macenrola
"""


from ase import Atoms, io, calculators
import ase.calculators.cp2k
from gpaw import GPAW, PW
from ase.optimize import BFGS
from ase.vibrations import Vibrations

mol = io.read("/home/macenrola/Documents/cb8_MV_electrochemistry_JIA/apbs_ref/benzene.pdb")
mol.center(vacuum=3.5)
print(mol.cell)
print(mol.positions)
c = ase.calculators.cp2k.CP2K()
CP2K.command="env OMP_NUM_THREADS=2 mpirun -np 4 cp2k"
print(mol.set_calculator(c))

BFGS(mol).run(fmax=0.01)
vib = Vibrations(mol)

vib.run()

print(vib.summary())
