def format_mol_for_vasp(atm_list, atm_coords, flags_T_F, name):
	"""
	:param rdkitmol: Takes a rdkitmol
	:return: prints a POSCAR file for use in VASP operations, box size will be 3 Angstrom larger than the max in each directions
	"""
	def get_box_size(atm_coords):
		"""
		:param atm_coords: take in the cartesian coordinates of atoms and
		:return: returns the size of the naive bounding box in x,y,z
		"""
		x_width, y_width, z_width = 0 , 0 , 0
		return x_width, y_width, z_width


def align_xylenes(mol1, mol2, mol3, core='[CH3]~[CH0]~[C;!+]~[C;!+]~[C]~[CH]~[CH]'):
	"""
	:param mol1: o xylene
	:param mol2: m xylene
	:param mol3: p xylene
	:param core: the common fragment upon which align
	:return: None but will write the aligned versions to mol files
	"""
	for m in [mol1, mol2, mol3]:
		Chem.SanitizeMol(m)
		m.UpdatePropertyCache(strict=False)
		pattpos = Chem.MolFromSmarts('[CH3][CH]')
		pos = m.GetSubstructMatches(pattpos)[0][-1]
		m.GetAtomWithIdx(pos).SetFormalCharge(1)
		for at in m.GetAtoms():
			at.SetNoImplicit(True)
	core = Chem.MolFromSmarts(core, True)
	match1 = mol1.GetSubstructMatch(core)
	match2 = mol2.GetSubstructMatch(core)
	match3 = mol3.GetSubstructMatch(core)

	print match1, match2, match3

	AllChem.AlignMol(mol2, mol1,atomMap=zip(match2, match1), maxIters=1000, reflect=True)  # <- m2 is aligned to m1, return value is the RMSD for the alignment
	AllChem.AlignMol(mol3, mol1,atomMap=zip(match3, match1), maxIters=1000, reflect=False)  # <- m3 is aligned to m1, return value is the RMSD for the alignment

	Chem.MolToMolFile(mol1, '/home/macenrola/Desktop/m1.sdf')
	Chem.MolToMolFile(mol2, '/home/macenrola/Desktop/m2.sdf')
	Chem.MolToMolFile(mol3, '/home/macenrola/Desktop/m3.sdf')

	return


if __name__ == "__main__":
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import numpy as np

	oxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/oXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	mxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/mXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol',removeHs=False)
	pxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/pXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol',	removeHs=False)

	align_xylenes(oxylene, mxylene, pxylene)
# Chem.MolToMolFile(align_protonated_xylene_according_to_methyl(mol), '/home/macenrola/Desktop/pxylene.sdf')