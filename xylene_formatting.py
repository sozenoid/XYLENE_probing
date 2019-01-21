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


def align_protonated_xylene_according_to_methyl(rdkitmol):
	"""
	:param rdkitmol: Aligns the protonated xylenes according to the two methyl groups in the molecule, the protonated hydrogen is placed in negative x domain
	:return: a rdkit mol with adapted geometries
	"""
	def get_atom_coords_for_projection(rdkitmol):
		"""
		:param rdkitmol: the very same rdkit mol
		:return: the atom coords that one would use for the projection are returned, that is the methyl carbon not involved in the reaction and the other sp2 carbons
		"""
		mol = rdkitmol
		Chem.SanitizeMol(mol)
		mol.UpdatePropertyCache(strict=False)
		pattpos = Chem.MolFromSmarts('[CH3][CH]')
		pos = mol.GetSubstructMatches(pattpos)[0][-1]
		mol.GetAtomWithIdx(pos).SetFormalCharge(1)

		pattalign = Chem.MolFromSmarts('[CH3]~[CH0]~[C;!+]~[C]')
		print mol.GetSubstructMatches(pattalign)
		for at in mol.GetAtoms():
			print at.GetExplicitValence(), at.GetImplicitValence(), at.GetFormalCharge(), at.GetIsAromatic(), at.GetIdx(), at.GetHybridization()
		return []
	get_atom_coords_for_projection(rdkitmol)
	mol = None
	return rdkitmol

if __name__ == "__main__":
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	mol = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/oXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	Chem.MolToMolFile(align_protonated_xylene_according_to_methyl(mol), '/home/macenrola/Desktop/lowl.sdf')