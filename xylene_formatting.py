def format_mol_for_vasp(rdkitmol, name):
	"""
	:param rdkitmol: Takes a rdkitmol
	:return: prints a POSCAR file for use in VASP operations, box size will be 3 Angstrom larger than the max in each directions
	"""
	def remap_molecules(rdkitmol, isomertype):
		"""
		:param carbon_list: a list of carbon coordinates
		:param hydrogen_list: a list of carbon coordinates
		:param isomertype: p is the reference or 0, m is 1, o is 2
		:return: the same rdkitmol but with the atom indexing that is matched with the para protonated version of xylene
		"""
		table_pmo_carbon = [ # atom indexes for PARA, META, ORTHO
			[1, 7, 1], # CARBON
			[2, 4, 4], # CARBON
			[3, 5, 5], # CARBON
			[4, 6, 6], # CARBON
			[5, 2, 7], # CARBON
			[6, 3, 2], # CARBON
			[7, 13, 3],# CARBON
			[8, 1, 8], # CARBON
			[9, 16, 9],   # HYDROGEN
			[10, 17, 10], # HYDROGEN
			[11, 18, 11], # HYDROGEN
			[12, 11, 14], # HYDROGEN
			[13, 12, 15], # HYDROGEN
			[14, 14, 19], # HYDROGEN
			[15, 19, 12], # HYDROGEN
			[16, 8, 18],  # HYDROGEN
			[17, 9, 17],  # HYDROGEN
			[18, 10, 16], # HYDROGEN
			[19, 15, 13]  # HYDROGEN
		]

		ordered_carbon = []
		ordered_hydrogen = []
		conf = rdkitmol.GetConformer(-1)
		for i, c in enumerate([x[isomertype] for x in table_pmo_carbon]):
			pos = list(conf.GetAtomPosition(c-1))
			if i<8:
				ordered_carbon.append(pos)
			else:
				ordered_hydrogen.append(pos)

		return ordered_carbon, ordered_hydrogen

	def get_box_size(atom_pos_dic):
		"""
		:param atm_coords: take in the cartesian coordinates of atoms and
		:return: returns the size of the naive bounding box in x,y,z
		"""
		all_coords = []
		for k in atom_pos_dic:
			for els in atom_pos_dic[k]:
				all_coords.append(els)
		x_width, y_width, z_width = max(all_coords, key = lambda x: x[0])[0]-min(all_coords, key = lambda x: x[0])[0], \
									max(all_coords, key = lambda x: x[1])[1]-min(all_coords, key = lambda x: x[1])[1], \
									max(all_coords, key = lambda x: x[2])[2]-min(all_coords, key = lambda x: x[2])[2]
		return x_width, y_width, z_width

	def make_dic_coords(rdkitmol):
		"""
		:param rdkitmol: takes a rdkit mol
		:return: returns a dic formatted as {C:[[x1, y1, z1], [x2,y2,z2]], H:[[x3,y3,z3],[x4,y4,z4]]}
		"""
		atom_pos_dic ={}
		for i, at in enumerate(rdkitmol.GetAtoms()):
			symbol = at.GetSymbol()
			id = at.GetIdx()
			pos = list(rdkitmol.GetConformer(-1).GetAtomPosition(id))
			if symbol not in atom_pos_dic:
				atom_pos_dic[symbol] = [pos]
			else:
				atom_pos_dic[symbol].append(pos)
		return atom_pos_dic

	pad = 3
	name_ = "{0!s}\n"
	lattice_constant = "{0:.5f}\n"
	vector = '     {0:.6f}    {1:.6f}    {2:.6f}\n'
	number_of_atoms = "     {0:d}"
	typeprocess = "{0!s}\n"
	positionline = "  {0: .6f}  {1: .6f}  {2: .6f}   {3!s} {4!s} {5!s}\n"

	atm_dic = make_dic_coords(rdkitmol)
	x_width, y_width, z_width = get_box_size(atm_dic)

	if name[0]=='p':
		carbon_list, hydrogen_list = remap_molecules(rdkitmol, 0)
	elif name[0]=='m':
		carbon_list, hydrogen_list = remap_molecules(rdkitmol, 1)
	elif name[0]=='o':
		carbon_list, hydrogen_list = remap_molecules(rdkitmol, 2)

	atm_types = ['C', 'H']
	atm_numbers = [len(carbon_list), len(hydrogen_list)]

	final_string= [name_.format(name),
					lattice_constant.format(2),
					vector.format((x_width+pad*2)/2.0,0,0),
					vector.format(0,(y_width+2*pad)/2.0,0),
					vector.format(0,0,(z_width+2*pad)/2.0)]
	for els in atm_numbers:
		final_string.append(number_of_atoms.format(els))

	final_string.append(typeprocess.format('\ncart'))

	for pos in carbon_list:
		final_string.append(positionline.format(pos[0], pos[1], pos[2], 'T', 'T', 'T'))
	for pos in hydrogen_list:
		final_string.append(positionline.format(pos[0], pos[1], pos[2], 'T', 'T', 'T'))

	print ''.join(final_string)

	with open('/home/macenrola/Desktop/{}'.format(name), 'wb') as w:
		w.write(''.join(final_string))

	return

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

	Chem.MolToMolFile(mol1, '/home/macenrola/Desktop/oXylene.sdf')
	Chem.MolToMolFile(mol2, '/home/macenrola/Desktop/mXylene.sdf')
	Chem.MolToMolFile(mol3, '/home/macenrola/Desktop/pXylene.sdf')
	return
# def align_protonated_xylene_according_to_methyl(rdkitmol):
# 	"""
# 	:param rdkitmol: Aligns the protonated xylenes according to the two methyl groups in the molecule, the protonated hydrogen is placed in negative x domain
# 	:return: a rdkit mol with adapted geometries
# 	"""
# 	def get_atom_coords_for_projection(rdkitmol):
# 		"""
# 		:param rdkitmol: the very same rdkit mol
# 		:return: the atom coords that one would use for the projection are returned, that is the methyl carbon not involved in the reaction and the other sp2 carbons
# 		"""
# 		mol = rdkitmol
# 		Chem.SanitizeMol(mol)
# 		mol.UpdatePropertyCache(strict=False)
# 		pattpos = Chem.MolFromSmarts('[CH3][CH]')
# 		pos = mol.GetSubstructMatches(pattpos)[0][-1]
# 		mol.GetAtomWithIdx(pos).SetFormalCharge(1)
# 		for at in mol.GetAtoms():
# 			at.SetNoImplicit(True)
# 			# print at.GetExplicitValence(), at.GetImplicitValence(), at.GetFormalCharge(), at.GetIsAromatic(), at.GetIdx(), at.GetHybridization()
# 		# pattalign = Chem.MolFromSmarts('[CH3]~[CH0]~[C;!+]~[C;!+]~[C]')
# 		# return mol.GetSubstructMatches(pattalign)[0]
# 		pattalign = Chem.MolFromSmarts('[CR]')
# 		return [x[0] for x in mol.GetSubstructMatches(pattalign)]
#
#
# 	def generate_mol_from_MDL(RDKIT_BLOCK_IN, coord_matrix):
# 		"""Will write the MDL of the mol file then replace the xyz coordinates from the coord_matrix"""
# 		RDKIT_BLOCK = [x+'\n' for x in RDKIT_BLOCK_IN.split('\n')]
# 		atm_number = int(RDKIT_BLOCK[3][:3])
# 		for i in range(0,atm_number):
# 			j = i+4
# 			RDKIT_BLOCK[j] = RDKIT_BLOCK[j].split()
# 			RDKIT_BLOCK[j][:3] = coord_matrix[i, :]
# 			RDKIT_BLOCK[j] = (' '*(3+int(np.sign(RDKIT_BLOCK[j][0])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][0])+
# 					' '*(3+int(np.sign(RDKIT_BLOCK[j][1])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][1])+
# 					' '*(3+int(np.sign(RDKIT_BLOCK[j][2])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][2])+
# 					' {}   '.format(RDKIT_BLOCK[j][3]) + '  '.join(RDKIT_BLOCK[j][4:]) + '\n'
# 					)
#
# 		RDKIT_BLOCK_OUT = ''.join(RDKIT_BLOCK)
#
# 		return RDKIT_BLOCK_OUT
#
# 	def get_fit_coords_all_coords(rdkitmol, fitIdx):
# 		"""
# 		:param rdkitmol: a regular rdkit mol with the Hydrogens atoms on
# 		:return: the coordinates to perform the pca alignment on and all the coordinates to be aligned, the coordinates are given in the same order the rdkit supplies the atoms (same as mol block)
# 		"""
# 		fitcoords = []
# 		allcoords = []
# 		for at in rdkitmol.GetAtoms():
# 			id = at.GetIdx()
# 			pos = list(rdkitmol.GetConformer(-1).GetAtomPosition(id))
# 			if id in fitIdx:
# 				fitcoords.append(pos)
# 			allcoords.append(pos)
# 		return fitcoords, allcoords
#
# 	def fit_and_transform_coords(fitcoords, allcoords):
# 		"""
# 		:param fitcoords: takes in the coordinates for the pca to be fitted on
# 		:param allcoords: takes in the coordinates on which the alignment needs to be performed
# 		:return: the transformed allcoords data according to the fit obtained from fitcoords
# 		"""
# 		from sklearn.decomposition import PCA
# 		pca = PCA(n_components=3)
# 		pca.fit(fitcoords)
# 		transformcoords = pca.transform(allcoords)
# 		transform_fitcoords = pca.transform(fit_coords)
# 		transform_coord_centered = transformcoords.copy()
# 		transform_coord_centered[:,0] = transformcoords[:,0] - np.mean(transform_fitcoords[:,0])
# 		transform_coord_centered[:,1] = transformcoords[:,1] - np.mean(transform_fitcoords[:,1])
# 		transform_coord_centered[:,2] = transformcoords[:,2] - np.mean(transform_fitcoords[:,2])
# 		return transformcoords
#
# 	coordsForFit = get_atom_coords_for_projection(rdkitmol)
# 	print coordsForFit
# 	fit_coords, allcoords = get_fit_coords_all_coords(rdkitmol, coordsForFit)
# 	transformcoords = fit_and_transform_coords(fit_coords, allcoords)
# 	print allcoords, transformcoords, fit_coords
# 	Chem.SanitizeMol(rdkitmol)
# 	rdkitmol.UpdatePropertyCache()
# 	molblock = generate_mol_from_MDL(Chem.MolToMolBlock(rdkitmol), transformcoords)
# 	mol = Chem.MolFromMolBlock(molblock, removeHs=False)
# 	return mol

if __name__ == "__main__":
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import numpy as np

	# oxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/pXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	# format_mol_for_vasp(oxylene, 'pXylene_protonated')
	# oxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/oXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	# mxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/mXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol',removeHs=False)
	# pxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/pXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol',	removeHs=False)
	# align_xylenes(oxylene, mxylene, pxylene)

	mxylene = Chem.MolFromMolFile('/home/macenrola/Desktop/mXylene.sdf', removeHs=False, sanitize=False)
	pxylene = Chem.MolFromMolFile('/home/macenrola/Desktop/pXylene.sdf', removeHs=False, sanitize=False)
	oxylene = Chem.MolFromMolFile('/home/macenrola/Desktop/oXylene.sdf', removeHs=False, sanitize=False)

	format_mol_for_vasp(mxylene, 'mXylene-Protonated')
	format_mol_for_vasp(pxylene, 'pXylene-Protonated')
	format_mol_for_vasp(oxylene, 'oXylene-Protonated')
	# Chem.MolToMolFile(align_protonated_xylene_according_to_methyl(mol), '/home/macenrola/Desktop/pxylene.sdf')