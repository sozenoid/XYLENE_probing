import rdkit 
from rdkit import Chem
from rdkit.Chem import AllChem
import scipy
from sklearn.decomposition import PCA
import numpy as np

"""
We're following this tutorial
http://ambermd.org/tutorials/advanced/tutorial21/index.html
"""

def get_atoms_coords(RDKIT_BLOCK):
	"""Takes as input an RDKIT BLOCK and returns a list of atoms with a numpy array containing the coordinates"""
	RDKIT_BLOCK = RDKIT_BLOCK.split('\n')
	atm_number = int(RDKIT_BLOCK[3][:3])
	RDKIT_BLOCK = [x.split() for x in RDKIT_BLOCK]
	atm_list = []
	coords_array = np.zeros([atm_number, 3], dtype=float)
	for i, line in enumerate(RDKIT_BLOCK[4:4+atm_number]):
	    coords_atm = line
	    atm_list.append(coords_atm[3])
	    coords_array[i, :] = coords_atm[:3]
	return atm_list, coords_array


def align_mol(RDKIT_BLOCK):
	"""
	PRE: Takes as input a RDKIT_BLOCK with valid 3D coordinates (basically a SDF file in text format)
	POST: Returns a RDKIT_BLOCK with the 3D coordinates rotated so that the main axis of the molecule (PCA wise) is aligned with the z axis and all coordinates are 0 centered
	"""
	def align_xyz(list_atoms, coord_matrix):
		"""This method uses a PCA on the atoms of the molecules as represented by the nuclear coordinates as points 
		and their electron cloud as a sphere of points around those nuclear coordinates
		The function returnCircleAsTuple is imported via a pybind11 generated module from c++ code"""
		sphere_points = get_sphere_points(list_atoms, coord_matrix)
		total_matrix = np.concatenate((coord_matrix, sphere_points),axis = 0)
		pca = PCA(n_components = 3)

		transform = pca.fit_transform(total_matrix)
		transform_coord = pca.transform(coord_matrix) 

		# point_cloud = zip(transform.T[1][:], transform.T[2][:])
		# height = np.max(transform.T[0][:]) - np.min(transform.T[0][:])
		# rad = make_circle(point_cloud)
		# rad = returnCircleAsTuple(point_cloud)

		transform_coord_centered = transform_coord.copy()
		transform_coord_centered[:,0] = transform_coord[:,0] - np.mean(transform_coord[:,0])
		transform_coord_centered[:,1] = transform_coord[:,1] - np.mean(transform_coord[:,1])
		transform_coord_centered[:,2] = transform_coord[:,2] - np.mean(transform_coord[:,2])
		return transform_coord_centered

	def get_sphere_points(list_atoms, coord_matrix):
		"""Find the thinnest cylinder dimesions where the molecule could fit add clouds of points around the atoms before the PCA (pi/3 in every direction with atomic radius)
		x = r*sin(theta)*cos(theta)
		y = r*sin(theta)*sin(theta)        theta is inclination from top (0 to pi) and phi is azimuth (0 to 2pi)
		z = r*cos(theta)
		This method uses list comprehension and return all the points representing the spheres around the atoms
		"""
		index_of_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06} 
		angle_nbs = 6 # for angle_nbs = 3 there will be points spaced by pi/3 so three sections in inclinations and six in azimutal 
		sphere_points = []
		for i in xrange(len(list_atoms)):
			radius = index_of_radii[list_atoms[i]]
			top_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(0)]
			bottom_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(np.pi)]
			sphere_points.append(top_sphere_point)
			sphere_points.append(bottom_sphere_point)


		new_sphere_points = [[coord_matrix[i][0]+index_of_radii[list_atoms[i]]*np.sin(inclination_angle)*np.cos(azymuth_angle), coord_matrix[i][1]+index_of_radii[list_atoms[i]]*np.sin(inclination_angle)*np.sin(azymuth_angle), coord_matrix[i][2]+index_of_radii[list_atoms[i]]*np.cos(inclination_angle)] 
							for azymuth_angle in np.linspace(0, np.pi *2, angle_nbs*2) for inclination_angle in np.linspace(0, np.pi, angle_nbs + 1)[1:-1] for i in range(len(list_atoms))]
		sphere_points.extend(new_sphere_points)
		return sphere_points

	def generate_mol_from_MDL(RDKIT_BLOCK_IN, coord_matrix):
		"""Will write the MDL of the mol file then replace the xyz coordinates from the coord_matrix"""
		RDKIT_BLOCK = [x+'\n' for x in RDKIT_BLOCK_IN.split('\n')]
		atm_number = int(RDKIT_BLOCK[3][:3])
		for i in range(0,atm_number):
			j = i+4
			RDKIT_BLOCK[j] = RDKIT_BLOCK[j].split()
			RDKIT_BLOCK[j][:3] = coord_matrix[i, :]
			RDKIT_BLOCK[j] = (' '*(3+int(np.sign(RDKIT_BLOCK[j][0])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][0])+
					' '*(3+int(np.sign(RDKIT_BLOCK[j][1])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][1])+
					' '*(3+int(np.sign(RDKIT_BLOCK[j][2])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][2])+
					' {}   '.format(RDKIT_BLOCK[j][3]) + '  '.join(RDKIT_BLOCK[j][4:]) + '\n'
					)

		RDKIT_BLOCK_OUT = ''.join(RDKIT_BLOCK)

		return RDKIT_BLOCK_OUT  

	atom_coords = get_atoms_coords(RDKIT_BLOCK)
	transformed_coords = align_xyz(atom_coords[0], atom_coords[1])
	return generate_mol_from_MDL(RDKIT_BLOCK, transformed_coords)


def make_pdb_with_named_residue(base_pdb, TAG, CAPS=False):
	"""
	PRE: Takes a RDKIT_BLOCK with valid 3D coordinates (basically a SDF file in text format), centered around 0 and with principal axis (in the PCA sense) aligned along the z axis
	POST: Creates three PDB files with properly named residues (default GST for the guest, CB7 for CB7) and connection records to be used in AMBER
	"""
	def remove_CONNECT_LINES(fname):
		"""
		PRE: fname contains a valid PDB file of a molecule with 3D coordinates
		POST: The lines containing CONNECT and MASTER are removed from the original file, the original file IS MODIFIED
		"""
		with open(fname, 'rb') as r:
			lines = r.readlines()
		with open(fname, 'wb') as w:
			w.writelines([x for x in lines if 'CONECT' not in x if 'MASTER' not in x][1:])

	def fix_PDB_spacing(fname):
		"""
		PRE: The PDB files is formated by rdkit MolToPDBFile flavour 28 without any MASTER, TER, CONECT or Charge flags, 
		POST: The file IS MODIFIED to be formated with the adequate amount of blank space to be read by AMBER, This method destroys the original file
		"""
		raw_spc = [7, 2, 4, 4, 4, 12, 8, 8, 6, 6, 12, 0]
		new_lines = []
		cb_yet = False
		with open(fname, 'rb') as r:
			for line in r:
				if 'CB7' in line and not cb_yet:
					cb_yet = True
					new_lines.append('TER')
				if 'ATOM' in line:
					line_content = line.split()
					line_content.insert(4, 'A')
					# print line_content
					ls = [len(x) for x in line_content]
					actual_spacing = [raw_spc[0]-ls[1], # after ATOM 
					raw_spc[1],  # after ATOM# (11)
					raw_spc[2]-ls[2], # after ATOM NAME (H41)
					raw_spc[3]-ls[3], # after RESIDUE NAME (CUC)
					raw_spc[4]-ls[4], # after CHAIN ID (A)
					raw_spc[5]-ls[6], # after RESIDUE NUMBER (1)
					raw_spc[6]-ls[7], # after X CART COORDINATE (6.171)
					raw_spc[7]-ls[8], # after Y CART COORDINATE (3.377)
					raw_spc[8]-ls[9], # after Z CART COORDINATE (21.096)
					raw_spc[9]-ls[10], # after enigmatic number (1.00)
					raw_spc[10]-ls[11], # after partial charge (0.00)
					raw_spc[11], # after ATOMIC SYMBOL (N)
					]
					if CAPS == True:
						## FIX CAPS CHAOS, antechamber uses full caps element names but rdkit cannot read these
						# Typical line like this ATOM      1  C1  GST A   1      -4.939  -0.132   1.572  1.00  0.00           C
						caps_name = list(line_content[2])
						caps_name[:len(line_content[-1])] = line_content[-1]
						line_content[2] = "".join(caps_name)
						##
					new_lines.append(''.join([x[0]+' '*x[1] for x in zip(line_content, actual_spacing)]))
					print new_lines[-1]

				else:
					new_lines.append(line.strip())
					print new_lines[-1]
		with open(fname, 'wb') as w:
			w.writelines('\n'.join(new_lines)+'\n')

	def renumber_pdb_file(pdbfile, TAG, outputpdb):
		"""

		Parameters
		----------
		pdbfile : TYPE PDB file
			DESCRIPTION. Takes  in a pdb file and 

		Returns: will produce a new pdb file 
			that is suited to work with amber
		-------
		None.

		"""
		mol = Chem.MolFromPDBFile(pdbfile, removeHs=False)
		atm_dic = {}
		for atom in mol.GetAtoms():
			if atom.GetSymbol() not in atm_dic:
				atm_dic[atom.GetSymbol()] = 1
			else: atm_dic[atom.GetSymbol()] += 1
			atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} {}'.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10), TAG)))
			atom.GetMonomerInfo().SetResidueNumber(1)
		Chem.MolToPDBFile(mol, outputpdb, flavor=flavour)
		return 
	
	flavour = 28
	output_pdb = base_pdb#[:-4]+"-named.pdb"
# =============================================================================
# 	with open(base_pdb, "r") as r:
# 		print r.readlines()
# =============================================================================
	mol = Chem.MolFromPDBFile(base_pdb, removeHs=False)
	frags = Chem.GetMolFrags(mol, asMols=True)
	if True:#len(frags)==1:
		renumber_pdb_file(base_pdb, TAG, output_pdb)
		fix_PDB_spacing(output_pdb)
	else:
		mol.SetProp('_Name', 'CMP')
		Chem.MolToPDBFile(mol, output_pdb, flavor=flavour)
		remove_CONNECT_LINES(output_pdb)
		fix_PDB_spacing(output_pdb)
		
# =============================================================================
# 	CB = Chem.MolFromMolBlock(get_CB_BLOCK(), removeHs=False)
# 	atm_dic = {}
# 	for atom in CB.GetAtoms():
# 		if atom.GetSymbol() not in atm_dic:
# 			atm_dic[atom.GetSymbol()] = 1
# 		else: atm_dic[atom.GetSymbol()] += 1
# 		atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} CB7'.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10))))
# 		atom.GetMonomerInfo().SetResidueNumber(2)
# 
# 	Chem.MolToPDBFile(CB, pdb_file_cb, flavor=flavour)
# 	fix_PDB_spacing(pdb_file_cb)
# 
# 
# 	complex_CB_host = Chem.CombineMols(guest, CB)
# 	yo = Chem.GetMolFrags(complex_CB_host)
# 	for els in yo:
# 		print els
# 	complex_CB_host.SetProp('_Name', 'CMP')
# 	Chem.MolToPDBFile(complex_CB_host, pdb_file_complex, flavor=flavour)
# 	remove_CONNECT_LINES(pdb_file_complex)
# 	fix_PDB_spacing(pdb_file_complex)
# =============================================================================

	return 


if __name__ == "__main__":
	pass
	#make_pdb_with_named_residue("/home/macenrola/Documents/ML/ChemTS/new_scoring_for_mcts/docking_targets/adamantanone-docked-named-opti.pdb", "CMP")

