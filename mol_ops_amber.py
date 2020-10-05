import rdkit 
from rdkit import Chem
from rdkit.Chem import AllChem
import scipy
from sklearn.decomposition import PCA
import numpy as np
import sys, glob

"""
We're following this tutorial
http://ambermd.org/tutorials/advanced/tutorial21/index.html
"""
def make_part_tleap_inputs(rdid, TAG, solvate=True, old_tleap_solvate_log=""):
	"""
	PRE  : there exist a frcmod and mol2 file that have been obtained using
	frcmod_cmd = "{0}/parmchk2 -i outputs/{1}.mol2 -f mol2 -o outputs/{1}.frcmod -f frcmod".format(antechamber_path, rdid).split()
	get_mol2_cmd = "{0}/antechamber -i outputs/{1}.pdb -fi pdb -o outputs/{1}.mol2 -fo mol2 -c bcc -nc {2} -s 2 -j 4 -pl 15".format(antechamber_path, rdid, charge)
	POST :
	"""
	solv_nums = {"CB7":11, "MVp":13, "GST":15}
	tleap_file = "{}.in".format(rdid)
	if solvate:
		
		solv_line = "solvatebox {} TIP5PBOX {} iso".format(TAG, solv_nums[TAG])
	else:
		solv_line =""
	if old_tleap_solvate_log=="":
		remove_water=""
	else:
		with open(old_tleap_solvate_log, "r") as r:
			for line in r:
				if "WAT" in line: 
					wat_nbr = int(line.split()[-1])
					print rdid, wat_nbr
		remove_water = "\n".join(["remove {0} {0}.{1}".format(TAG, x) for x in range(1501, wat_nbr+1)[::-1]])
	prmtop_script="""
source leaprc.gaff
source oldff/leaprc.ff14SB
loadAmberParams frcmod.tip5p
loadAmberParams {0}.mol2.frcmod
{1} = loadMol2 {0}.mol2
check {1}
loadamberparams gaff.dat
{2}
{3}
saveamberparm {1} {0}.prmtop {0}.rst7
saveoff {1} {0}.lib
quit
""".format(rdid.split("/")[-1], TAG, solv_line, remove_water)
	with open(tleap_file, "w") as w:
		w.write(prmtop_script)

	return tleap_file
def make_complex_tleap_script(rdids_tags, complex_pdb, solvate=True, old_tleap_solvate_log=""):
	"""
	Parameters
	----------
	rdids_tags : tuple
		of the shape [(rdid1, TAG1), (rdid2, TAG2), ...]
		it supposes there are lib, frcmod, and mol2 files associated with these 

	Returns
	-------
	a tleap file to merge them together

	"""
	if solvate:
		solv_line = "solvatebox {} TIP5PBOX 11 iso".format("CMP")
	else:
		solv_line =""
	if old_tleap_solvate_log=="":
		remove_water=""
	else:
		with open(old_tleap_solvate_log, "r") as r:
			for line in r:
				if "WAT" in line: wat_nbr = int(line.split()[-1])
		remove_water = "\n".join(["remove {0} {0}.{1}".format("CMP", x) for x in range(1501, wat_nbr+1)[::-1]])
	tleap_complex_file = "{}-CMP.in".format(rdids_tags[0][0])

	combine_script = """
source leaprc.gaff
source oldff/leaprc.ff14SB
loadAmberParams frcmod.tip5p
{0}
{1}
{2}
CMP = loadpdb {4}
savemol2 CMP {3}.mol2 1
{5}
{6}
saveamberparm CMP {3}.prmtop {3}.rst7
quit
""".format("\n".join(["loadoff {0}.lib".format(x[0].split("/")[-1]) for x in rdids_tags]), 
			"\n".join(["loadamberparams {0}.mol2.frcmod".format(x[0].split("/")[-1]) for x in rdids_tags]),
			"\n".join(["{1} = loadmol2 {0}.mol2".format(x[0].split("/")[-1], x[1]) for x in rdids_tags]), rdids_tags[-1][0].split("/")[-1].replace("MV", "CMP"), complex_pdb.split("/")[-1],
			solv_line, remove_water)
	
	with open(tleap_complex_file, "w") as w:
		w.write(combine_script)
		
	return tleap_complex_file


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


def make_pdb_with_named_residues(base_pdb, TAG, CAPS=False, residue_number={"CB7":1, "MVp":3, "GST":2}):
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

	def fix_PDB_spacing(fname, residue_number=1):
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
					line_content[5]=str(residue_number)
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
	output_pdb = base_pdb[:-4]+"-named.pdb"
# =============================================================================
# 	with open(base_pdb, "r") as r:
# 		print r.readlines()
# =============================================================================
	mol = Chem.MolFromPDBFile(base_pdb, removeHs=False)
	print mol
	frags = Chem.GetMolFrags(mol, asMols=True)
	if len(frags)==1:
		renumber_pdb_file(base_pdb, TAG, output_pdb)
		fix_PDB_spacing(output_pdb, residue_number[TAG])
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

def make_complex_named_pdb_from_named_parts(pdbfiles_to_stack):
	"""
	PRE: Takes in n files with named pdbs and will stack them without connect lines and call it a complex 
	"""
	with open(pdbfiles_to_stack[0]+"-COMPLEX.pdb", 'w') as w:
		w.write("COMPND    {}\n".format(pdbfiles_to_stack[0].split("/")[-1]))
		for els in pdbfiles_to_stack:
			with open(els , 'r') as r:
				for line in r:
					if "ATOM" in line:
						w.write(line)
		w.write("END\n")
	return 

def get_equilibration_script(rdid_path, bin_path=""):
	"""
	PRE: There exist all the relevant min files 
	POST: Generates the script to run the minimizations
	"""
	rdid = rdid_path.split("/")[-1]
	script = []
	script.append('# Minimization with restraints of the different parts\n')
	script.append('{0}pmemd.MPI -O -i restrain.in -p {1}.prmtop -c {1}_iso.rst7 -o {1}_min.out -r {1}_min.rst7 -ref {1}_iso.rst7 -x {1}_iso.nc'.format(bin_path, rdid ))

	#pmemd.MPI or sander 
	script.append('# Minimization without restraints of the different parts\n')
	script.append('{0}pmemd.MPI -O -i all_min.in -p {1}.prmtop -c {1}_min.rst7 -o {1}_min2.out -r {1}_min2.rst7 -x {1}_all.nc'.format(bin_path, rdid))

	script.append('# Minimization with NVT\n')
	script.append('{0}pmemd.MPI -O -i nvt_md.in -p {1}.prmtop  -c {1}_min2.rst7 -o {1}_nvt.out -r {1}_nvt.rst7 -x {1}_nvt.nc'.format(bin_path, rdid))

	script.append('# Minimization with MD1\n')
	script.append('{0}pmemd.MPI -O -i npt_md.in -p {1}.prmtop  -c {1}_nvt.rst7  -o {1}_npt.out -r {1}_npt.rst7 -x {1}_npt.nc'.format(bin_path, rdid))

	script.append('# Minimization with MD2\n')
	script.append('{0}pmemd.MPI -O -i npt_md2.in -p {1}.prmtop -c {1}_npt.rst7  -o {1}_npt2.out -r {1}_npt2.rst7 -x {1}_npt2.nc'.format(bin_path, rdid))

	script.append('# Average the volumes\n')
	script.append('grep VOL {}_npt2.out | head -n -2 | awk \'{{sum+=$9}}END{{printf \"%10.7f\\n\",(sum/NR)^(1/3)}}\''.format(rdid))

	with open(rdid_path+"_prepa.sh", "w") as w:
		w.write("\n".join(script))
		

if __name__ == "__main__":
	pass
	#make_pdb_with_named_residue("/home/macenrola/Documents/ML/ChemTS/new_scoring_for_mcts/docking_targets/adamantanone-docked-named-opti.pdb", "CMP") # JUST AN EXAMPLE

	"""
	1) get the named residues, otherwise amber might crash unexpectedly 
	2) from the named pdbs obtain the mol2 and frcmod using antechamber and parmchk2
	3) get tleap input files to create topologies for the parts of each complex
	4) SINCE we TARGET 1500 molecules in each case, run tleap to get solvated structures
	5) Use the log from this tleap run to know how many water molecules to remove
	6) rerun tleap by removing a predicable amount of water
	7) build the complex as a stack, unconnected amount of named pdbs
	8) repeat steps 3) to 6) for the complex
	"""
# =============================================================================
# 	for f in glob.glob("/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/*.pdb"):
# 		if "MV." in f:
# 			tag = 'MVp'
# 		if "CB." in f :
# 			tag = "CB7"
# 		if "G2." in f:
# 			tag = "GST"
# 		make_pdb_with_named_residues(f, tag) # replace CB7, GST, MVp
# =============================================================================
# =============================================================================
# 	# GET MOL2, FRCMOD
# 	for f in glob.glob("/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/*named.pdb"):
# 		if "MV-named." in f:
# 			tag = 'MVp'
# 		if "CB-named." in f :
# 			tag = "CB7"
# 		if "G2-named." in f:
# 			tag = "GST"
# 		make_part_tleap_inputs(f, tag, old_tleap_solvate_log="") # replace CB7, GST, MVp
# 		#make_part_tleap_inputs(f, tag, old_tleap_solvate_log=f+".in.log") # replace CB7, GST, MVp
# =============================================================================
# =============================================================================
# 	# GET prmtop, rst7, lib for parts
# 	pdbs = sorted(glob.glob("/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/*named.pdb"))
# 	for els in range(0, len(pdbs), 3):
# 		print pdbs[els:els+3]
# 		make_complex_named_pdb_from_named_parts(pdbs[els:els+3])
# 		#make_complex_tleap_script(zip(pdbs[els:els+3], ["CB7",  "GST", "MVp",]), pdbs[els]+"-COMPLEX.pdb")
# 		make_complex_tleap_script(zip(pdbs[els:els+3], ["CB7",  "GST", "MVp",]), pdbs[els]+"-COMPLEX.pdb", old_tleap_solvate_log=pdbs[els]+"-CMP.in.log")
# 	# GET prmtop, rst7, lib for complex
# =============================================================================
	
# =============================================================================
# 	# Now build the system for just CB and MVp 
# 	#make_complex_named_pdb_from_named_parts(["/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_AN_sol_3_CB-named.pdb",
# 	#									  "/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_AN_sol_3_MV-named.pdb"])
# 	#make_complex_tleap_script(zip(["/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_AN_sol_3_CB-named.pdb",
# 	#									  "/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_AN_sol_3_MV-named.pdb"], ["CB7", "MVp",]), "/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_CB_BINARY-COMPLEX.pdb")
# 	#make_complex_tleap_script(zip(["/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_AN_sol_3_CB-named.pdb",
# 	#									  "/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_AN_sol_3_MV-named.pdb"], ["CB7", "MVp",]), "/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_CB_BINARY-COMPLEX.pdb",
# 	#					   old_tleap_solvate_log="/home/macenrola/Documents/water_in_cavity_CB8/CB8-MV-G2_CPCM/CB8_MV_AN_sol_3_CB-named.pdb-CMP.in.log")
# =============================================================================
	for f in glob.glob("/home/macenrola/Documents/water_in_cavity_CB8/tests/*.prmtop"):
		rdid = f[:-7]
		get_equilibration_script(rdid)
