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
					# vector.format((x_width+pad*2)/2.0,0,0),
					# vector.format(0,(y_width+2*pad)/2.0,0),
					# vector.format(0,0,(z_width+2*pad)/2.0)
				   vector.format(13/2., 0, 0),
				   vector.format(0, 12. / 2.0, 0),
				   vector.format(0, 0, 10 / 2.0)

				   ]
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
	for i,m in enumerate([mol1, mol2, mol3]):
		print i
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

	Chem.MolToMolFile(mol1, '/home/macenrola/Desktop/oXylene_align.sdf')
	Chem.MolToMolFile(mol2, '/home/macenrola/Desktop/mXylene_align.sdf')
	Chem.MolToMolFile(mol3, '/home/macenrola/Desktop/pXylene_align.sdf')
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

def get_CB_guest_atomxyz(rdkitmol):
	"""
	:param rdkitmol: takes in a rdkit_mol
	:return: returns the atomic coordinates for the cb and the guest
	"""
	try:
		complex, guest = Chem.GetMolFrags(rdkitmol, asMols=True)
		if complex.GetNumAtoms()<guest.GetNumAtoms():
			complex, guest = guest, complex

		#### Get poits for the the hull convex
		hull_points = []
		complexc = complex.GetConformer(-1)
		for i in range(complexc.GetNumAtoms()):
			cc = [complex.GetAtomWithIdx(i).GetSymbol()]
			cc.extend(list(complexc.GetAtomPosition(i)))
			hull_points.append(cc)
			# print cc
	except:
		hull_points = []
		guest = rdkitmol
	#### Guest the guest points
	guest_points = []
	guestc = guest.GetConformer(-1)
	for i in range(guestc.GetNumAtoms()):
		cc = [guest.GetAtomWithIdx(i).GetSymbol()]
		cc.extend(list(guestc.GetAtomPosition(i)))
		guest_points.append(cc)
		# print cc

	return hull_points, guest_points


def make_nw_paramfile(inpdbfile):
	"""
	:param inpdbfile: Takes in a PDB file and will print the corresponding nwchem file for minimization
	:return:
	"""
	mol = Chem.MolFromPDBFile(inpdbfile, removeHs=False)
	cb_points, guest_points = get_CB_guest_atomxyz(mol)
	print cb_points, guest_points

	start_script = """
echo
start {0}
charge 1
geometry units angstrom nocenter noautosym noautoz\n"""
	positionline = " {0} {1: .6f}  {2: .6f}  {3: .6f}\n"
	end_script = """
basis
 C library 6-31g*
 H library 6-31g*
end
dft
 xc pbe0
 disp vdw 3
 iterations 300
end
driver
 maxiter 300 
end
task dft optimize
"""
	fname = inpdbfile[:-4]+'.nw'
	name = inpdbfile.strip().split('/')[-1][:-4]+'_nw'
	script = start_script.format(name)
	for p in cb_points+guest_points:
		script += positionline.format(*p)
	script += end_script
	with open(fname, 'wb') as w:
		w.write(script)
	script_launch = """
#!/bin/bash -l

# Batch script to run an MPI NWChem job on Legion with the upgraded 
# software stack under SGE. Oct 2015

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request thirty minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=12:00:0

# 3. Request 1 gigabyte of RAM.
#$ -l mem=1G

# 4. Set the name of the job.
#$ -N {0}

# 5. Select the MPI parallel environment and 8 processors.
#$ -pe mpi 64

# 7. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to your $HOME.
#
# NOTE: this directory must exist.
#
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -cwd

# 8. Now we need to set up and run our job. 

module load python/2.7.12
module load nwchem/6.8-47-gdf6c956/intel-2017

module list

# $NSLOTS will get the number of processes you asked for with -pe mpi.
mpirun -np $NSLOTS -machinefile $TMPDIR/machines nwchem {1}
	"""
	print script_launch.format(name, fname.strip().split('/')[-1])
	with open(fname[:-3]+'.sh', 'wb') as w:
		w.write(script_launch.format(name, fname.strip().split('/')[-1]))


def make_gaussian_input_files_for_molgroup(list_of_pdbs):
	"""
	:param list_of_pdbs: take in a list of pdb files located into the very same directory
	:return: produce a com file for each of the pdb files, a paramfile and a sh file to launch gaussian
	"""
	path = '/'.join(list_of_pdbs[0].split('/')[:-1]) + '/'
	genname = 'GUEST'
	paramfile = path + '{}paramfile'.format(genname)
	nproc = 12
	memperproc = 2
	shfile = path + '{}.sh'.format(genname)
	comstring = "\
%Chk={0}.chk\n\
%NProcShared={1}\n\
%mem={2}gb\n\
#n wB97XD/6-31G* Opt\n\
\n\
{0}\n\
\n\
0 1\n\
"
	shstring = "#!/bin/bash -l\n\
#$ -S /bin/bash\n\
#$ -cwd\n\
#$ -l h_rt=3:0:0\n\
#$ -l mem=%iG\n\
#$ -l tmpfs=100G\n\
#$ -N %s\n\
#$ -pe smp %i\n\
#$ -t 1-%i\n\
number=$SGE_TASK_ID\n\
paramfile=$(pwd)/%s\n\
index=$(sed -n ${number}p $paramfile | awk '{print $1}')\n\
g16infile=$(sed -n ${number}p $paramfile | awk '{print $2}')\n\
g16outfile=$g16infile'_OUT.out'\n\
module load gaussian/g16-a03/pgi-2016.5\n\
source $g16root/g16/bsd/g16.profile\n\
mkdir -p $GAUSS_SCRDIR\n\
# echo $g16root $GAUSS_SCRDIR $GAUSS_EXEDIR $lindaConv\n\
# echo \"GAUSS_SCRDIR = $GAUSS_SCRDIR\"\n\
# echo \"Running: lindaConv $g16infile $JOB_ID $TMPDIR/machines\"\n\
# echo ""\n\
# $lindaConv $g16infile $JOB_ID $TMPDIR/machines\n\
# echo \"Running: g16 < $g16infile > $g16outfile\"\n\
# export GAUSS_MEMDEF=76MW\n\
# export GAUSS_LFLAGS='-vv -opt \"Tsnet.Node.lindarsharg: ssh\"'\n\
# time g16 < job$JOB_ID.com > $g16outfile\n\
time g16 < $g16infile > $g16outfile\n\
	"
	posstring = '{0}         {1: .5f}       {2: .5f}       {3: .5f}\n'
	print paramfile
	with open(paramfile, 'wb') as w: pass
	with open(paramfile, 'ab') as a:
		for i, fname in enumerate(sorted(list_of_pdbs)):
			comname = fname[:-4] + '.com'
			a.write('{0:04d}\t{1}\n'.format(i+1, comname.split('/')[-1]))

			cmol = Chem.MolFromMolFile(fname, removeHs=False)
			ccomplex, cguest = get_CB_guest_atomxyz(cmol)
			ccomstring = comstring.format(fname.split('/')[-1][:-4], nproc, nproc*memperproc)
			for coord in ccomplex+cguest:
				ccomstring = ccomstring+posstring.format(*coord)
			ccomstring += '\n'

			with open(comname, 'wb') as w:
				w.write(ccomstring)

	with open(shfile, 'wb') as w:
		w.write(shstring%(memperproc, genname,nproc,i+1, paramfile.split('/')[-1]))


if __name__ == "__main__":
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import numpy as np
	import glob

	# oxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/pXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	# format_mol_for_vasp(oxylene, 'pXylene_protonated')
	# oxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/oXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	# mxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/mXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol',removeHs=False)
	# pxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/pXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol',	removeHs=False)
	# align_xylenes(oxylene, mxylene, pxylene)

	# oxylene = Chem.MolFromPDBFile('/home/macenrola/Documents/nwchem/OXYLENE_OUT.pdb', removeHs=False, sanitize=False)
	# pxylene = Chem.MolFromPDBFile('/home/macenrola/Documents/nwchem/PXYLENE_OUT.pdb', removeHs=False, sanitize=False)
	# mxylene = Chem.MolFromPDBFile('/home/macenrola/Documents/nwchem/MXYLENE_OUT.pdb', removeHs=False, sanitize=False)
	# print oxylene, pxylene, mxylene
	# align_xylenes(oxylene, mxylene, pxylene)
	flist = glob.glob('/home/macenrola/Documents/amberconvergedmols/top50/MINE/*GUEST.sdf')
	make_gaussian_input_files_for_molgroup(flist)
