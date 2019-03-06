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
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
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
		print 'ayayay, this is not a complex'
	#### Guest the guest points
	guest_points = []
	guestc = guest.GetConformer(-1)
	for i in range(guestc.GetNumAtoms()):
		cc = [guest.GetAtomWithIdx(i).GetSymbol()]
		cc.extend(list(guestc.GetAtomPosition(i)))
		guest_points.append(cc)
		# print cc

	return hull_points, guest_points


def make_nw_paramfile(list_of_sdfs):
	"""
	:param inpdbfile: Takes in a list of SDF files
	:return: will write a .nw and an .sh file for all of them to be run on Stratus
	"""
	path = '/'.join(list_of_sdfs[0].split('/')[:-1]) + '/'
	start_script = "\
echo\n\
start {0}\n\
charge 0\n\
geometry units angstrom noautoz\n"
	positionline = " {0} {1: .6f}  {2: .6f}  {3: .6f}\n"
	end_script = "\
end\n\
basis\n\
 * library 3-21g\n\
end\n\
dft\n\
 xc pbe0\n\
 disp vdw 3\n\
 iterations 300\n\
 direct\n\
 grid NODISK xfine\n\
end\n\
driver\n\
 maxiter 300\n\
end\n\
task dft optimize\n\
	"
	script_launch = """
#!/bin/bash
## -N = Job name							##
## select = No of nodes required, max 144 nodes				##
## ncpus = No of cores/node required, max 64cores/node			##
##         Request all 64 for exclusive use of the nodes for the job	##
## mpiprocs = How many MPI processes on each node			##
## ompthreads = How many OpenMP threads on each MPI process		##
##              Leave this to 1 for a pure MPI job			##
## walltime = HH:MM:SS							##
##									##
## Output and error file ##						##
## -o testjob.out							##
## -e testjob.err							##
##									##
## -q normal --> default queue, will route to (express,long or all)	##
## express = max 128 core/job, walltime 60 mins                 	##
## long = max 128 core/job, walltime 14 days                    	##
## all = max 9216 core/job, walltime 3 days                     	##
##									##
##									##
## Below example shows a (2 x 64) cores job using the Intel MPI library	##

#PBS -N {0}
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#PBS -l walltime=24:0:00
#PBS -o {0}.out
#PBS -e {0}.err
#PBS -q normal
#PBS -l software="nwchem"

cd $PBS_O_WORKDIR

module load intel/17.0.1.132
module load impi/2017_Update_1
module load nwchem/6.6-rev29223
# Now run the program
time mpirun nwchem {1} > {2}
	"""
	shstring = "\
#!/bin/bash -l\n\
#$ -S /bin/bash\n\
#$ -cwd\n\
#$ -l h_rt=12:10:0\n\
#$ -l mem=4G\n\
#$ -l tmpfs=100G\n\
#$ -N {}\n\
#$ -pe mpi 12\n\
module load python/2.7.12\n\
module unload compilers mpi\n\
module load compilers/intel/2017/update4\n\
module load mpi/intel/2017/update3/intel\n\
module load nwchem/6.8-47-gdf6c956/intel-2017\n\
module list\n\
mpirun -np $NSLOTS -machinefile $TMPDIR/machines nwchem {}\n\
	"
	for i, fname in enumerate(sorted(list_of_sdfs)):
		nwname = fname[:-4] + '.nw'
		shname = fname[:-4] + '.sh'
		genname = fname.split('/')[-1][:-4]
		cmol = Chem.MolFromMolFile(fname, removeHs=False)
		print cmol
		ccomplex, cguest = get_CB_guest_atomxyz(cmol)
		cscript = start_script.format(genname)
		for coord in ccomplex + cguest:
			cscript = cscript + positionline.format(*coord)
		cscript += end_script

		with open(nwname, 'wb') as w:
			w.write(cscript)
		with open(shname, 'wb') as w:
			w.write(shstring.format(genname, genname+'.nw',genname+'_LOGOUT'))


def make_neb_input_for_nwchem(pdb1, pdb2):
	"""
	:param pdb1: first pdb file for the initial configuration usually the fully docked version
	:param pdb2: second pdb file for the final configuration
	:note: try adding cgmin to the dft route if convergence sucks, for neutral molecules it seems to kind of okay though
	:return: a nwchem file and its launching scrip for a neb computation
	"""
	# path = '/'.join(pdb2.split('/')[:-1]) + '/'
	# genname = 'step_in'
	start_script = "\
echo\n\
start {0}\n\
charge 0\n\
geometry units angstrom nocenter noautosym noautoz\n"
	positionline = " {0} {1: .6f}  {2: .6f}  {3: .6f}\n"
	endgeom = "end\n\
geometry endgeom nocenter noautosym noautoz\n"
	end_script = "\
end\n\
basis\n\
 * library 3-21g\n\
end\n\
dft\n\
 xc pbe0\n\
 disp\n\
 iterations 300\n\
 direct\n\
 grid xfine NODISK\n\
 cgmin\n\
end\n\
driver\n\
 maxiter 300\n\
end\n\
neb\n\
  nbeads 20\n\
  kbeads 1.0\n\
  maxiter 100\n\
  stepsize 0.10\n\
  print_shift 1\n\
end\n\
task dft neb ignore\n\
	"
	script_launch = """
#!/bin/bash
## -N = Job name							##
## select = No of nodes required, max 144 nodes				##
## ncpus = No of cores/node required, max 64cores/node			##
##         Request all 64 for exclusive use of the nodes for the job	##
## mpiprocs = How many MPI processes on each node			##
## ompthreads = How many OpenMP threads on each MPI process		##
##              Leave this to 1 for a pure MPI job			##
## walltime = HH:MM:SS							##
##									##
## Output and error file ##						##
## -o testjob.out							##
## -e testjob.err							##
##									##
## -q normal --> default queue, will route to (express,long or all)	##
## express = max 128 core/job, walltime 60 mins                 	##
## long = max 128 core/job, walltime 14 days                    	##
## all = max 9216 core/job, walltime 3 days                     	##
##									##
##									##
## Below example shows a (2 x 64) cores job using the Intel MPI library	##

#PBS -N {0}
#PBS -l select=5:ncpus=64:mpiprocs=64:ompthreads=1
#PBS -l walltime=24:0:00
#PBS -o {0}.out
#PBS -e {0}.err
#PBS -q normal
#PBS -l software="nwchem"

cd $PBS_O_WORKDIR

module load intel/17.0.1.132
module load impi/2017_Update_1
module load nwchem/6.6-rev29223
# Now run the program
time mpirun nwchem {1} > {2}
	"""
	shstring = "\
#!/bin/bash -l\n\
#$ -S /bin/bash\n\
#$ -cwd\n\
#$ -l h_rt=24:10:0\n\
#$ -l mem=4G\n\
#$ -l tmpfs=100G\n\
#$ -N {}\n\
#$ -pe mpi 32\n\
module load python/2.7.12\n\
module unload compilers mpi\n\
module load compilers/intel/2017/update4\n\
module load mpi/intel/2017/update3/intel\n\
module load nwchem/6.8-47-gdf6c956/intel-2017\n\
module list\n\
mpirun -np $NSLOTS -machinefile $TMPDIR/machines nwchem {}\n\
	"
	nwname = pdb2[:-4] + '.nw'
	shname = pdb2[:-4] + '.sh'
	genname = pdb2.split('/')[-1][:-4]
	cmol = Chem.MolFromPDBFile(pdb1, removeHs=False)
	cmol_final = Chem.MolFromPDBFile(pdb2, removeHs=False)
	print cmol
	ccomplex, cguest = get_CB_guest_atomxyz(cmol)
	fincomplex, finguest = get_CB_guest_atomxyz(cmol_final)
	cscript = start_script.format(genname)
	for coord in ccomplex + cguest:
		cscript+= positionline.format(*coord)
	cscript+=endgeom
	for coord in fincomplex+finguest:
		cscript+= positionline.format(*coord)
	cscript += end_script.format(genname)

	with open(nwname, 'wb') as w:
		w.write(cscript)
	with open(shname, 'wb') as w:
		w.write(shstring.format(genname, genname+'.nw',genname+'_LOGOUT'))


def make_gaussian_input_files_for_molgroup(list_of_pdbs):
	"""
	:param list_of_pdbs: take in a list of sdf files located into the very same directory
	:return: produce a com file for each of the pdb files, a paramfile and a sh file to launch gaussian
	"""
	path = '/'.join(list_of_pdbs[0].split('/')[:-1]) + '/'
	genname = 'COMPLEX'
	paramfile = path + '{}paramfile'.format(genname)
	nproc = 16
	memperproc = 2
	shfile = path + '{}.sh'.format(genname)
	comstring = "\
%Chk={0}.chk\n\
%NProcShared={1}\n\
%mem={2}gb\n\
#n wB97XD/3-21G Opt\n\
\n\
{0}\n\
\n\
1 1\n\
"
	shstring = "#!/bin/bash -l\n\
#$ -S /bin/bash\n\
#$ -cwd\n\
#$ -l h_rt=12:0:0\n\
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
time g16 < $g16infile > $g16outfile\n\
	"
	posstring = '{0}         {1: .5f}       {2: .5f}       {3: .5f}\n'
	print paramfile
	with open(paramfile, 'wb') as w: pass
	with open(paramfile, 'ab') as a:
		for i, fname in enumerate(sorted(list_of_pdbs)):
			comname = fname[:-4] + '.com'
			a.write('{0:04d}\t{1}\n'.format(i+1, comname.split('/')[-1]))

			cmol = Chem.MolFromPDBFile(fname, removeHs=False)
			ccomplex, cguest = get_CB_guest_atomxyz(cmol)
			ccomstring = comstring.format(fname.split('/')[-1][:-4], nproc, nproc*memperproc)
			for coord in ccomplex+cguest:
				ccomstring = ccomstring+posstring.format(*coord)
			ccomstring += '\n'

			with open(comname, 'wb') as w:
				print comname
				w.write(ccomstring)

	with open(shfile, 'wb') as w:
		w.write(shstring%(memperproc, genname,nproc,i+1, paramfile.split('/')[-1]))


def get_xyz_list_from_xyz(fxyz):
	"""
	:param fxyz: takes in an xyz file
	:return: a list of list containing the atom types and the atomic positions for that xyz
	"""
	with open(fxyz, 'rb') as r:
		xyzlist = [[y[0], float(y[1]), float(y[2]), float(y[3])] for y in [x.strip().split() for x in r.readlines()[2:]]]
		print xyzlist
	return xyzlist, []

def make_gaussian_input_files_for_xyzgroup(list_of_xyz):
	"""
	:param list_of_pdbs: take in a list of sdf files located into the very same directory
	:return: produce a com file for each of the pdb files, a paramfile and a sh file to launch gaussian
	"""
	path = '/'.join(list_of_xyz[0].split('/')[:-1]) + '/'
	genname = 'prot-unbound'
	paramfile = path + '{}paramfile'.format(genname)
	nproc = 16
	memperproc = 4
	shfile = path + '{}.sh'.format(genname)
	# n wB97XD/3-21G opt=(calcall,tight,ts,cartesian,noeigentest)\n\
	comstring = "\
%Chk={0}.chk\n\
%NProcShared={1}\n\
%mem={2}gb\n\
#n wB97XD/6-31g* Opt\n\
\n\
{0}\n\
\n\
1 1\n\
"
	shstring = "#!/bin/bash -l\n\
#$ -S /bin/bash\n\
#$ -cwd\n\
#$ -l h_rt=6:10:0\n\
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
time g16 < $g16infile > $g16outfile\n\
	"
	posstring = '{0}         {1: .5f}       {2: .5f}       {3: .5f}\n'
	print paramfile
	with open(paramfile, 'wb') as w: pass
	with open(paramfile, 'ab') as a:
		for i, fname in enumerate(sorted(list_of_xyz)[::-1]):
			comname = fname[:-4] + '.com'
			a.write('{0:04d}\t{1}\n'.format(i+1, comname.split('/')[-1]))

			ccomplex, cguest = get_xyz_list_from_xyz(fname)
			ccomstring = comstring.format(fname.split('/')[-1][:-4], nproc, nproc*memperproc)
			for coord in ccomplex+cguest:
				ccomstring = ccomstring+posstring.format(*coord)
			ccomstring += '\n\n'

			with open(comname, 'wb') as w:
				w.write(ccomstring)

	with open(shfile, 'wb') as w:
		w.write(shstring%(memperproc, genname,nproc,i+1, paramfile.split('/')[-1]))



def doc_pdb_in_cb6(cb6fline, listofmols):
	"""
	:param cb6fline: file name showing the pdb file for cb6
	:param listofmols: list of pdb file names for merging with cb6
	:return: the merged version of these loose complexes with added hydrogens
	"""
	cb6mol = Chem.MolFromPDBFile(cb6fline, removeHs=False)
	benzene = Chem.MolFromSmiles('c1ccccc1')
	mxylene = Chem.MolFromMolFile('/home/macenrola/Desktop/0.sdf', sanitize=False, removeHs=False)
	oxylene = Chem.MolFromMolFile('/home/macenrola/Desktop/1.sdf', sanitize=False, removeHs=False)
	pxylene = Chem.MolFromMolFile('/home/macenrola/Desktop/2.sdf', sanitize=False, removeHs=False)
	kind = {'7237':oxylene, '7809': pxylene, '7929':mxylene}
	# mxylene = Chem.MolFromSmiles('Cc1cc(C)ccc1')
	# oxylene = Chem.MolFromSmiles('Cc1c(C)cccc1')
	# pxylene = Chem.MolFromSmiles('Cc1ccc(C)cc1')
	# for i, m in enumerate([mxylene, oxylene, pxylene]):
	# 	m=Chem.AddHs(m)
	# 	AllChem.EmbedMolecule(m)
	# 	Chem.MolToMolFile(m, '/home/macenrola/Desktop/{}.sdf'.format(i))
	sp2 = benzene.GetAtoms()[0].GetHybridization()
	for f in listofmols:
		tmpmol = Chem.MolFromPDBFile(f, sanitize=False)
		num = f.split("/")[-1][:4]
		if num[0]=='C':
			continue
		pattern = kind[num]
		Chem.rdMolAlign.AlignMol(pattern, tmpmol, atomMap=zip(range(pattern.GetNumAtoms()),range(tmpmol.GetNumAtoms())))
		# for i in tmpmol.GetAtoms():
		# 	if i.IsInRing():
		# 		i.SetIsAromatic(True)
		# 		i.SetHybridization(sp2)
		# 		i.SetNumExplicitHs(1)
				# i.SetNoImplicit(True)
		# tmpmol.UpdatePropertyCache()
		# Chem.SanitizeMol(tmpmol)
		# print Chem.MolToSmiles(tmpmol)
		# tmpmol = Chem.AddHs(tmpmol, addCoords=True)

		cmol = Chem.CombineMols(pattern, cb6mol)
		Chem.MolToPDBFile(cmol, f[:-4]+'-complex.pdb')

def reformat_avogadro_file_for_reaxff(fname):
	"""
	:param fname: a lammps data file generated by avogadro
	:return: the same file + reformmated with the atom positions reformatted to match lammps requirements
	"""
	pasatoms = False
	pasbonds = False
	outf = fname+'reformated'
	with open(outf, 'wb') as w:
		with open(fname, 'rb') as r:
			for line in r:
				if 'Bonds' in line:
					pasbonds = True
				if line.strip() == '': continue
				if not pasatoms and not pasbonds:
					w.write(line)
				if pasatoms and not pasbonds:
					l = line.strip().split()
					sx, sy, sz = l[4:7]
					tmp = []
					for s in [sx, sy, sz]:
						if len(str(int(float(s)))) ==1:
							tmp.append(' ')
						else:
							tmp.append('')
					sx, sy, sz = tmp
					print l
					s = '{0:>5} {1:} 0.0 {6}{2:>+2.6f}  {7}{3:>+2.6f}  {8}{4:>+2.6f} {5}\n'.format(l[0], l[2], float(l[4]),float(l[5]), float(l[6]), ' '.join(l[-2:]),
																								 sx, sy, sz)
					print s
					w.write(s)
				if 'Atoms' in line:
					pasatoms = True


def make_min_script_cp2k_for_xyzlist(xyzlist):
	"""
	:param xyzlist: Takes in a list of xyz files
	:return: returns a list of cp2k inputs to run the desired minimization, the charge by default will be zero, modify by hand if needed
	"""
	script = """
 &GLOBAL
   PRINT_LEVEL  LOW
   PROJECT_NAME {0}
   RUN_TYPE  GEO_OPT
 &END GLOBAL
 &MOTION
   &GEO_OPT
     TYPE  MINIMIZATION
     OPTIMIZER  BFGS
     MAX_ITER  2000
     MAX_DR     1.0000000000000000E-03
     MAX_FORCE     1.0000000000000000E-03
     RMS_DR     1.0000000000000000E-03
     RMS_FORCE     1.0000000000000000E-03
     STEP_START_VAL  170
   &END GEO_OPT
 &END MOTION
 &FORCE_EVAL
   METHOD  QS
   &DFT
     CHARGE 0
     BASIS_SET_FILE_NAME ./basis-set-popple-3-21gd
     &SCF
       MAX_SCF  200
       EPS_SCF     1.0000000000000001E-05
       SCF_GUESS  ATOMIC
       &DIAGONALIZATION  T
         ALGORITHM  STANDARD
       &END DIAGONALIZATION
       &MIXING  T
         METHOD  PULAY_MIXING
         ALPHA     5.0000000000000000E-01
         NBUFFER  5
       &END MIXING
       &PRINT
         &RESTART  OFF
         &END RESTART
       &END PRINT
     &END SCF
     &QS
       EPS_DEFAULT     9.9999999999999995E-08
       METHOD  GAPW
     &END QS
     &MGRID
       NGRIDS  4
       CUTOFF     2.0000000000000000E+02
       REL_CUTOFF     3.0000000000000000E+01
     &END MGRID
     &XC
       DENSITY_CUTOFF     1.0000000000000000E-10
       GRADIENT_CUTOFF     1.0000000000000000E-10
       TAU_CUTOFF     1.0000000000000000E-10
       &XC_GRID
         XC_SMOOTH_RHO  NN50
         XC_DERIV  NN50_SMOOTH
       &END XC_GRID
       &XC_FUNCTIONAL  NO_SHORTCUT
         &PBE  T
         &END PBE
       &END XC_FUNCTIONAL
       &VDW_POTENTIAL
         POTENTIAL_TYPE  PAIR_POTENTIAL
         &PAIR_POTENTIAL
           R_CUTOFF     1.6000000000000000E+01
           TYPE  DFTD3
           PARAMETER_FILE_NAME dftd3.dat
           REFERENCE_FUNCTIONAL PBE
         &END PAIR_POTENTIAL
       &END VDW_POTENTIAL
     &END XC
     &POISSON
       POISSON_SOLVER  WAVELET
       PERIODIC  NONE
     &END POISSON
   &END DFT
   &SUBSYS
     &CELL
       A     2.5000000000000014E+01    0.0000000000000000E+00    0.0000000000000000E+00
       B     0.0000000000000000E+00    2.5000000000000014E+01    0.0000000000000000E+00
       C     0.0000000000000000E+00    0.0000000000000000E+00    2.5000000000000014E+01
       PERIODIC  NONE
       MULTIPLE_UNIT_CELL  1 1 1
     &END CELL
     &COORD
{1}
     &END COORD
     &KIND H
       BASIS_SET 3-21Gx
       POTENTIAL ALL
       &POTENTIAL
1 0 0
0.2000000000000000E+00
       &END POTENTIAL
     &END KIND
     &KIND O
       BASIS_SET 3-21Gx
       POTENTIAL ALL
       &POTENTIAL
4 4 0
0.2476208600000000E+00
       &END POTENTIAL
     &END KIND
     &KIND N
       BASIS_SET 3-21Gx
       POTENTIAL ALL
       &POTENTIAL
4 3 0
0.2891792300000000E+00
       &END POTENTIAL
     &END KIND
     &KIND C
       BASIS_SET 3-21Gx
       POTENTIAL ALL
       &POTENTIAL
4 2 0
0.3488304500000000E+00
       &END POTENTIAL
     &END KIND
     &TOPOLOGY
       NUMBER_OF_ATOMS  172
       MULTIPLE_UNIT_CELL  1 1 1
       &CENTER_COORDINATES  T
       &END CENTER_COORDINATES
     &END TOPOLOGY
   &END SUBSYS
 &END FORCE_EVAL"""

	shscript = """
#!/bin/bash -l

#$ -S /bin/bash

#$ -l h_rt=2:10:00

#$ -l mem=3G

#$ -N {0}

#$ -pe mpi 36

#$ -cwd 

# 8. Load required module alongside default module
module unload compilers mpi
module load compilers/gnu/4.9.2
module load mpi/openmpi/1.8.4/gnu-4.9.2
module load openblas/0.2.14/gnu-4.9.2
module load fftw/3.3.4/gnu-4.9.2
module load libxc/2.2.2/gnu-4.9.2
module load libint/1.1.4/gnu-4.9.2
module load quip/18c5440/gnu-4.9.2
module load cp2k/4.1/ompi/gnu-4.9.2

# Gerun is our mpirun wrapper which sets number of cores and machinefile for you
gerun cp2k.popt -in {1}
	"""

	for f in xyzlist:
		outf = f[:-4] + '.inp'
		outsh = f[:-4] + '.sh'
		with open(f, 'rb') as r:
			with open(outf, 'wb') as w:
				lines = r.readlines()
				print lines
				w.write(script.format(f.split('/')[-1][:-4], '           '+'           '.join(lines[2:])))
				with open(outsh, 'wb') as wsh:
					wsh.write(shscript.format(f.split('/')[-1][:-4], f.split('/')[-1][:-4]+'.inp'))


if __name__ == "__main__":
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import numpy as np
	import glob

	# oxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/pXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	# format_mol_for_vasp(oxylene, 'pXylene_protonated')
	# oxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/oXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	# mxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/mXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol', removeHs=False)
	# pxylene = Chem.MolFromMolFile('/home/macenrola/Documents/vasp/xylene/alignments/pXyleneProtonated_wB97XD_631Gd_small_complexes.com_OUT.mol',	removeHs=False)
	# align_xylenes(oxylene, mxylene, pxylene)

	# oxylene = Chem.MolFromPDBFile('/home/macenrola/Documents/nwchem/OXYLENE_OUT.pdb', removeHs=False, sanitize=False)
	# pxylene = Chem.MolFromPDBFile('/home/macenrola/Documents/nwchem/PXYLENE_OUT.pdb', removeHs=False, sanitize=False)
	# mxylene = Chem.MolFromPDBFile('/home/macenrola/Documents/nwchem/MXYLENE_OUT.pdb', removeHs=False, sanitize=False)
	# print oxylene, pxylene, mxylene
	# align_xylenes(oxylene, mxylene, pxylene)


	#
	# flist = sorted(glob.glob('/home/macenrola/Desktop/intermediate/try-stepper-popping-ts/oxylene-cb6-neutral.xyz'))
	# make_gaussian_input_files_for_xyzgroup(flist)
	# print flist


	flist = sorted(glob.glob('/home/macenrola/Desktop/prepareinputs/cp2kneb/raw-G16-Outputs/xyz-oxylene/*.xyz'))
	make_min_script_cp2k_for_xyzlist(flist)


	# flist =glob.glob('/home/macenrola/Documents/XYLENE/popping_guests/nwchem_hf/*0.pdb')
	# for f in flist:
	# 	make_neb_input_for_nwchem(f, f[:-5]+'9.pdb')
	# 	make_neb_input_for_nwchem(f, f[:-4] + '_up.pdb')

	# flist = glob.glob('/home/macenrola/Documents/XYLENE/popping_guests/popping-neutral-xyl/stepper-ts/oxylene-cb6-neutral.sdf')
	# make_nw_paramfile(flist)


	# reformat_avogadro_file_for_reaxff('/home/macenrola/Documents/XYLENE/reaxff-md/prot-xylene-bulk/mxylene/prot-mXylene.lmpdat')



	# doc_pdb_in_cb6('/home/macenrola/Documents/XYLENE/neutral_cb6/CB6.pdb', glob.glob('/home/macenrola/Documents/XYLENE/neutral_cb6/*.pdbqt.pdb'))

