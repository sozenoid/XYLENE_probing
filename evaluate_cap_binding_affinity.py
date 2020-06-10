def build_mol_from_smiles(SMI="c1ccccc1"):
	"""

	Parameters
	----------
	SMI : TYPE a string representing a molecule, needs to be valid otherwise returns None
		DESCRIPTION.

	Returns a 3D version of the molecule
	-------
	Also produces a pdb file and a mol2 file from antechamber
	"""
	#BUILDS 3D structure
	mol = Chem.MolFromSmiles(SMI)
	mol = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol)
	AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
	
	#GETS FORMAL CHARGE
	charge = Chem.GetFormalCharge(mol)
	#ASSIGN RANDOM NAME
	rdstring = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(32)])
	#SAVE IN OUTPUTS
	Chem.MolToPDBFile(mol, "outputs/{}.pdb".format(rdstring), flavor =12)	
	
	#BUILDS THE MOL2 with correct charges and residue names for the docking
	make_pdb_with_named_residue("outputs/{}.pdb".format(rdstring), "GST")
	pdb2mol2(rdstring, charge)
	
	return mol, charge, rdstring

def pdb2mol2(rdid, charge):
	"""
	Parameters
	----------
	rdstring : TYPE a string that represent a molecule for which a pdb is available
		DESCRIPTION. The reference of the molecule
	charge : TYPE integer as string
		DESCRIPTION. provides the charge of the molecule being treated

	Returns
	-------
	None.
	Will produce the corresponding mol2 file for the provided rdid string

	"""
	get_mol2_cmd = "{0}/antechamber -i outputs/{1}.pdb -fi pdb -o outputs/{1}.mol2 -fo mol2 -c bcc -nc {2} -s 2 -j 4".format(antechamber_path, rdid, charge).split()
	proc = subprocess.call(get_mol2_cmd, shell=False)
	return 

def mol22frcmod(rdid):
	"""
	Parameters
	----------
	mol2file : TYPE takes in a mol2 file 
		DESCRIPTION.

	Returns a frcmod file to be used for topology production
	-------
	None.

	"""
	frcmod_cmd = "{0}/parmchk2 -i outputs/{1}.mol2 -f mol2 -o outputs/{1}.frcmod -f frcmod".format(antechamber_path, rdid).split()
	proc = subprocess.call(frcmod_cmd, shell=False)
	return 

def make_docking_script(docking_targetPDB):
	"""

	Parameters
	----------
	docking_targetPDB : a pdb file that is the docking target for the docking procedure
		DESCRIPTION.

	Returns
	-------
	None.
	"""
	docking_script="""
Receptor
{}
RMSD
1.0
Binding pocket
-20 20
-20 20
-20 20
Number of binding poses
20
Ligands list
ligands.list
	""".format(docking_targetPDB)
	
	with open("ledock.config", "w") as w:
		w.write(docking_script)
	return 

def dock_mol_to_host(rdid, pdbtarget):
	"""
	Parameters
	----------
	mol : TYPE the mol id that is assumed to be in wdir/outputs/rdid.mol2
		DESCRIPTION.

	Returns 
	-------
	None.

	"""
	make_docking_script(pdbtarget)
	with open("ligands.list", "w") as w:
		w.write("outputs/{}.mol2".format(rdstring))
		
	ledock_cmd = "./ledock_linux_x86 ledock.config".format(rdstring).split()
	print ledock_cmd
	proc = subprocess.call(ledock_cmd, shell=False)

	return 

def dok2pdb(rdid, charge, reconstructing_pdbtarget, n=0):
	"""
	PRE: TAKES in a rdid with a corresponding mol2 and dok file for the guest
	POST: Will produce the n best complexes that correspond to the docking in pdb and mol2 formats
	"""
	summaryfile = "summary_file"
	with open(summaryfile, 'wb'): pass
	results = []
	
	# ERASE OLD AND SPLITS IN INDIVIDUAL PDB FILES
	for f in glob.glob("outputs/{}-{}.pdb".format(rdid,"*")):
		with open(f, 'wb') as w: pass
	
	i=0
	#SPLITS THE DOK INTO INDIVIDUAL PDBs
	guest_list = []
	try:
		with open("outputs/{}.dok".format(rdid), 'rb') as r:
			for line in r:
				if i>n: break
				curcomplex = "outputs/{}-{}.pdb".format(rdid, i)
				with open(curcomplex, "ab") as a:
					a.write(line)
				if "END" in line:
					make_pdb_with_named_residue(curcomplex, "GST")
					guest_list.append(curcomplex)
					i=i+1
				elif "Cluster" in line:
					with open(summaryfile, 'ab') as a:
						#a.write('{}\t{}\t{}'.format(name, i, line))
						results.append((float(line.split()[-2]), rdid, i))
	except:
		pass
	
	# CREATES THE COMPLEXES
	cb = Chem.MolFromPDBFile(reconstructing_pdbtarget, removeHs=False) #HETATM tags rather than ATOM are needed for safe reading
	complex_list = []
	for f in sorted(glob.glob("outputs/{}-{}.pdb".format(rdid, "*"))):
		print(f)
		try:
			guest = Chem.MolFromPDBFile(f, removeHs=False)
			cbguest = Chem.CombineMols(cb, guest)
			complex_pdb_fname = f.replace(rdid, rdid+'-CB')
			Chem.MolToPDBFile(cbguest, complex_pdb_fname)
			complex_list.append(complex_pdb_fname)
		except:
			print "Complex reassembling failed"
	with open(summaryfile, 'ab') as a:
		for res in sorted(results):
			a.write("{}\t{}\t{} kcal/mol\n".format(res[1], res[2], res[0]))
	complex_list = [x.split("/")[-1][:-4] for x in complex_list]
	guest_list = [x.split("/")[-1][:-4] for x in guest_list]
	for gu in guest_list:
		pdb2mol2(gu, charge)
	return complex_list[0], guest_list[0]

def make_tleap_inputs(rdid, complex_id):
	tleap_gu_file = "outputs/{}.in".format(rdid)
	tleap_complex_file = "outputs/{}.in".format(complex_id)
	prmtop_script="""
source leaprc.gaff
source oldff/leaprc.ff14SB
loadAmberParams {0}.frcmod
GST = loadMol2 {0}.mol2
check GST
loadamberparams gaff.dat
saveamberparm GST {0}.prmtop {0}.rst7
saveoff GST {0}.lib
quit
""".format("outputs/{}".format(rdid))
	with open(tleap_gu_file, "w") as w:
		w.write(prmtop_script)
	combine_script = """
source leaprc.gaff
loadoff {0}.lib
loadoff {1}.lib
loadoff {2}.lib
loadamberparams {0}.frcmod
loadamberparams {1}.frcmod
loadamberparams {2}.frcmod
ADA = loadmol2 {0}.mol2
CB7 = loadmol2 {1}.mol2
GST = loadmol2 {2}.mol2
COMPLEX = loadPDB {3}.pdb
savemol2 COMPLEX {3}.mol2 1
saveamberparm COMPLEX {3}.prmtop {3}.rst7
quit
""".format(adamantanoneMOL2[:-5], CBMOL2[:-5], "outputs/{}".format(rdid), "outputs/{}".format(complex_id))
	
	with open(tleap_complex_file, "w") as w:
		w.write(combine_script)
		
	return tleap_gu_file, tleap_complex_file


def get_topology_files(rdid, best_complex_pdb):
	"""
	Parameters
	----------
	rdid : string
		the rdid of the best guest molecule, there is a corresponding pdb
	best_complex_pdb : string
		the rdid of the best complex made from the best guest molecule, there also is a corresponding pdb

	Returns None, the topology files for the complex and the guest as well are just produced by appending .prmtop to their rdid
	-------
	NOTE: here three molecules are combined, it is important that all residue are named (CB7, GST, ADA, ...) and in the same order 
	and namimg convention in the mol2 and pdb files. If not there will be a crash. I suggest building the binary docking target
	Then name all the atoms accordingly. Then split it into the CB7 and inclusion guest. Then the individual molecules of the target
	should have the correct cartesian coordinates and the correct naming convention. Then to use the complex pdb as a docking target,
	with the topology files built from the split individual molecules files.  
	"""
	mol22frcmod(rdid)
	tleap_gu_file, tleap_complex_file = make_tleap_inputs(rdid, best_complex_pdb)
	tleap_cmd = "{0}/tleap -s -f {1} > {1}.log"
	tleap_gu_cmd, tleap_complex_cmd = tleap_cmd.format(antechamber_path, tleap_gu_file), tleap_cmd.format(antechamber_path, tleap_complex_file)
	
	proc = subprocess.call(tleap_gu_cmd.split(), shell=False)
	proc = subprocess.call(tleap_complex_cmd.split(), shell=False)
	
	return 


def get_apbs_files(mol=None, pdbfile=None, charge=0):
	"""

	Parameters
	----------
	mol : TYPE, optional
		DESCRIPTION. The default is None.
	pdbfile : TYPE, optional
		DESCRIPTION. The default is None.

	Returns the apbs files to run the solvation test
	-------
	None.

	"""
	return 


def get_cp2k_file(rdid, cp2k_padron, motion):
	"""
	Parameters
	----------
	rdid : a string
		there are two files outputs/rdid.pdb and outputs/rdid.prmtop files that exist and are valid
	cp2k_padron : TYPE
		a file to be formatted with the two files above
		it will be named outputs/rdid-motion.inp

	Returns
	-------
	a cp2k file that will be written to outputs/rdid.inp
	"""
	cp2k_input_file = "outputs/{}-{}.inp".format(rdid, motion)
	with open(cp2k_padron, "r") as r:
		cp2k_input = r.read()
	
	if motion == "VIBRATIONAL_ANALYSIS": # then the pdb and the prmtop do not have the same rdid, the prmtop does not need the OPTI tag
		cp2k_input_formatted=cp2k_input.format("{}.pdb".format(rdid),"{}.prmtop".format(rdid[:-5]), rdid, motion)
	else:
		cp2k_input_formatted=cp2k_input.format("{}.pdb".format(rdid),"{}.prmtop".format(rdid), rdid, motion)

	with open(cp2k_input_file, "w") as w:
		w.write(cp2k_input_formatted)
	return cp2k_input_file.split("/")[-1]


def get_cp2k_optimised_mol(rdid):
	"""
	PRE  : a molecule outputs/rdid-pos-1.xyz exists an contains a tightly optimsed structure with the 
	amber parameters provided
	POST : Will return a file PDB containing the optimised version
	
	NOTE: CP2K will output an xyz file or a PDB with broken naming, 
	if openbabel is used to get a pdb using an xyz, obviously the naming will be broken as well
	HOWEVER as long as the prmtop is still intact and the atom number not changed, the new pdb with broken naming
	and its old prmtop will still work, so there is that.
	"""
	opti_traj = "outputs/{}-pos-1.xyz".format(rdid)
	opti_best = "outputs/{}-OPTI.xyz".format(rdid)
	with open(opti_traj, "r") as r:
		lines = r.readlines()
	indices = [i for i, l in enumerate(lines) if "i" in l ]
	with open(opti_best, "w") as w:
		w.writelines(lines[indices[-1]-1:])
	convert_to_pdb_cmd = "{} {} -O {}.pdb".format(obabel_path, opti_best, opti_best[:-4])
	proc = subprocess.call(convert_to_pdb_cmd.split(), shell=False)
	return opti_best.split("/")[-1][:-4]


def run_cp2k_file(cp2k_input_file):
	"""
	Parameters
	----------
	cp2k_input_file : string
		DESCRIPTION.

	Returns None: will just run the cp2k_input_file
	-------
	None.
	"""
	run_cmd = "{} -i {}".format(cp2k_path, cp2k_input_file, cp2k_input_file[:-4])
	os.chdir("outputs")
	proc = subprocess.call(run_cmd.split(), shell=False)
	os.chdir("..")
	return 


def merge_2mols(pdb1, pdb2):
	"""
	Parameters
	----------
	pdb1 : TYPE a pdb file 
		DESCRIPTION.
	pdb2 : TYPE a second pdb file
		DESCRIPTION.

	Returns
	-------
	None. will produce a third pdb file having merged the two others
	"""
	mol1 = Chem.MolFromPDBFile(pdb1, removeHs=False)
	mol2 = Chem.MolFromPDBFile(pdb2, removeHs=False)
	mol3 = Chem.CombineMols(mol1, mol2)
	Chem.MolToPDBFile(mol3, pdb2[:-4]+"-docked-named.pdb", flavor=28)
	return guest_list

if __name__=="__main__":
	"""
	ACCESS to the RDKIT, ledock, antechamber, apbs and cp2k are required. Access to a docked 
	The workflow goes as follows.
	From a molecule represented as smiles, a 3D configuration is created with the RDKIT
	The 3D molecule is docked using ledock
	antechamber is used to create topology files 
	
	additional lines
	/home/macenrola/anaconda3/envs/chemts/bin/parmchk2 -i $guest.mol2 -f mol2 -o $guest.frcmod -f frcmod

	additional notes, the names of the atoms like "C14" need to be strictly conserved in the pdb, mol2 and complex pdb
	ideally you should have all the pdbs in the correct position, correct labelling and then get the mol2, frcmod, lib, prmtops and finally call all for the complex
	"""
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import random, string
	import subprocess, os, glob
	from mol_ops_amber import make_pdb_with_named_residue
	
	wdir = "/home/macenrola/Documents/ML/ChemTS/new_scoring_for_mcts"
	adamantanoneMOL2, CBMOL2, docking_targetPDB, cp2k_opti_file = "docking_targets/adamantanone-GOOD.mol2", "docking_targets/CB7-GOOD.mol2", "docking_targets/adamantanone-docked-named.pdb", "opti_vib_cp2k.inp"
	obabel_path = "/home/macenrola/anaconda3/envs/chemts/bin/obabel"
	ledock_path = "" # Ledock should be in the same folder as wdir
	apbs_path = ""
	antechamber_path = "/home/macenrola/anaconda3/envs/chemts/bin/"
	cp2k_path = "/home/macenrola/anaconda3/pkgs/cp2k-6.1.0-hc6cf775_3/bin/cp2k.sopt" # double check that cp2k is executing on a single core as it should
	
	os.chdir(wdir)
	
	mol, charge, rdstring = build_mol_from_smiles() # creates the molecule from smiles
	dock_mol_to_host(rdstring, docking_targetPDB) # will dock the molecule
	best_complex, best_guest = dok2pdb(rdstring, charge, docking_targetPDB) # will convert the dok file to a pdb again
	print best_complex
	get_topology_files(best_guest, best_complex) # produce the topology files for the cap and the ternary complex
	opt_cp2k_file = get_cp2k_file(best_complex, cp2k_opti_file, "GEO_OPT") # produces an rdkit opti file
	run_cp2k_file(opt_cp2k_file)
	opti_pdb = get_cp2k_optimised_mol(best_complex)
	vib_cp2k_file = get_cp2k_file(opti_pdb, cp2k_opti_file, "VIBRATIONAL_ANALYSIS") # produces an rdkit vibrational file
	run_cp2k_file(vib_cp2k_file)

