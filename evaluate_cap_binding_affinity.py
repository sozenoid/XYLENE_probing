# -*- coding: utf-8 -*-
def build_mol_from_smiles(SMI=None, pdbfile=None, mol = None):
	"""

	Parameters
	----------
	SMI : TYPE a string representing a molecule, needs to be valid otherwise returns None
		DESCRIPTION.
		 trop is "C1=CC=C[CH+]C=C1"
	Returns a 3D version of the molecule
	-------
	Also produces a pdb file and a mol2 file from antechamber
	"""
	#BUILDS 3D structure
	if SMI is not None:
		mol = Chem.MolFromSmiles(SMI)
		mol = Chem.AddHs(mol)
		AllChem.EmbedMolecule(mol)
		
	elif pdbfile is not None:
		mol = Chem.MolFromPDBFile(pdbfile, removeHs=False)
	elif mol is not None:
		pass 
	else:
		print "No valid input. Provide valid SMI, pdbfile or mol"
		return 
	AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
	#GETS FORMAL CHARGE
	charge = Chem.GetFormalCharge(mol)
	#ASSIGN RANDOM NAME
	rdstring = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(32)])
	#SAVE IN OUTPUTS
	Chem.MolToPDBFile(mol, "outputs/{}.pdb".format(rdstring), flavor =12)	
	
	#BUILDS THE MOL2 with correct charges and residue names for the docking
	make_pdb_with_named_residue("outputs/{}.pdb".format(rdstring), "GST")
	#pdb2mol2(rdstring, charge)
	convert_cmd = "{0} outputs/{1}.pdb -O outputs/{1}.mol2".format(obabel_path, rdstring)
	proc = subprocess.check_output(convert_cmd.split(), shell=False)

	
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
	proc = subprocess.check_output(get_mol2_cmd, shell=False)
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
	proc = subprocess.check_output(frcmod_cmd, shell=False)
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
		w.write("outputs/{}.mol2".format(rdid))
		
	ledock_cmd = "./ledock_linux_x86 ledock.config".format(rdid).split()
	print ledock_cmd
	proc = subprocess.check_output(ledock_cmd, shell=False)

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
	#try:
	with open("outputs/{}.dok".format(rdid), 'rb') as r:
		for line in r:
			if i>n: break
			curcomplex = "outputs/{}-{}.pdb".format(rdid, i)
			with open(curcomplex, "ab") as a:
				#a.write(line)
				if "ATOM" in line:
					pts = line.strip().split()
					# EXAMPLE:  ATOM      2  C1  LIG     0      -6.550  -3.810  -2.641
					rec_line = ("{0: <8}{1: >3}  {2: <4}{3}     {4: >1}      {5: >6}  {6: >6}  {7: >6}\n".format(pts[0], pts[1], pts[2].upper(), pts[3], pts[4], pts[5], pts[6], pts[7]))
				else:
					rec_line = line
				a.write(rec_line)
			if "END" in line:
				make_pdb_with_named_residue(curcomplex, "GST")
				guest_list.append(curcomplex)
				i=i+1
			elif "Cluster" in line:
				with open(summaryfile, 'ab') as a:
					#a.write('{}\t{}\t{}'.format(name, i, line))
					results.append((float(line.split()[-2]), rdid, i))
# =============================================================================
# 	except:
# 		pass
# =============================================================================
	
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
source oldff/leaprc.ff14SB
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
	
	proc = subprocess.check_output(tleap_gu_cmd.split(), shell=False)
	proc = subprocess.check_output(tleap_complex_cmd.split(), shell=False)
	
	return 


def get_apbs_files(rdid):
	"""

	Parameters
	----------
	rdid : string
		this is the random id representing the molecule

	
	Returns the apbs files to run the solvation test
	Prepares also the pqr file that serves as basis to the apbs optimization
	Better to use a mol2 file with correct charges to obtain the PQR file
	Indeed there is no way to supply and distribute the charge after that
	
	THERE NEEDS TO BE A MOL2 file for this rdid 
	
	the file will be written to outputs/rdid-apbs.inp 
	-------
	None.

	"""
	outf = "outputs/{}-apbs.inp".format(rdid)
	mol22pqr_cmd = "{0} outputs/{1}.mol2 -O outputs/{1}.pqr".format(obabel_path, rdid)
	proc = subprocess.check_output(mol22pqr_cmd.split(), shell=False, stderr=subprocess.STDOUT) # the last argument is to suppress the obabel error messages
	with open(apbs_inp, "r") as r:
		apbs_input = r.read()
	with open(outf, "w") as w:
		w.write(apbs_input.format(rdid))
		
	return outf.split("/")[-1] 

def compute_solvation_from_apbs_file(apbsfile):
	"""
	Parameters
	----------
	apbsfile : string
		there exist a file outputs/apbsfile ending in .inp

	Returns
	-------
	Will return the apolar solvation energy and the polar solvation energy in kcal/mol
	"""
	
	apbs_cmd = "{} {}".format(apbs_path, apbsfile)
	os.chdir("outputs")
	apbs_output = subprocess.check_output(apbs_cmd.split(), shell=False)
	os.chdir("..")
	output_file = "outputs/{}.out".format(apbsfile[:-4])
	E_apol, E_pol = 0,0
	with open(output_file, "w") as w:
		w.write(apbs_output)
	for l in apbs_output.split("\n")[::-1]:
		if "Total non-polar energy =" in l:
			E_apol = float(l.split()[-2])/4.18 #returns kcal/mol
		elif "Total electrostatic energy =" in l:
			E_pol = float(l.split()[-2])/4.18 #returns kcal/mol
	return E_pol, E_apol 

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
	a cp2k file that will be written to outputs/rdid-motion.inp
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
	proc = subprocess.check_output(convert_to_pdb_cmd.split(), shell=False)
	return opti_best.split("/")[-1][:-4]


def run_cp2k_file(cp2k_input_file):
	"""
	Parameters
	----------
	cp2k_input_file : string
		DESCRIPTION.

	Returns None: will just run the cp2k_input_file and return the energy value for GEO_OPT and the S_VIB value for VIBRATIONAL_ANALYSIS
	-------
	None.
	"""
	run_cmd = "{} -i {}".format(cp2k_path, cp2k_input_file, cp2k_input_file[:-4])
	os.chdir("outputs")
	result = subprocess.check_output(run_cmd.split(), shell=False)
	outf = cp2k_input_file[:-4]+".out"
	with open(outf, "w") as w:
		w.write(result)
	if "GEO_OPT" in cp2k_input_file:
		for l in result.split("\n")[::-1]:
			if "ENERGY| Total FORCE_EVAL ( FIST ) energy (a.u.):" in l:
				energy = float(l.split()[-1])
				os.chdir("..")
				return energy*627.5 #returns kcal/mol
	elif "VIBRATIONAL_ANALYSIS" in cp2k_input_file:
		for l in result.split("\n")[::-1]:
			if " VIB|              Entropy [kJ/(mol K)]:" in l:
				s_vib = float(l.strip().split()[-1])
				os.chdir("..")
				return s_vib*1000/4.18 # Returns the total entropy as cal/mol
	else:
		os.chdir("..")
		return
	
	 


# =============================================================================
# # NOT NEEDED AS CP2K computes the whole entropy straight away
# def get_rotational_and_translational_entropy(rdid):
# 	"""
# 	Parameters
# 	----------
# 	rdid : string
# 		Mol that we are targetting
# 
# 	Returns
# 	-------
# 	The rotational and translational entropy for that molecule
# 	They are computed and cross checked using formulas and data from page 137 MacQuarrie 1976, Statistical Mechanics.
# 	Values in cal/K for rotational and translational entropy can be found here (https://doi.org/10.1021/j150379a015   J. Phys. Chem. 1937, 41, 1, 149–158) for benzene
# 	"""
# 	pressure = 1e5 # Pascals or 1 bar
# 	kbolz, h, e, pi, T, amu, R = 1.38e-23, 6.63e-34, np.exp(1), 3.141592, 300, 1.66e-27, 8.3145 # all Si units
# 	# make the rdkit mol
# 	mol = Chem.MolFromPDBFile("outputs/{}.pdb".format(rdid), removeHs=False)
# 	# Get the inertia moments in amu*angstrom**2
# 	Ia = Chem.Descriptors3D.PMI1(mol)
# 	Ib = Chem.Descriptors3D.PMI2(mol)
# 	Ic = Chem.Descriptors3D.PMI3(mol)
# 	
# 	#Get the molecular mass
# 	M = Descriptors.ExactMolWt(mol) # in amu
# 	
# 	#Get the translational entropy 
# 	S_trans = np.log((2*pi*(M*amu)*kbolz*T/h**2)**1.5*e**2.5*kbolz*T/pressure) # This is in fact S_trans/Nk or S_trans/R
# 	
# 	#Get the sigma value (first gets an xyz as pymatgen reads xyz files)
# 	xyz_cmd = "{0} {1}/outputs/{2}.pdb -O {1}/outputs/{2}.xyz".format(obabel_path, wdir, rdid)
# 	proc = subprocess.check_output(xyz_cmd.split(), shell=False)
# 
# 	sigma_cmd = "{} {} {}/outputs/{}.xyz ".format(sigma_python_executable, sigma_finding_script, wdir, rdid)
# 	output = subprocess.check_output(sigma_cmd.split(), shell=False)
# 	print output
# 	sigma = int(output.split()[-1])
# 
# 	S_rot = np.log(pi**.5 *e**1.5/sigma*(T**3*(8*pi**2*kbolz/h**2)**3*Ia*Ib*Ic*(amu*1e-20)**3)**.5) # actually S_rot/Nk or S_rot/R
# 	
# 	
# 	#print 1/(Ia*(8*pi**2*kbolz/h**2)*(amu*1e-20)), 1/(Ib*(8*pi**2*kbolz/h**2)*(amu*1e-20)), 1/(Ic*(8*pi**2*kbolz/h**2)*(amu*1e-20)) # provides accurate rotational temperatures for ammonia as taken from MacQuarrie 1976, Statistical Mechanics
# 	
# 	return S_trans*R/4.18, S_rot*R/4.18 # Returns the entropy in cal/K
# =============================================================================

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

def build_the_reference_dictionary(cb = "CB7-GOOD", incl = "adamantanone-GOOD", bi_compl="adamantanone-docked-named"):
	"""
	Returns
	-------
	A dictionary that contains the values of entropy and energy for the binary complex
	upon which the docking is performed 
	
	The building blocs are contained in the docking_targets file and the mol2, pqr, pdb and prmtop files exist 
	Everything is contained in the output file
	"""

	# last_time_ref={'E_apol_adamantanone-docked-named': 4.300049187995215, 'E_pol_CB7-GOOD': 10282.201941279905, 'E_tot_CB7-GOOD': 1582.9666829580815, 'S_tot_adamantanone-GOOD': 78.41123444976077, 'S_tot_CB7-GOOD': 302.7110406698565, 'E_apol_adamantanone-GOOD': 0.5367150016944977, 'E_apol_CB7-GOOD': 4.157477875138756, 'E_tot_adamantanone-docked-named': 1558.8121405502288, 'E_pol_adamantanone-docked-named': 10617.540482114833, 'E_pol_adamantanone-GOOD': 301.1580062090909, 'S_tot_adamantanone-docked-named': 354.3125167464115, 'E_tot_adamantanone-GOOD': 11.125611172566936}
	ref_dic = {}

	for system in [cb, incl, bi_compl]:		
		print "OPTIMISE {}".format(system)
		system_opt_cp2k_file = get_cp2k_file(system, cp2k_opti_file, "GEO_OPT") # produces an rdkit opti file
		system_energy = run_cp2k_file(system_opt_cp2k_file)
		system_opti_pdb = get_cp2k_optimised_mol(system)
	
		
		print "PERFORM THE NMA FOR THE {} AND EXTRACT ENTROPY VALUES".format(system)
		system_vib_cp2k_file = get_cp2k_file(system_opti_pdb, cp2k_opti_file, "VIBRATIONAL_ANALYSIS") # produces an rdkit vibrational file
		S_tot_system = run_cp2k_file(system_vib_cp2k_file)
		#S_trans_complex, S_rot_complex = get_rotational_and_translational_entropy(system_opti_pdb)
	
		
		print "OBTAINS POLAR AND APOLAR CONTRIBUTIONS TO THE SOLVATION FOR {}".format(system)
		apbs_input_file_system = get_apbs_files(system)
		E_pol_system, E_apol_system = compute_solvation_from_apbs_file(apbs_input_file_system)
		
		ref_dic["S_tot_{}".format(system)]  = S_tot_system
		ref_dic["E_tot_{}".format(system)]  = system_energy
		ref_dic["E_apol_{}".format(system)] = E_apol_system
		ref_dic["E_pol_{}".format(system)]  = E_pol_system

	
	return ref_dic
	
def get_ref_dic():
	return {'E_apol_adamantanone-docked-named': 4.3000, 'E_pol_CB7-GOOD': 10282.2019, 'E_tot_CB7-GOOD': 1582.9667, 'S_tot_adamantanone-GOOD': 78.4112, 'S_tot_CB7-GOOD': 302.7110, 'E_apol_adamantanone-GOOD': 0.5367, 'E_apol_CB7-GOOD': 4.1574, 'E_tot_adamantanone-docked-named': 1558.8121, 'E_pol_adamantanone-docked-named': 10617.5404, 'E_pol_adamantanone-GOOD': 301.1580, 'S_tot_adamantanone-docked-named': 354.3125, 'E_tot_adamantanone-GOOD': 11.1256}

def estimate_dG(mol_representation, ref_dic):
	"""
	Parameters
	----------
	SMI : smiles
		smiles as an input for the computation of dG for a cap binding with a binary inclusion complex

	Returns
	-------
	the binding affinity in kcal/mol
	"""
	# print build_the_reference_dictionary().__repr__()
	t0 = time.time()
	print "BUILDING THE MOL, DOCKING IT AND RECONSTRUCTING THE RESULTING PDBs ({}s)".format(time.time()-t0)
	mol, charge, rdstring = build_mol_from_smiles(SMI=mol_representation) # creates the molecule from smiles
	dock_mol_to_host(rdstring, docking_targetPDB) # will dock the molecule
	best_complex, best_guest = dok2pdb(rdstring, charge, docking_targetPDB) # will convert the dok file to a pdb again

	print "GET THE TOPOLOGY FILES FOR THE BEST COMPLEX: {}".format(best_complex)
	get_topology_files(best_guest, best_complex) # produce the topology files for the cap and the ternary complex
	
	print "OPTIMISE THE RESULTING COMPLEX USING CP2K (no charge, already included in the prmtop) ({}s)".format(time.time()-t0)
	complex_opt_cp2k_file = get_cp2k_file(best_complex, cp2k_opti_file, "GEO_OPT") # produces an rdkit opti file
	complex_energy = run_cp2k_file(complex_opt_cp2k_file)
	complex_opti_pdb = get_cp2k_optimised_mol(best_complex)

	
	print "OPTIMISE THE BEST GUEST ({}s)".format(time.time()-t0)
	guest_opt_cp2k_file = get_cp2k_file(best_guest, cp2k_opti_file, "GEO_OPT") # produces an rdkit opti file
	guest_energy = run_cp2k_file(guest_opt_cp2k_file)
	guest_opti_pdb = get_cp2k_optimised_mol(best_guest)	
	
	print "PERFORM THE NMA FOR THE COMPLEX AND EXTRACT ENTROPY VALUES ({}s)".format(time.time()-t0)
	complex_vib_cp2k_file = get_cp2k_file(complex_opti_pdb, cp2k_opti_file, "VIBRATIONAL_ANALYSIS") # produces an rdkit vibrational file
	S_tot_complex = run_cp2k_file(complex_vib_cp2k_file)
	#S_trans_complex, S_rot_complex = get_rotational_and_translational_entropy(complex_opti_pdb)

	print "PERFORM THE NMA FOR THE GUEST AND EXTRACT ENTROPY VALUES ({})".format(time.time()-t0)
	guest_vib_cp2k_file = get_cp2k_file(guest_opti_pdb, cp2k_opti_file, "VIBRATIONAL_ANALYSIS") # produces an rdkit vibrational file
	S_tot_guest = run_cp2k_file(guest_vib_cp2k_file)
	#S_trans_guest, S_rot_guest = get_rotational_and_translational_entropy(guest_opti_pdb)
	
	print "OBTAINS POLAR AND APOLAR CONTRIBUTIONS TO THE SOLVATION FOR COMPLEX ({})".format(time.time()-t0)
	apbs_input_file_complex = get_apbs_files(best_complex)
	E_pol_complex, E_apol_complex = compute_solvation_from_apbs_file(apbs_input_file_complex)
	
	print "OBTAINS POLAR AND APOLAR CONTRIBUTIONS TO THE SOLVATION FOR GUEST ({}s)".format(time.time()-t0)
	apbs_input_file_guest = get_apbs_files(best_guest)
	E_pol_guest, E_apol_guest = compute_solvation_from_apbs_file(apbs_input_file_guest)
	
	print "="*30
	print "CAP"
	cap_summary = "E_tot: {0:4.4f}, S_tot: {1:4.4f}, E_pol: {2:4.4f}, E_apol: {3:4.4f} (all kcal/mol except S_tot cal/mol/K)".format(
		guest_energy,
		S_tot_guest,
		E_pol_guest,
		E_apol_guest
		)
	print cap_summary
	print "TERNARY COMPLEX"
	complex_summary = "E_tot: {0:4.4f}, S_tot: {1:4.4f}, E_pol: {2:4.4f}, E_apol: {3:4.4f}".format(
		complex_energy,
		S_tot_complex,
		E_pol_complex,
		E_apol_complex
		)
	print complex_summary
	print "BINARY COMPLEX (target)"
	binary_summary = "E_tot: {0:4.4f}, S_tot: {1:4.4f}, E_pol: {2:4.4f}, E_apol: {3:4.4f}".format(
		ref_dic["E_tot"],
		ref_dic["S_tot"],
		ref_dic["E_pol"],
		ref_dic["E_apol"]
		)
	print binary_summary
	print "difference"
	difference_summary = "E_tot: {0:4.4f},S_tot: {1:4.4f}, E_pol: {2:4.4f}, E_apol: {3:4.4f} |||| dG= {4:4.4f} kcal/mol".format(
		complex_energy-guest_energy-ref_dic["E_tot"],
		S_tot_complex-S_tot_guest-ref_dic["S_tot"],
		E_pol_complex-E_pol_guest-ref_dic["E_pol"],
		E_apol_complex-E_apol_guest-ref_dic["E_apol"],
		(complex_energy-guest_energy-ref_dic["E_tot"]-298/1000*(S_tot_complex-S_tot_guest-ref_dic["S_tot"])+
	    E_pol_complex-E_pol_guest-ref_dic["E_pol"]+E_apol_complex-E_apol_guest-ref_dic["E_apol"])
		)
	print difference_summary
	print "="*30

	return rdstring, {"CAP":cap_summary, "TER_CMP": complex_summary, "diff":difference_summary}

def split_sdf_file(fname, n):
	"""
	Parameters
	----------
	fname : str
		a fname that contains a number of molecules in sdf file
	n: int
		the number of molecule in each subfile that will be produced  

	Returns
	-------
	files named fname-part-i.sdf

	"""
	supp = Chem.SDMolSupplier(fname, removeHs=False)
	for i, mol in enumerate(supp):
		if i%n==0:
			if i!=0:
				w.close()
			k = i/n+1
			w = Chem.SDWriter("{}-part-{}.sdf".format(fname[:-4], i))
		w.write(mol)
			
	return 

def tar_it_all(rdid):
	"""
	rdid : a molecule id value
	will tar all related files to outputs/rdid-ALL-SAVED.tar
	"""
	print "SAVING AND DELETING {} RELATED FILES".format(rdid)
	os.chdir("outputs")
	to_tar = ' '.join(glob.glob("{}*".format(rdid)))
	tar_cmd = "tar -cvzf {0}-ALL-SAVED.tar {1}".format(rdid, to_tar)
	result = subprocess.check_output(tar_cmd.split(), shell=False)
	rm_cmd = "rm {}".format(to_tar)
	result = subprocess.check_output(rm_cmd.split(), shell=False)
	os.chdir("..")
	return 

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
	from rdkit.Chem import AllChem, Descriptors3D, Descriptors 
	import random, string
	import subprocess, os, glob, sys
	from mol_ops_amber import make_pdb_with_named_residue
	import numpy as np
	import pickle 
	import time
	adamantanoneMOL2, CBMOL2, docking_targetPDB, cp2k_opti_file, apbs_inp = "docking_targets/adamantanone-GOOD.mol2", "docking_targets/CB7-GOOD.mol2", "docking_targets/adamantanone-docked-named.pdb", "opti_vib_cp2k.inp", "apbs_inp"
	#WORKSATION
# =============================================================================
# 	wdir = "/home/macenrola/Documents/ML/ChemTS/new_scoring_for_mcts"
# 	obabel_path = "/home/macenrola/anaconda3/envs/chemts/bin/obabel"
# 	ledock_path = "" # Ledock should be in the same folder as wdir
# 	apbs_path = "/usr/bin/apbs"
# 	antechamber_path = "/home/macenrola/anaconda3/envs/chemts/bin"
# 	cp2k_path = "/home/macenrola/anaconda3/pkgs/cp2k-6.1.0-hc6cf775_3/bin/cp2k.sopt" # double check that cp2k is executing on a single core as it should
# =============================================================================
	
	#MYRIAD
	wdir = "/home/uccahcl/Scratch/FIND_CAP/FOR_SALE"
	obabel_path = "/home/uccahcl/OB232/bin/obabel"
	ledock_path = "" # Ledock should be in the same folder as wdir
	apbs_path = "/home/uccahcl/apbs-pdb2pqr/bin/apbs"
	antechamber_path = "/home/macenrola/anaconda3/envs/chemts/bin"
	cp2k_path = "/home/uccahcl/cp2k/exe/local/cp2k.sopt" # double check that cp2k is executing on a single core as it should
	
	# sigma_finding_script="/home/macenrola/Documents/XYLENE_probing/find_symmetry_number_and_point_group.py"
	# sigma_python_executable = "/home/macenrola/anaconda3/envs/chemvae/bin/python"
	
	binary_complex_values = {"E_tot":1558.8121, "S_tot": 354.3125, "E_pol":10617.5404, 'E_apol':4.3} #all kcal/mol except S_tot in cal/mol/K
	
	print "MOVING TO {}".format(wdir)
	os.chdir(wdir)
	#fname = sys.argv[1]
	fname = "outputs/for_sale_smi_part00"
	tot_dic = {}
	with open(fname, "r") as r:
		for line in r:
			smi, name = line.strip().split()
			rdid, res_dic = estimate_dG(smi, binary_complex_values)
			tot_dic["{}-{}".format(name, rdid)] = res_dic
			tar_it_all(rdid)
	with open("{}.resdic".format(fname[:-4]), "w") as w:
		pickle.dump(tot_dic, w)
	with open("{}.resdic".format(fname[:-4]), "r") as r:
		rec = pickle.load(r)
		for k in rec:
			print k, rec[k] 
# =============================================================================
# 	dok2pdb("sQNLnWfSgIo7QVnuyaHoAWkSscNzqhHQ", 0, docking_targetPDB)
# =============================================================================
		
# =============================================================================
# 	estimate_dG("c1ccccc1", binary_complex_values)
# =============================================================================
