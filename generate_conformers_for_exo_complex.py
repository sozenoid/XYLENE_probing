def generate_conformers_for_exo_complex(sdfile, nconf=50):
	"""
	PRE: Takes in a sdf file
	POST: Will write nconfs conformers for that exo complex
	"""
	base_complex = Chem.MolFromMol2File(sdfile, removeHs=False)

	frags = Chem.GetMolFrags(base_complex, asMols=True)

	if frags[0].GetNumAtoms()<frags[1].GetNumAtoms():
		guest, cb = frags[0], frags[1]
	else:
		guest, cb = frags[1], frags[0]

	print Chem.MolToMolBlock(guest)

	cb.AddConformer(cb.GetConformer(),1)

	tmp_guest=Chem.Mol(guest)
	AllChem.EmbedMultipleConfs(guest, nconf)
	AllChem.MMFFOptimizeMoleculeConfs(guest)

	for i in range(guest.GetNumConformers()):

		ff = AllChem.MMFFGetMoleculeForceField(guest, pyMMFFMolProperties =AllChem.MMFFGetMoleculeProperties(guest),ignoreInterfragInteractions=False, confId=i)
		print ff.CalcEnergy()
		guest.SetProp("_Name", "{}th-conformer-energy-is-({})kcal/mol".format(i, ff.CalcEnergy()))
		Chem.MolToMolFile(guest, sdfile+"_{}.sdf".format(i), confId=i)


		tmp_mol=Chem.Mol(tmp_guest)
		tmp_mol.GetConformers()
		tmp_mol.AddConformer(guest.GetConformer(id=i),1)
		Chem.rdMolAlign.AlignMolConformers(tmp_mol, atomIds=[6,7,8])

		tmp_comp = Chem.CombineMols(tmp_mol, cb)

		# AllChem.MMFFOptimizeMolecule(tmp_comp, ignoreInterfragInteractions=False)

		Chem.MolToMolFile(tmp_comp, sdfile+"_COMP_{}.sdf".format(i), confId=1)
		tmp_comp = Chem.MolFromMolFile(sdfile+"_COMP_{}.sdf".format(i), removeHs=False)
		print AllChem.MMFFOptimizeMolecule(tmp_comp, ignoreInterfragInteractions=False, maxIters=10000)

		ff = AllChem.MMFFGetMoleculeForceField(tmp_comp, pyMMFFMolProperties =AllChem.MMFFGetMoleculeProperties(tmp_comp),ignoreInterfragInteractions=False)
		print ff.CalcEnergy()
		tmp_comp.SetProp("_Name", "{}th-conformer-energy-is-({})kcal/mol".format(i, ff.CalcEnergy()))
		Chem.MolToMolFile(tmp_comp, sdfile+"_COMP_{}.sdf".format(i))

		tmp_mol=None
		# for k in range(tmp_mol.GetNumConformers()):
		# 	print Chem.MolToMolBlock(tmp_mol, confId=k)
		# print tmp_mol.GetNumConformers()
		# print Chem.MolToMolBlock(tmp_mol)
		# print Chem.MolToMolBlock(guest, confId=i)


def rename_mols(file_list):
	"""
	PRE: Takes in a file list
	POST: Will rename the molecules inside
	"""
	mol_list=[]
	comp_list=[]
	for f in file_list:
		mol = Chem.MolFromMolFile(f, removeHs=False)
		e = mol.GetProp("_Name")
		if "COMP" in f:
			comp_list.append((float(e.split("(")[-1].split(")")[0]), mol))
		else:
			mol_list.append((float(e.split("(")[-1].split(")")[0]), mol))
	s_comp_list =sorted(comp_list)
	s_mol_list =sorted(mol_list)

	for i in range(len(mol_list)):
		s_mol_list[i][1].SetProp("_Name", "Produced by HL from TCL's group, conformer number {}; energy is ({}) kcal/mol".format(i, s_mol_list[i][0]))
		s_comp_list[i][1].SetProp("_Name", "Produced by HL from TCL's group, conformer number {}; energy is ({}) kcal/mol".format(i, s_comp_list[i][0]))
		Chem.MolToMolFile(s_mol_list[i][1],  "/home/macenrola/Documents/heptylamine/HEPTYLAMINE_GUEST_CONF_{}.sdf".format(i))
		Chem.MolToMolFile(s_comp_list[i][1], "/home/macenrola/Documents/heptylamine/HEPTYLAMINE_COMPLEX_CONF_{}.sdf".format(i))
		# mol.SetProp("_Name", "Produced by HL from TCL group {}".format(e))
		# Chem.MolToMolFile(mol, f)

if __name__=="__main__":
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import glob

	rename_mols(glob.glob("/home/macenrola/Documents/heptylamine/CB6_C7NH3+_MMFF94.mol2*.sdf"))
	# generate_conformers_for_exo_complex("/home/macenrola/Documents/heptylamine/CB6_C7NH3+_MMFF94.mol2", 500)
