def interpolate_poscars_keep_selective_flags(poscar1, poscar2, nimages, path=None):
	"""
	:param poscar1 and poscar2: Take poscar1 and poscar2 that correspond to the reagent and product of a chemical reaction
					Both poscar have:
										The same cell parameters
										The same number, type and atom orderning
										The same keywords cart selective
										Consistent flags T T T of F F F or T F T for respective atoms, only the flags, vectors and keywords from poscar1 will be used though (overwriting poscar2 in the last file even)
	:nimages: The number of images, poscar1 and poscar2 are not included
	:return: a POSCAR file in each of the files; for example for 3 images; there will be 00 01 02 03 04 files with poscar1 in 00 and poscar2 in 04; the files 01 02 03 containt the POSCARs of the images
	"""

	def build_poscar_dic_representation(poscar):
		"""
		:param poscar: take in a poscar file
		:return: a dictionary the keys are name, scale, v1, v2, v3, atmnum, keywords, and numbers for the atomic positions
		"""
		posdic = {}
		with open(poscar, 'rb') as r:
			lines = [x.strip().split() for x in r.readlines()]
		print lines
		posdic['name'] = lines[0][0]
		posdic['scale'] = float(lines[1][0])
		posdic['v1'] = [float(x) for x in lines[2]]
		posdic['v2'] = [float(x) for x in lines[3]]
		posdic['v3'] = [float(x) for x in lines[4]]
		posdic['atmnum'] = [int(x) for x in lines[5]]
		posdic['keywords'] = lines[6]
		for i, line in enumerate(lines[7:]):
			posdic[i] = [float(x) for x in line[:3]] + line[3:]
		recenter_dic_coords(posdic)
		return posdic

	def build_poscar_from_dic_repr(dic):
		"""
		:param dic: from a dic  produced by build_poscar_dic_representation
		:return: a string that correspond to a poscar
		"""
		print dic
		name = "{0!s}\n"
		lattice_constant = "{0:.5f}\n"
		vector = '     {0:.6f}    {1:.6f}    {2:.6f}\n'
		number_of_atoms = "     {0:d}"
		typeprocess = "{0!s} "
		positionline = "  {0: .6f}  {1: .6f}  {2: .6f}   {3!s} {4!s} {5!s}\n"
		postring = ''.join([
			name.format(dic['name']),
			lattice_constant.format(dic['scale']),
			vector.format(*tuple(dic['v1'])),
			vector.format(*tuple(dic['v2'])),
			vector.format(*tuple(dic['v3'])),
			''.join([number_of_atoms.format(x) for x in dic['atmnum']])+'\n',
			''.join([typeprocess.format(x) for x in dic['keywords']]) + '\n'
		])
		for i in range(len(dic.keys())-7):
			postring = postring + positionline.format(*tuple(dic[i]))
		return postring

	def recenter_dic_coords(dic):
		"""
		:param dic: Takes in a dic as provided by the function  build_poscar_dic_representation
		:return: will modify its cartesian coordinates so that its center appears at the middle of the cell (otherwise visualization becomes very confusing)
		"""
		x_center_ref, y_center_ref, z_center_ref = dic['v1'][0]*dic['scale']/2., dic['v2'][1]*dic['scale']/2., dic['v3'][2]*dic['scale']/2.
		x_center, y_center, z_center = [],[],[]
		for i in range(len(dic.keys())-7):
			cx,cy,cz = dic[i][:3]
			x_center.append(cx)
			y_center.append(cy)
			z_center.append(cz)

		xc = sum(x_center)*1.0/len(x_center)
		yc = sum(y_center) * 1.0 / len(y_center)
		zc = sum(z_center) * 1.0 / len(z_center)

		print xc, yc, zc

		for i in range(len(dic.keys())-7):
			tc = dic[i][:3]
			dic[i][:3] = [x[0]+x[1]-x[2] for x in zip(tc, [x_center_ref, y_center_ref, z_center_ref],[xc, yc, zc])]

	def build_interpolated_dics(posdic1, posdic2, nimages):
		"""
		:param posdic1 posdic2: two dics corresponding to poscar1 and poscar2 as returned by build_poscar_dic_representation
		:param nimages: the number of images
		:return: a list of dics corresponding to the images poscars, there will be nimages+2 dics with the first one corresponding to poscar1 and the last one corresponding to poscar2
		"""
		import numpy as np
		diclist= []
		for k in range(nimages+2):# Creates the diclist
			if k==0:
				diclist.append(posdic1)
			# elif k==nimages+1:
			# 	diclist.append(posdic2)
			# 	pass
			else:
				diclist.append({'name':posdic1['name']+'->'+posdic2['name'],
								'scale':posdic1['scale'],
								'v1':posdic1['v1'],
								'v2': posdic1['v2'],
								'v3': posdic1['v3'],
								'atmnum':posdic1['atmnum'],
								'keywords':posdic1['keywords']
								})

		for atom_line in range(len(posdic1.keys())-7): # will iterate over the atom numbers (there are 6 fields that are not line numbers)
			i = atom_line
			tmp_pos1 = posdic1[i]
			tmp_pos2 = posdic2[i]
			flags1 = tmp_pos1[3:]
			x_range = np.linspace(tmp_pos1[0], tmp_pos2[0], nimages+2)
			y_range = np.linspace(tmp_pos1[1], tmp_pos2[1], nimages+2)
			z_range = np.linspace(tmp_pos1[2], tmp_pos2[2], nimages+2)
			for j, intermediate_point_coords in enumerate(zip(x_range, y_range, z_range)):
				if j==0:
					continue
				elif j==nimages+1:
					pass
				diclist[j][i] = list(intermediate_point_coords) + flags1
		return diclist


	posdic1 = build_poscar_dic_representation(poscar1)
	posdic2 = build_poscar_dic_representation(poscar2)

	interpolated_dics = build_interpolated_dics(posdic1, posdic2, nimages)

	# if path==None:
	# 	import os
	# 	path = os.getcwd()
	#
	for i, dic in enumerate(interpolated_dics):
		s = build_poscar_from_dic_repr(dic)
		with open(path+'0{}/POSCAR'.format(i), 'wb') as w:
			w.write(s)


def recenter_direct_poscar(poscar):
	"""
	:param poscar: Takes in a direct poscar and recenters it
	:return: another file named poscar_recentered that will contain the recentered data
	"""
	with open(poscar, 'rb') as r:
		with open(poscar+'_recentered', 'wb') as w:
			for i, line in enumerate(r):
				if i<8 or line==' \n':
					w.write(line)
				else:
					coords = [float(x) for x in line.strip().split()]
					coords_center = [x if abs(x)<0.5 else x-1.0 for x in coords]
					w.write('  {0: .16f}  {1: .16f}  {2: .16f}\n'.format(*coords_center))



def plot_interp_neb_traj(fin='/home/macenrola/Documents/nwchem/NEB_MP/final_E_MP_ISOLATED'):
	"""
	:param fin: takes in a final neb path with only the energies
	:return: will plot an image, its interpolation and the value of the extremas
	"""
	from matplotlib import pyplot as plt
	from scipy import interpolate
	import numpy as np
	coords = []
	Es = []
	with open(fin, 'rb') as r:
		for line in r:
			if line[0] == '#':
				continue
			else:
				cur_coords, cur_Es = line.strip().split()
				coords.append(float(cur_coords))
				Es.append(float(cur_Es))
	print coords
	print Es
	min_E = min(Es)
	print min_E
	Es = [(x-min_E)*27.2116*23.06 for x in Es]
	print Es
	left = Es[0]
	right = Es[-1]
	f = interpolate.interp1d(coords, Es, 'cubic')
	x_new = np.linspace(0, 1, 200)
	y_new = f(x_new)
	max_E = max(y_new)
	plt.plot(coords, Es, 'o', x_new, y_new, '-')
	plt.xlabel("Reaction Coordinate")
	plt.ylabel("Energy (kcal/mol)")
	plt.title('M->O: left={0:.4f} kcal/mol; right={1:.4f} kcal/mol; right-left:{2:.4f}kcal/mol; max-left:{3:.4f} kcal/mol'.format(left, right, right-left, max_E-left))
	plt.show()


def pop_xylenes_out_of_CB(rdkitcomplex, lengthA=10, nimages=2):
	"""
	:param rdkitcomplex: given a complex of cb + a molecule, lengthA is the distance in Angstrom that is to be applied to the guest during offsetting
	:nimages: is the number of intermeediate steps, minimum 2
	:return: two additional molecules, one where the guest is popped upward and another where the guest is popped downward, the direction of popping is according to the least principal
			axis of the CB
	"""
	import numpy as np
	from sklearn.decomposition import PCA
	cb_pt, guest_pt = get_CB_guest_atomxyz(rdkitcomplex)

	pca = PCA(n_components=3)
	pca.fit([x[1:] for x in cb_pt])

	complex, guest = Chem.GetMolFrags(rdkitcomplex, asMols=True)
	if complex.GetNumAtoms() < guest.GetNumAtoms():
		complex, guest = guest, complex

	# dummy_atom = Chem.MolFromSmarts('[#2]')
	# AllChem.EmbedMolecule(dummy_atom)
	# print Chem.MolToMolBlock(dummy_atom)

	imgs = []
	for k in range(0, nimages):
		offpt = rdkit.Geometry.rdGeometry.Point3D()
		offvect = pca.components_[2]*lengthA*float(k)/(nimages-1)
		offpt.x, offpt.y, offpt.z = offvect

		offdummypt = rdkit.Geometry.rdGeometry.Point3D()
		offdummyvect = -pca.components_[2]*lengthA*2
		offdummypt.x, offdummypt.y, offdummypt.z = offdummyvect

		tmp_comp = Chem.CombineMols(guest, complex, offset=offpt)
		# tmp_comp_and_dummy = Chem.CombineMols( tmp_comp,dummy_atom, offset=offdummypt) # uncomment to activate dummy atom
		# imgs.append(tmp_comp_and_dummy)
		imgs.append(tmp_comp)
	# offptup = rdkit.Geometry.rdGeometry.Point3D()
	# offptdown = rdkit.Geometry.rdGeometry.Point3D()
	# offvectup = pca.components_[2]*lengthA
	# offvectdown = -pca.components_[2]*lengthA
	# offptup.x, offptup.y, offptup.z = offvectdown
	# offptdown.x, offptdown.y, offptdown.z = offvectup
	#
	# offptmidup = rdkit.Geometry.rdGeometry.Point3D()
	# offptmiddown = rdkit.Geometry.rdGeometry.Point3D()
	# offvectmidup = pca.components_[2]*4
	# offvectmiddown = -pca.components_[2]*4
	# offptmidup.x, offptmidup.y, offptmidup.z = offvectmiddown
	# offptmiddown.x, offptmiddown.y, offptmiddown.z = offvectmidup

	#
	# molup = Chem.CombineMols( guest,complex, offset = offptmidup)
	# moldown = Chem.CombineMols( guest,complex, offset = offptmiddown)
	#
	# molmidup = Chem.CombineMols(guest, complex, offset=offptmidup)
	# molmiddown = Chem.CombineMols(guest, complex, offset=offptmiddown)

	return imgs



def invert_xz_xyz(xyzfile):
	"""
	:param xyzfile: takes in an xyz file
	:return: will return the same xyz file with inversion of the x and z components
	IDEALLY SHOULD OPTIMIZE A SPLIT COMPLEX WITH CONSTRAINTS BECAUSE RVDW=12A
	"""
	with open(xyzfile, 'rb') as r:
		lines = r.readlines()
		for line in lines[2:]:
			line = line.strip().split()
			s = '{}        {}        {}        {}'.format(line[0], line[2], line[3], 10+float(line[1]))
			print s

def print_just_xyz(xyzfile):
	"""
	:param xyzfile: takes in a xyz file
	:return:  will print out only the xyz coordinates
	"""
	with open(xyzfile, 'rb') as r:
		lines = r.readlines()
		for f in lines[2:]:
			print '{}     {}      {}'.format(*f.strip().split()[1:])

def align_molecules(molA, molB):
	"""
	:param mola: molecule a as sdf
	:param molb: molecule b as sdf
	:return: will attempt to align them and return a third molecule aligned
	"""
	core = '[#7]1[#6][#7][#6][#7][#6][#7][#6]1'
	core = Chem.MolFromSmarts(core, True)


	mola = Chem.MolFromMolFile(molA, removeHs=False)
	molb = Chem.MolFromMolFile(molB, removeHs=False)
	cb = Chem.GetMolFrags(mola, asMols=True)[0]
	print Chem.MolToSmiles(cb)

	print Chem.MolToMolBlock(mola)
	print Chem.MolToMolBlock(molb)
	match1 = mola.GetSubstructMatch(core)
	match2 = molb.GetSubstructMatch(core)

	print match2, match1
	AllChem.AlignMol(mola, molb, atomMap=zip(match2, match1), maxIters=1000, reflect=False)  # <- m2 is aligned to m1, return value is the RMSD for the alignment

	Chem.MolToMolFile(mola, molA+'-aligned.sdf')



if __name__ == "__main__":
	from xylene_formatting import get_CB_guest_atomxyz
	import sys
	import rdkit
	from rdkit import Chem
	import glob
	from rdkit.Chem import AllChem
	# interpolate_poscars_keep_selective_flags(sys.argv[1], sys.argv[2], sys.argv[3])
	# interpolate_poscars_keep_selective_flags('/home/macenrola/Documents/vasp/mXylene-Protonated', '/home/macenrola/Documents/vasp/oXylene-Protonated',
	# 										 8, '/home/macenrola/Documents/vasp/NEB_moXylene/')
	# recenter_direct_poscar('/home/macenrola/Documents/vasp/mXylene-Protonated')
	# plot_interp_neb_traj('/home/macenrola/Documents/nwchem/NEB_MO/final_E_MO_ISOLATED')
	# flist = glob.glob('/home/macenrola/Documents/XYLENE/charged_cb7/converged/h-on-cb/*.pdb')
	# for f in flist:
	# 	mol = Chem.MolFromPDBFile(f, removeHs=False)
	# 	imgs = pop_xylenes_out_of_CB(mol)
	# 	for i, k in enumerate(imgs):
	# 		Chem.MolToPDBFile(k, f[:-4]+'_{}.pdb'.format(i))

	# invert_xz_xyz('/home/macenrola/Desktop/prepareinputs/cp2kneb/raw-G16-Outputs/oxylene-transitions/neutral-o-out-neutral-cb/neutral-oxylene.xyz')


	# align_molecules('/home/macenrola/Desktop/prepareinputs/cp2kneb/raw-G16-Outputs/oxylene-transitions/proton-exchange-cb-oxylene/prot-oxylene-cb7.xyz.sdf',
	# 				'/home/macenrola/Desktop/prepareinputs/cp2kneb/raw-G16-Outputs/oxylene-transitions/proton-exchange-cb-oxylene/h-single-oxylene-cb7.xyz.sdf')


	# align_molecules('/home/macenrola/Desktop/prepareinputs/cp2kneb/raw-G16-Outputs/oxylene-transitions/proton-exchange-cb-oxylene/prot-oxylene-cb6.xyz.sdf',
	# 				'/home/macenrola/Desktop/prepareinputs/cp2kneb/raw-G16-Outputs/oxylene-transitions/proton-exchange-cb-oxylene/h-single-oxylene-cb6.xyz.sdf')

	print_just_xyz('/home/macenrola/Desktop/prepareinputs/cp2kneb/raw-G16-Outputs/oxylene-transitions/proton-exchange-cb-oxylene/prot-oxylene-cb6-aligned.xyz')
		# Chem.MolToPDBFile(up, f[:-4]+'_up.pdb')
		# Chem.MolToPDBFile(down, f[:-4]+'_down.pdb')
		# Chem.MolToPDBFile(midup, f[:-4]+'_midup.pdb')
		# Chem.MolToPDBFile(middown, f[:-4]+'_middown.pdb')
