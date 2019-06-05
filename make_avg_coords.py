
################
#
# This part deals with averaging coordinates using z matrices, I don't think
# that it works well for molecules with multiple intricated cycles like CBs
#
#
################

def parse_z_matrix(z_matrix):
	"""
	PRE: Takes in a z matrix as produced by obabel as a text block
	POST: returns a dictionary representation of the matrix
		one key holds the structure ("struct") the other keys are named after the variable names that are expected to be consistent

	"""
	if len(z_matrix.split())<3:
		return None

	zmatdic={}
	lines=z_matrix.split('\n')
	pastvariables=False
	for l in lines:
		if "Variables:" in l:
			pastvariables=True
		if pastvariables and '=' in l:
			temp_line=l.strip().split('=')
			zmatdic[temp_line[0]]=temp_line[1].strip()
		if 'structure' not in zmatdic and not pastvariables:
			zmatdic['structure']=[l]
		elif not pastvariables:
			zmatdic['structure'].append(l)
	zmatdic["structure"]="\n".join(zmatdic["structure"])
	return zmatdic

def average_z_matrix(list_of_z_matrix_as_dic):
	"""
	PRE: Takes a list of dics representing z matrices
	POST: Returns a single dic with all values averaged and another dic with the deviations
	"""
	example_dic=dict.fromkeys(list_of_z_matrix_as_dic[0].keys())
	clustered_dic={}
	avg_clustered_dics=[]
	##### SPLIT THE ORIGINAL LIST BY CLUSTER NUMBER
	for dic in list_of_z_matrix_as_dic:
		if dic["cluster"] not in clustered_dic:
			clustered_dic[dic["cluster"]] = [dic]
		else:
			clustered_dic[dic["cluster"]].append(dic)

	#### COMPUTE THE AVERAGE Z MATRIX CLUSTER WISE
	print [len(x) for x in clustered_dic.values()]
	for c in clustered_dic:
		res_dic = dict.fromkeys(clustered_dic[c][0].keys())
		res_dic["structure"] = clustered_dic[c][0]["structure"]

		for k in res_dic:
			temp_value=[]
			if k=='structure': continue
			for dic in clustered_dic[c]:
				temp_value.append(float(dic[k]))

			# print k,temp_value, statistics.stdev(temp_value)
			# temp_value=[x if x<180 else abs(x-360) for x in temp_value]
			# print k,temp_value, statistics.stdev(temp_value)

			res_dic[k]=1.0*sum(temp_value)/len(temp_value)
			# if res_dic[k]<0: res_dic[k]=res_dic[k]+360
		avg_clustered_dics.append(res_dic)
	return avg_clustered_dics

def write_z_matrix_from_dic(z_matrix_dic, fname):
	"""
	PRE: takes in a z matrix written as dic
	POST: will produce a file formatted as obabel z matrix
	"""
	with open(fname, 'wb') as w:
		w.write(z_matrix_dic["structure"]+"\n")
		w.write("Variables:\n")
		for k in sorted(z_matrix_dic.keys()):
			if k =="structure":
				continue
			w.write("{}={}\n".format(k, z_matrix_dic[k]))

def cluster_list_of_zdics(list_of_z_matrix_as_dic):
	"""
	PRE: Takes in a list of z matrices as cluster_list_of_zdics
	POST: will use PCA or kmeans to cluster the coordinates.
	"""
	full_list=[]
	for dic in list_of_z_matrix_as_dic:
		tmplist=[]
		for k in sorted(dic.keys()):
			if k=="structure":
				continue
			tmplist.append(dic[k])
		full_list.append(tmplist)
		# print tmplist
	full_list_array= np.array([np.array(x) for x in full_list])
	# print(full_list_array.shape)

	######## COMPUTE THE PCA
	# pca=PCA(n_components=2)
	# transformed_coords =  pca.fit_transform(full_list_array)
	# transformed_coords_pd = (pd.DataFrame(transformed_coords))
	# transformed_coords_pd.columns = ["PC1" ,"PC2"]

	####### PLOT AS A SANITY CHECK
	# ax = transformed_coords_pd.plot(kind='scatter', x='PC2', y='PC1', figsize=(16,8))
	# plt.show()
	# print transformed_coords
	# print(pca.explained_variance_ratio_)
	# plt.scatter(transformed_coords.T[0], transformed_coords.T[1])
	# plt.show()

	#### COMPUTE THE KMEANS AND BUILD THE ESTIMATOR
	kmeans = KMeans(n_clusters=100)
	clusters = kmeans.fit_predict(full_list_array)

	#### PLOT AS A SANITY CHECK
	# transformed_coords_pd['cluster'] = pd.Series(clusters.labels_, index=transformed_coords_pd.index)
	# axk =transformed_coords_pd.plot(
    # kind='scatter',
    # x='PC2',y='PC1',
    # c=transformed_coords_pd.cluster.astype(np.float),
    # figsize=(16,8))
	# plt.show()

	####### ADD A KEY TO EACH DICTIONARY INDICATING WHICH CLUSTER
	for i, dic in enumerate(list_of_z_matrix_as_dic):
		dic["cluster"] = clusters[i]
	return list_of_z_matrix_as_dic

def process_z_matrix_trajectory(zmatrix_traj_file):
	"""
	PRE : takes in a zmatrix_traj_file
	POST: will extract all the indivudual blocks and process them using parse_z_matrix
	"""
	temp_z=[]
	zmatlist=[]
	current_snapshot=0
	with open(zmatrix_traj_file, 'rb') as r:
		for i,lines in enumerate(r):
			# if i==100000: break
			if '#' in lines:
				current_snapshot+=1
				if current_snapshot%100==0:
					print "current snapshot is {}".format(current_snapshot)
				if temp_z!=[]:
					temp_z_mat=parse_z_matrix(''.join(temp_z))
					if temp_z_mat is not None:
						zmatlist.append(temp_z_mat)
				temp_z=[]
				temp_z.append(lines)
			else:
				temp_z.append(lines)

		if temp_z!=[]:
			# temp_z.append(lines)
			temp_z_mat=parse_z_matrix(''.join(temp_z))
			if temp_z_mat != []:
				if zmatlist!=[]:
					nkeys = len(zmatlist[-1].keys())
					if nkeys!=len(temp_z_mat.keys()):
						print "corrupted file, expected {} keys but found {}.\n Total indivudual snapshots found {}".format(nkeys, len(temp_z_mat.keys()), len(zmatlist))
					else:
						zmatlist.append(temp_z_mat)
			temp_z=[]

	clustered_zmatlist = cluster_list_of_zdics(zmatlist)
	########## PRINTING THOSE WINNERS
	avg_z_matlist = average_z_matrix(clustered_zmatlist) # will contain as many matrices as there are clusters
	print len(avg_z_matlist)
	for i, avg_z_mat in enumerate(avg_z_matlist):
		write_z_matrix_from_dic(avg_z_mat, "avg_resdic-{}.gzmat".format(i))

	return

#################
#
# This part is averaging an rms aligned trajectory in xyz coordinates
# It may just work better
#
#
##################
def parse_xyz_block(xyz_block):
	"""
	PRE: Takes in an xyz block as a list of lines
	POST: Will return an array of array as np.array(np.array(), np.array(),...) each sub array contains [x,y,z] coord of a specific atom
	"""
	# let go of the atom names, the order is supposed to keep that information
	coords = np.array([np.array([float(y) for y in x.strip().split()[1:]]) for x in xyz_block[2:]], dtype=np.float64) # skip the two first lines that hold atom number and energy info
	return coords

def process_xyz_coordinates(xyzfile):
	"""
	PRE: Takes in an xyz file  that contains an xyz coordinate, each snapshot must have an equal number of atoms in the same order. They need to be pre aligned RMS wise.
	POST: Will process each snapshot and produce an average coordinate
	"""
	natoms = -1
	snapshot_list=[]
	with open(xyzfile, 'rb') as r:
		tempxyz=[]
		for i, line in enumerate(r):
			# if i==10000: break
			if i==0:
				natoms=line.strip()
			if line.strip()==natoms and tempxyz!=[]:
				if i<1000:
					atm_list=[x.strip().split()[0] for x in tempxyz[2:]]
				xyz_array=parse_xyz_block(tempxyz)
				snapshot_list.append(xyz_array)
				tempxyz=[]
				nsnaps=len(snapshot_list)
				if nsnaps%10==0: print "{}th snapshot processed".format(nsnaps)
			tempxyz.append(line)

	snapshot_array=np.array(snapshot_list[10000:])
	average_coord=average_xyz_coordinates(snapshot_array)
	generate_mol_from_xyz_and_pattern(average_coord, atm_list, fname=xyzfile+'_average')


def average_xyz_coordinates(xyz_array):
	"""
	PRE: takes in an array (xyz of a trajectory) of arrays (xyz of a single molecule) of arrays (xyz of single atom)
		THE SNAPSHOTS NEED TO BE ALIGNED, FIRST TEST USING OPENBABEL
	POST: returns the average array in xyz coords
	"""

	return np.mean(xyz_array, axis=0)

def generate_mol_from_xyz_and_pattern(xyz_array, atm_list, fname="average_mol.xyz"):
	"""
	PRE  : Takes in an array with xyz coordinates and an atom list
	POST : Will produce a molecule using both information
	"""
	with open(fname, "wb") as w:
		w.write(str(len(atm_list))+"\n")
		w.write("AVERAGE MOLECULE COORDINATES\n")
		for i, atm in enumerate(atm_list):
			w.write("{}      {}      {}      {}\n".format(atm, xyz_array[i][0], xyz_array[i][1], xyz_array[i][2]))












#################3
#
#
# This part is to average a trajectory of a cb complex as cylindrical distribution
# The trajectories need to be aligned according to CB then converted to cylindrical coordinates along the CB axis
# The trajectory is initially converted form xyz to sdf using openbabel
#
#################

def edit_xyz_file_to_add_format_charge(sdffile, chargeTAG="M  CHG  1   8   1"):
	"""
	PRE : XYZ trajs from cp2k don't contain the charge so converting them to openbabel don't actually add the charges which raises an error in rdkit
	POST : This method will modify the sdf file in question to add a charge tag at the bottom of each conformation to make it readable by the rdkit
	"""
	with open(sdffile, 'rb') as r:
		with open(sdffile+"-w-charge.sdf", 'wb') as w:
			for line in r:
				if "M  END" in line:
					w.write(chargeTAG +"\n"+line)
				else:
					w.write(line)

	return None

def align_trajectory_with_axes_according_to_cb(sdftraj):
	"""
	PRE:  Takes an xyz trajectory
	POST: Returns an rdkit mol with each conformer storing a trajectory snapshot, will align the conformers based on the
	"""

	supp = Chem.SDMolSupplier(sdftraj, removeHs=False)
	w = Chem.SDWriter(sdftraj+"-aligned.sdf")
	### GET THE BASE CONFORMER
	traj = supp[0]

	#### Get the alignment with the axes transform
	frags = Chem.GetMolFrags(traj, asMols=True)
	if frags[0].GetNumAtoms()<frags[1].GetNumAtoms():
		guest, cb = frags[0], frags[1]
	else:
		guest, cb = frags[1], frags[0]
	cantransform=Chem.rdMolTransforms.ComputeCanonicalTransform(cb.GetConformer())
	#### ALIGNS traj
	AllChem.TransformMol(traj, cantransform)
	w.write(traj)

	for mol in supp:
		AllChem.TransformMol(mol, cantransform)
		w.write(mol)
		# traj.AddConformer(mol.GetConformer())

def get_xyz_from_mol(mol):
		"""
		PRE: Takes in a mol rdkit that has one conformer (or at least the coordinats will be taken from the first one)
		POST: returns an ndarray with the xyz coordinates
		"""
		atom_list=[]
		conf=mol.GetConformer(0)
		for i in range(mol.GetNumAtoms()):
			atom_list.append(np.array(list(conf.GetAtomPosition(i))))
		return np.array(atom_list)

def xyz_to_cylindrical(xyzarray):
	"""
	PRE: Takes in a xyz array
	POST: returns an array that corresponds the cylindrical transform of those coordinates with the cylindrical axis being along the CB axis
	"""
	pass

def compute_cylindrical_transform_for_sdf_trajectory(sdffile):
	"""
	PRE:  Takes in sdf trajectory that has been pre aligned along with the axes according to the CB
	POST: Will produce an array of snapshots in cylindrical coordinates and probably plot it
	"""
	supp = Chem.SDMolSupplier(sdffile, removeHs=False)
	get_xyz_from_mol(supp[0])

def make_mass_matrix_from_atom_list(atom_list):
	"""
	PRE: Takes in an atom list
	POST: Will return the mass matrix associated with it
	the matrix is square and has 3n*3n elements and is diagonal for n atoms, the square root of the matrix is also returned
	"""
	dic_weights={'C':12,'H':1,'N':14, 'O':16}
	diag=[]
	for ats in atom_list:
		diag.extend([dic_weights[ats]*1.66e-27]*3) # those masses are atomic masses so the factor 1.66e-27

	mass_matrix = np.diag(diag)
	# sqrt_mass_matrix = scipy.linalg.sqrtm(mass_matrix) # in the case of diagonal matrices, it is just the square root of the elements element wise
	sqrt_mass_matrix=np.diag([x**.5 for x in diag])

	return mass_matrix, sqrt_mass_matrix


def make_covariance_matrix(xyzfile):
	"""
	PRE: Takes in an xyz file containing a trajectory
	POST: Will return the covariance of the atomic coordinates: For n atoms the matrix will be 3n*3n elements are each component of speed is an individual variable
	"""
	print xyzfile
	natoms = -1
	snapshot_list=[]
	with open(xyzfile, 'rb') as r:
		tempxyz=[]
		for i, line in enumerate(r):
			# if i==10000: break
			if i==0:
				natoms=line.strip()
			if line.strip()==natoms and tempxyz!=[]:
				if i<1000:
					atm_list=[x.strip().split()[0] for x in tempxyz[2:]]
				xyz_array=parse_xyz_block(tempxyz)
				xyz_vec=np.reshape(xyz_array, 3*int(natoms))
				snapshot_list.append(xyz_vec)
				tempxyz=[]
				nsnaps=len(snapshot_list)
				# if nsnaps%10==0: print "{}th snapshot processed".format(nsnaps)
			tempxyz.append(line)

	mass_matrix, sqrt_mass_matrix= make_mass_matrix_from_atom_list(atm_list)
	snapshot_array=np.array(snapshot_list).T
	covmat=np.cov(snapshot_array*1e-10) # Factor 1e-10 there because we fed angstrom in the procedure

	mass_weighted_covmat = np.matmul(np.matmul(sqrt_mass_matrix,covmat), sqrt_mass_matrix)

	eigenvalues, eigenvectors = scipy.linalg.eig(mass_weighted_covmat)
	k,T,hbar = 1.38e-23, 300, 1.054e-34 # those are the masses of
	frequencies = [np.sqrt(k*T/np.real(x)) for x in eigenvalues if x>0] # the frequencies are normally given in Hz but here we divide by the speed of light in cm-1/s to get cm-1 (2.9979e10)

	get_entropy_from_frequency_list(frequencies, k, T, hbar)


def get_entropy_from_frequency_list(frequency_list, k,T, hbar ):
	"""
	PRE  : Takes in a list in frequencies in hertz
	POST : returns a value of the entropy as given by the formulae of statistical thermodynamic
	"""
	alpha=hbar/k/T
	print 0.00198588*sum([alpha *x/(np.exp(alpha*x) - 1) -  np.log(1-np.exp(-alpha*x)) for x in frequency_list[6:]])  # The prefactor instead of k is the value of R in kcal/mol

if __name__ == "__main__":
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import statistics
	import numpy as np
	import numpy.core.multiarray
	from sklearn.decomposition import PCA
	from sklearn.cluster import KMeans
	import matplotlib.pyplot as plt
	import pandas as pd
	import scipy
	import glob

	# process_z_matrix_trajectory('cb6.inp-pos-1-aligned.gzmat')
	# for f in ['/home/macenrola/Documents/XYLENE/inputs/for-reaction-frozen-cb/MO-CB6.inp-pos-1-aligned-just-CB6.xyz',
	# 			'/home/macenrola/Documents/XYLENE/inputs/for-reaction-frozen-cb/MO-CB7.inp-pos-1-aligned-just-CB6.xyz',
	# 			'/home/macenrola/Documents/XYLENE/inputs/for-reaction-frozen-cb/MP-CB6.inp-pos-1-aligned-just-CB6.xyz',
	# 			'/home/macenrola/Documents/XYLENE/inputs/for-reaction-frozen-cb/MP-CB7.inp-pos-1-aligned-just-CB6.xyz']:
	# process_xyz_coordinates('/home/macenrola/Documents/XYLENE/base-systems-equilibrated/starting-trajectories-for-frozen/CB6-aligned-from-MP-centered.xyz')
	for f in glob.glob("/home/macenrola/Documents/XYLENE/Entropy-from-covariance/*xyz"):
		make_covariance_matrix(f)

	# align_trajectory_with_axes_according_to_cb("/home/macenrola/Documents/heptylamine/TRAJS/700-heptylamine.inp-pos-1.sdf-w-charge.sdf")
	# edit_xyz_file_to_add_format_charge("/home/macenrola/Documents/heptylamine/TRAJS/700-heptylamine.inp-pos-1.sdf")
	# compute_cylindrical_transform_for_sdf_trajectory("/home/macenrola/Documents/heptylamine/TRAJS/700-heptylamine.inp-pos-1.sdf-w-charge.sdf-aligned.sdf")
