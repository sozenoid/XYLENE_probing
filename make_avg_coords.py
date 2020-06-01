
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

def return_list_of_snapshots_and_atom_list(xyzfile):
	"""
	PRE: Takes in an xyz file
	POST: Returns a list of atoms and a list of snapshots of the trajectory as arrays
	"""
	natoms = -1
	snapshot_list=[]
	with open(xyzfile, 'rb') as r:
		tempxyz=[]
		for i, line in enumerate(r):
# =============================================================================
# 			if i==1000000: break
# =============================================================================
			if i==0:
				natoms=line.strip()
			if line.strip()==natoms and tempxyz!=[]:
				if i<1000:
					atm_list=[x.strip().split()[0] for x in tempxyz[2:]]
				xyz_array=parse_xyz_block(tempxyz)
				snapshot_list.append(xyz_array)
				tempxyz=[]
				nsnaps=len(snapshot_list)
				if nsnaps%10000==0: print "{}th snapshot processed".format(nsnaps)
			tempxyz.append(line)

	return atm_list, np.array(snapshot_list)[:]

def process_xyz_coordinates(xyzfile):
	"""
	PRE: Takes in an xyz file  that contains an xyz coordinate, each snapshot must have an equal number of atoms in the same order. They need to be pre aligned RMS wise.
	POST: Will process each snapshot and produce an average coordinate
	"""
	atm_list, snapshot_array= return_list_of_snapshots_and_atom_list(xyzfile)
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
	PRE: Takes in a xyz array where the cb has it symmetry axis along z or the last of the three xyz coordinates
	POST: returns an array that corresponds the cylindrical transform of those coordinates with the cylindrical axis being along the CB axis: it will return (r, theta, z)
	"""
	def convertxyzcyl(point):
		return np.array([np.sqrt(point[0]**2 + point[1]**2), np.arctan(point[1]/point[0]), point[2]])
	cylarray = np.array(map(convertxyzcyl, xyzarray))

	return cylarray

def get_probability_density_by_histogram(rlist, zlist):
	"""
	PRE : Takes in two lists that define a cloud of points, first list is x and the second y or (r, z) here
	POST: will plot the density
	"""
	from matplotlib.colors import LogNorm
	xedges=np.linspace(0,8, 100)
	yedges=np.linspace(-15,5, 100)

	H, xedges, yedges = np.histogram2d(rlist, zlist, bins=(xedges, yedges))
	H = H.T  # Let each row list bins with common y range.

	# fig = plt.figure(figsize=(7, 3))
	# ax = fig.add_subplot(131, title='imshow: square bins')
	plt.imshow(H, interpolation='nearest', origin='low'
	, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap="jet", norm=LogNorm())

	plt.colorbar()

	# plt.contour(H)
	#### RAW SCATTER IS TOO SHADY, TOO MANY POINTS
	# plt.scatter(rlist, zlist, s=1)
	plt.xlabel("RADIUS FROM CB SYMMETRY AXIS (ANGSTROM)")
	plt.ylabel("Z DISTANCE FROM ORIGIN (ANGSTROM)")

def compute_cylindrical_transform_for_sdf_trajectory(sdffile):
	"""
	PRE:  Takes in sdf trajectory that has been pre aligned along with the axes according to the CB
	POST: Will produce an array of snapshots in cylindrical coordinates and probably plot it
	"""

	########## TO PROCESS
	# supp = Chem.SDMolSupplier(sdffile, removeHs=False)
	# rlist = []
	# zlist = []
	# for i, mol in enumerate(supp):
	# 	print i
	# 	# if i==100: break
	#
	# 	temp_cylarray = xyz_to_cylindrical(get_xyz_from_mol(mol)).T
	# 	rlist.extend(temp_cylarray[0])
	# 	zlist.extend(temp_cylarray[2])

	###### FROM SAVED
	rlist, zlist = cPickle.load(open(sdffile, 'rb'))
	get_probability_density_by_histogram(rlist, zlist )

	plt.savefig(sdffile+'-rz.png')
	plt.show()
	
	############

	# cPickle.dump((rlist,zlist), open(sdffile+'-rzlists', 'wb'))
	# plt.show()

def make_mass_matrix_from_atom_list(atom_list):
	"""
	PRE: Takes in an atom list
	POST: Will return the mass matrix associated with it
	the matrix is square and has 3n*3n elements and is diagonal for n atoms, the square root of the matrix is also returned
	"""
	dic_weights={'C':12.,'H':1.,'N':14., 'O':16., "Cl":35., "F":19., "Br":79.}
	diag=[]
	for ats in atom_list:
		diag.extend([dic_weights[ats]*1.660540e-27]*3) # those masses are atomic masses so the factor 1.66e-27
# =============================================================================
# 		diag.extend([dic_weights[ats]]*3) # those masses are atomic masses so the factor 1.66e-27
# =============================================================================

	mass_matrix = np.diag(diag)
# =============================================================================
# 	sqrt_mass_matrix = scipy.linalg.sqrtm(mass_matrix) # in the case of diagonal matrices, it is just the square root of the elements element wise
# =============================================================================
	sqrt_mass_matrix=np.diag([x**.5 for x in diag])

	return mass_matrix, sqrt_mass_matrix


def make_covariance_matrix(xyzfile):
	"""
	PRE: Takes in an xyz file containing a trajectory, the trajectory needs to be centered and RMS aligned, otherwise spurious frequencies pop up
	POST: Will return the covariance of the atomic coordinates: For n atoms the matrix will be 3n*3n elements are each component of speed is an individual variable
	"""
	print xyzfile
	atm_list, snapshot_array= return_list_of_snapshots_and_atom_list(xyzfile)

	reshaped_snapshots = [np.reshape(x, len(atm_list)*3) for x in list(snapshot_array)]
	mass_matrix, sqrt_mass_matrix= make_mass_matrix_from_atom_list(atm_list)
	snapshot_array=np.array(reshaped_snapshots).T
	covmat=np.cov(snapshot_array*1e-10) # Factor 1e-10 there because we fed angstrom in the procedure
# =============================================================================
# 	covmat=np.cov(snapshot_array) # Factor 1e-10 there because we fed angstrom in the procedure
# =============================================================================

	mass_weighted_covmat = np.dot(np.dot(sqrt_mass_matrix,covmat), sqrt_mass_matrix)
	print mass_weighted_covmat
	eigenvalues, eigenvectors = scipy.linalg.eig(mass_weighted_covmat)
	k,T,hbar = 1.3807e-23, 300, 1.054e-34 # those are the masses of kB T and hbar
	negeigcount = sum(n < 0 for n in eigenvalues) # Count the negative frequencies to remove them
	print sorted(eigenvalues)
	frequencies = [(k*T/np.real(x))**.5 for x in eigenvalues if np.real(x)>0] # the frequencies are normally given in Hz but here we can divide by the speed of light in cm-1/s to get cm-1 (2.9979e10)
	print "Frequences are in cm-1", sorted([x/2.9979e10 for x in frequencies[6-negeigcount:]])
	get_entropy_from_frequency_list(frequencies[6-negeigcount:], k, T, hbar)
	get_energy_from_frequency_list(frequencies[6-negeigcount:], k, T, hbar)


def get_frequencies_quasiharmonic(xyzfile):
	"""
	PRE: Takes an xyz file with units in Angstrom, spaced by 1 fs 
	POST: Will compute the frequencies using the approach described by 
	Strachan, Alejandro. "Normal modes and frequencies from covariances in molecular dynamics or Monte Carlo simulations." The Journal of chemical physics 120.1 (2004): 1-4.
	"""
	c = 2.99e10 #cm/s
	print xyzfile
	atm_list, snapshot_array= return_list_of_snapshots_and_atom_list(xyzfile)
	reshaped_snapshots = [np.reshape(x, len(atm_list)*3) for x in list(snapshot_array)][-30000:]
	mass_matrix, sqrt_mass_matrix= make_mass_matrix_from_atom_list(atm_list)
	snapshot_array = np.array(reshaped_snapshots).T*1e-10
	snapshot_speed = np.gradient(snapshot_array, 1e-15, edge_order=1, axis=1)
	snapshot_acc   = np.gradient(snapshot_speed, 1e-15, edge_order=1, axis=1)

	mass_rescaled_pos = [(x-np.mean(x))*y/2**.5 for x, y in zip(snapshot_array, np.diag(sqrt_mass_matrix))]
	#mass_rescaled_pos = [(x)*y/2**.5 for x, y in zip(snapshot_array, np.diag(sqrt_mass_matrix))]
	mass_rescaled_speed = [(x)*y/2**.5 for x, y in zip(snapshot_speed, np.diag(sqrt_mass_matrix))]
	mass_rescaled_acceleration = [(x)*y/2**.5 for x, y in zip(snapshot_acc, np.diag(sqrt_mass_matrix))]
	
# =============================================================================
# 	covmat_pos   = np.cov(np.matmul(mass_matrix, mass_rescaled_pos))
# 	covmat_speed = np.cov(np.matmul(mass_matrix,mass_rescaled_speed))
# 	covmat_acc   = np.cov(np.matmul(mass_matrix,mass_rescaled_acceleration))
# =============================================================================
	covmat_pos   = np.cov(mass_rescaled_pos)
	covmat_speed = np.cov(mass_rescaled_speed)
	covmat_acc   = np.cov(mass_rescaled_acceleration)
	
	eigenvalues_pos, eigenvectors_pos     = scipy.linalg.eig(covmat_pos)
	eigenvalues_speed, eigenvectors_speed = scipy.linalg.eig(covmat_speed)
	eigenvalues_acc, eigenvectors_acc     = scipy.linalg.eig(covmat_acc)
	

# =============================================================================
# 	for e in eigenvalues_pos:
# 		print (1.38e-23*300/e)**.5/c/2/math.pi
# =============================================================================
		
# =============================================================================
# 	mass_rescaled_speed_snap = np.array(mass_rescaled_speed).T
# 	mass_rescaled_pos_snap = np.array(mass_rescaled_pos).T
# 	for e in eigenvectors_speed[0:]:
# 		e_projected_speed_snap = []
# 		e_projected_pos_snap = []
# 		for snap in mass_rescaled_speed_snap:
# 			e_projected_speed_snap.append([x*y for x, y in zip(e, snap)])
# 		for snap in mass_rescaled_pos_snap:
# 			e_projected_pos_snap.append([x*y for x, y in zip(e, snap)])
# 		
# 		print 1 / (2*math.pi)/c * (np.mean(e_projected_speed_snap)**2 / np.var(e_projected_pos_snap))**.5 
# =============================================================================
		
	
# =============================================================================
# 	print len(eigenvectors_pos)
# 	print eigenvectors_speed
# 	print eigenvectors_acc
# =============================================================================
	
	#frequencies = sorted([(x/y)**.5 / (2*math.pi) for x,y in zip(eigenvalues_speed, eigenvalues_pos)])[:]
	frequencies = sorted([(x/y)**.25 / (2*math.pi) for x,y in zip(eigenvalues_acc, eigenvalues_pos)])[:]
	#frequencies = sorted([(x/y)**.5 / (2*math.pi) for x,y in zip(eigenvalues_acc, eigenvalues_speed)])[:]
	
	#print [x/c for x in frequencies]
	
	k,T,hbar = 1.3807e-23, 300, 1.054e-34 # those are the masses of kB T and hbar

	get_entropy_from_frequency_list(frequencies[:], k, T, hbar*(2*math.pi))
	
def get_entropy_from_frequency_list(frequency_list, k,T, hbar ):
	"""
	PRE  : Takes in a list in frequencies in hertz
	POST : returns a value of the entropy as given by the formulae of statistical thermodynamic
	"""
	alpha=hbar/k/T
# =============================================================================
# 	print k*sum([alpha *x/(np.exp(alpha*x) - 1) -  np.log(1-np.exp(-alpha*x)) for x in sorted(frequency_list)])  # The prefactor instead of k is the value of R in kcal/mol, let go of rotational and translational frequencies
# =============================================================================
 	print 0.00198588*sum([alpha *x/(np.exp(alpha*x) - 1) -  np.log(1-np.exp(-alpha*x)) for x in sorted(frequency_list)])  # The prefactor instead of k is the value of R in kcal/mol, let go of rotational and translational frequencies

def get_energy_from_frequency_list(frequency_list, k,T, hbar ):
	"""
	PRE  : Takes in a list in frequencies in hertz
	POST : returns a value of the entropy as given by the formulae of statistical thermodynamic
	"""
	alpha=hbar/k/T
# =============================================================================
# 	print sum([(hbar*x/2 + hbar*x/(np.exp(alpha*x) - 1)) for x in sorted(frequency_list)])  # The prefactor instead of k is the value of R in kcal/mol, let go of rotational and translational frequencies
# =============================================================================
	print sum([(hbar*x/(np.exp(alpha*x) - 1)) for x in sorted(frequency_list)])  # The prefactor instead of k is the value of R in kcal/mol, let go of rotational and translational frequencies


def convert_position_snapshots_to_velocity_snapshots(xyztraj, dt=1e-15):
	"""
	PRE  : Takes in an xyz trajectory where the snaps are separated by a small time step (for the velocity to be "instant")
	POST : will compute the velocity based on the time step size and also provide an estimate of the kinetic energy at each step
	"""
	print xyztraj
	atm_list, snapshot_pos= return_list_of_snapshots_and_atom_list(xyztraj) # in Angstroms
	print snapshot_pos.shape
	snapshots_vels = np.diff(snapshot_pos, axis=0) # in Angstrom/fs
	print snapshots_vels.shape
	dic_weights={'C':12,'H':1,'N':14, 'O':16, 'F':19.}
	conv=(1e-10/1e-15)**2*1.66e-27*6.022e23/4.184/1e3 # (Ang/fs)**2*amu*mol/kcal
	ekin = [0.5*conv*sum([y[0]**2*y[1] for y in zip(np.reshape(x, x.shape[0]*x.shape[1]), [dic_weights[i] for k in atm_list for i in [k]*3])]) for x in snapshots_vels]
	plt.plot(ekin)
	plt.show()
	
def power_spectrum_from_velocityxyz(xyzfile,l):
	"""
	PRE : Takes in an xyz file containing velocity (most likely in bohr/au_time)
	POST: Will produce the power spectrum
	"""
	print "Getting the snaps"
	atm_list, snapshot_array= return_list_of_snapshots_and_atom_list(xyzfile)
	#
	
	print "Got the snaps, now reshaping"
	snapshot_array = [np.reshape(x, len(atm_list)*3) for x in list(snapshot_array)][-25000:]
	snapshot_array = [x*1e-10 for x in snapshot_array]
	
	#### CONVERT SNAPS OF POSITION TO SPEED
	snapshot_speed = np.gradient(snapshot_array, 1e-15, edge_order=2, axis=0)
	snapshot_array = snapshot_speed
	mass_matrix, sqrt_mass_matrix= make_mass_matrix_from_atom_list(atm_list)
	#### END CONVERSION
	
	print "reshaping level 1"
	snapshot_array=np.array(snapshot_array).T
	print "Reshaping level 2"
	total_correlation = None
	print "Computing the correlation"
	for i, component in enumerate(snapshot_array):
# =============================================================================
# 		print "Step {} in the correlation computation".format(i)
# =============================================================================
# =============================================================================
# 		if i==0: total_correlation=np.correlate(component, component, mode='full')*np.diag(mass_matrix)[i]
# 		else: total_correlation = total_correlation + np.correlate(component, component, mode='full')*np.diag(mass_matrix)[i]
# =============================================================================
		
# =============================================================================
# 		temp = np.array([np.abs(A)**2 for A in numpy.fft.fftshift(numpy.fft.fft(component))])*np.diag(mass_matrix)[i]
# 		if i==0: total_correlation = temp
# 		else: total_correlation = total_correlation + temp
# 	print "Summing the correlation"
# =============================================================================
		n=10000
		temp = scipy.signal.welch(component, fs=1e15, nperseg=n)[1]*np.diag(mass_matrix)[i]
		if i==0: total_correlation = temp
		else: total_correlation = total_correlation + temp
	print "Summing the correlation"	
	
	
	
# =============================================================================
# 	total_correlation = np.sum(np.array(total_correlation), axis=0)
# =============================================================================
	#
# =============================================================================
# 	cPickle.dump(total_correlation, open('/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/totcor', 'wb'))
# =============================================================================
	# total_correlation=cPickle.load(open('/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/totcor', 'rb'))
# =============================================================================
# 	print total_correlation
# =============================================================================
# =============================================================================
# 	plt.plot(total_correlation)
# =============================================================================
	magnitude  =[np.abs(A)**2 for A in numpy.fft.fftshift(numpy.fft.fft(total_correlation))]
	angle      =[np.angle(A) for A in numpy.fft.fftshift(numpy.fft.fft(total_correlation))]
	
	#### WELCH METHOD
	freqs = np.linspace(0, 1e15/2.99e10/2, 1+(n)/2)
# =============================================================================
# 	plt.plot(freqs, total_correlation, label=l) # gives the spectrum with x axis in cm-1
# 	plt.legend()
# 	plt.show()
# =============================================================================
	#### END WELCH METHOD
	#### DIVIDE BY TWO FOR 1 fs, divide by ONE for .5 fs, return only magnitude by default 
# =============================================================================
# 	freqs = np.linspace(-(1e-15)**(-1) /2.9979e10/2, (1e-15)**(-1)/2.9979e10/2, len(total_correlation))
# =============================================================================
	
# =============================================================================
# 	plt.plot(freqs, total_correlation, label=l) # gives the spectrum with x axis in cm-1
# 	plt.legend()
# 	plt.show()
# =============================================================================

	return  freqs, magnitude, total_correlation 
# =============================================================================
# 	significant = [x for x in zip(freqs, magnitude, angle) if x[1]>0.001]
# 	print significant
# 	plt.plot(freqs, magnitude, color='red')
# 	plt.stem([x[0] for x in significant], [x[2] for x in significant]) # gives the spectrum with x axis in cm-1
# =============================================================================
# =============================================================================
# 	plt.savefig('/home/macenrola/Desktop/last_20000.pdf')
# 	plt.show()
# =============================================================================
# =============================================================================
# 	
# def split_velocity_file(xyz_vel, natoms=(127,108,19)):
# 	"""
# 	=========
# 	GO FOR IT; there was a bug in the file splitting part i'd say, probably better to keep trajectories short 
# 	=========
# 	PRE  : Takes in a velocity file formatted as an open babel xyz file, the atom breakdown is given in natoms, the smaller fragment is on top of the stack
# 	POST : returns two xyz_velocity files and adapts the headers for the two molecule components
# 	"""
# 	atm_list, snapshot_pos= return_list_of_snapshots_and_atom_list(xyz_vel) # in Angstroms
# 	with open(xyz_vel+"-small-frag.xyz", "wb") as w:
# 		for i, els in enumerate(snapshot_pos):
# 			w.write(str(natoms[-1])+"\n")
# 			w.write("i =   {}; and this is the small fragment\n".format(i))
# 			for coord in zip(atm_list[:natoms[-1]], els):
# 				w.write("{}    {}    {}    {}\n".format(coord[0], *coord[1]))
# 				
# 
# 	with open(xyz_vel+"-large-frag.xyz", "wb") as w:
# 		for i, els in enumerate(snapshot_pos):
# 			w.write(str(natoms[1])+"\n")
# 			w.write("i =   {}; and this is the small fragment\n".format(i))
# 			for coord in zip(atm_list[natoms[-1]:], els):
# 				w.write("{}    {}    {}    {}\n".format(coord[0], *coord[1]))
# =============================================================================
				

def split_velocity_file(xyz_vel, natoms_small_frag=19):
	"""
	=========
	GO FOR IT; there was a bug in the file splitting part i'd say, probably better to keep trajectories short 
	=========
	PRE  : Takes in a velocity file formatted as an open babel xyz file, the atom breakdown is given in natoms, the smaller fragment is on top of the stack
	POST : returns two xyz_velocity files and adapts the headers for the two molecule components
	"""
	small_frag = xyz_vel+"-small_frag.xyz"
	large_frag = xyz_vel+"-large_frag.xyz"
	marker=None
	line_list=[]
	with open(small_frag, "wb") as ws:
		with open(large_frag, "wb") as wb:
			with open(xyz_vel, "rb") as r:
				for i, line in enumerate(r):
					if i%10000==0: print i
					line_list.append(line)
					if i==0:
						marker=line
					if marker == line and i>0:
						ws.writelines("     {}\n".format(natoms_small_frag))
						ws.writelines(line_list[1])
						ws.writelines(line_list[2:natoms_small_frag+2])
						wb.writelines("     {}\n".format(int(line_list[0]) - natoms_small_frag))
						wb.writelines(line_list[1])
						wb.writelines(line_list[natoms_small_frag+2:-1])
						line_list=[line]

def split_xyz_file_in_chuncks(xyz_file, nsnaps=10000):
	"""
	PRE  : Will take in an xyz trajectory formatted as in openbabel like
	=====================
	     127
 i =   954829, time =   954829.000, E =      -508.3449994676
    =====================
	It is assumed that the very first line of the file serves as a marker to separate the various snaps
	POST : Will produce a series of files named consecutively with each having nsnaps snapshots, the total number is the number of snaps in the 
	original file divided by nsnaps
	"""
	parts_name = xyz_file+"-part-{}-{}.xyz"
	csnap=-1
	marker=None
	prev_line=""
	with open(xyz_file, "rb") as r:
		for i, line in enumerate(r):
# =============================================================================
# 			if i==10000: break
# =============================================================================
			if i==0:
				marker=line
			if marker == line:
				csnap+=1
				if csnap%nsnaps==0:
					if i>0:
						w.close()
					w = open(parts_name.format(csnap, csnap+nsnaps), "wb")
			w.write(line)
				

def plot_sequence_of_ffts(MAG_LIST):
	"""
	PRE  : Takes in a list of correlation as a list of lists, assuming they correspond to the same frequency x axis 
	POST : Will plot a log scale evolution of the major frequency present at the last snapshot
	"""
	MAG_LIST = cPickle.load(open(MAG_LIST, "rb"))
	target_frequencies = [x>0.1 for x in MAG_LIST[-2]]
	MAG_ARRAY = np.asarray(MAG_LIST[10:])
	print MAG_ARRAY.shape
	mag_target = []
	for i, (target, mag) in enumerate(zip(target_frequencies, list(MAG_ARRAY.T))):
		if target:
			mag_target.append(mag)
	mag_target = list(np.array(mag_target))
	
	print mag_target
	
	for els in mag_target:
		plt.semilogy(els)
	plt.show()
	

def get_kinetic_energy_from_velocity_file(xyz_vel):
	"""
	PRE : Takes in a velocity file formatted as per open babel
	POST:
	"""
	atm_list, snapshot_vel= return_list_of_snapshots_and_atom_list(xyz_vel) # in Angstroms
	dic_weights={'C':12,'H':1,'N':14, 'O':16}
	conv=911.447*627.509 # (Ang/fs)**2*amu*mol/kcal the 911 factor is INCLUDING THE FACTOR 2, the proper conversion is amu/m_e
	ekin = [conv*sum([y[0]**2*y[1] for y in zip(np.reshape(x, x.shape[0]*x.shape[1]), [dic_weights[i] for k in atm_list for i in [k]*3])]) for x in snapshot_vel]
	cPickle.dump(ekin, open(xyz_vel+"-EKIN","wb"))
	plt.plot(ekin, linewidth=0.2)
	plt.show()
	
def get_kinetic_energy_from_position_file(xyz_pos):
	"""
	PRE : Takes in a velocity file formatted as per open babel
	POST:
	"""
	atm_list, snapshot_array= return_list_of_snapshots_and_atom_list(xyz_pos) # in Angstroms

	snapshot_array = [np.reshape(x, len(atm_list)*3) for x in list(snapshot_array)][:]
	snapshot_array = [x*1e-10 for x in snapshot_array]
	
	#### CONVERT SNAPS OF POSITION TO SPEED
# =============================================================================
# 	snapshot_speed = np.gradient(snapshot_array, 1e-15, edge_order=1, axis=1)
# =============================================================================
	snapshot_speed = np.diff(snapshot_array, axis=0)/1e-15
	mass_matrix, sqrt_mass_matrix= make_mass_matrix_from_atom_list(atm_list)
	#### END CONVERSION
	conv = 6.02e23/4.18/1e3# J/molecule to kcal/mol
	ekin = [conv*sum([x**2*y/2. for x,y in zip(k, np.diag(mass_matrix))]) for k in snapshot_speed]
	ekin = [x for x in ekin if x<1000]
	print xyz_pos, "has variance of ekin of {}".format(np.var(ekin[10000:]))
# =============================================================================
# 	ekin = [conv*sum([y[0]**2*y[1] for y in zip(np.reshape(x, x.shape[0]*x.shape[1]), [dic_weights[i] for k in atm_list for i in [k]*3])]) for x in snapshot_vel]
# =============================================================================
	cPickle.dump(ekin, open(xyz_pos+"-EKIN","wb"))
# =============================================================================
# 	plt.plot(ekin, linewidth=0.2)
# =============================================================================
# =============================================================================
# 	plt.show()
# =============================================================================

def plot_E_kin(fekin, legend, symb):
	"""
	PRE : Takes in an Ekin file
	POST: Plots it
	"""
	
	ekin_all=cPickle.load(open(fekin, "rb"))
	speedup=10
	ekin_sampled = ekin_all[0::speedup]
	x = np.linspace(0, len(ekin_all)/1000.0, len(ekin_all))[0::speedup]
	plt.plot(x,ekin_sampled, linewidth=0.2, alpha=0.2)
	width=1000
	moving_av = np.convolve(ekin_sampled, np.ones(width)/width)
	plt.plot(x[width:], moving_av[width-1:-width], linewidth=1, linestyle=symb, label=legend)
	plt.axhline(1,linestyle='--', linewidth=0.4, color='k', xmin=0.5)
	plt.axhline(0,linestyle='--', linewidth=0.4, color='k', xmin=0.5)
	plt.xlabel(r"Time [ps]")
	plt.ylabel(r"$E_{kin}$ [kcal/mol]")
	plt.legend()
	plt.savefig(fekin+'.pdf')
# =============================================================================
# 	plt.close()
# =============================================================================


def plot_single_spectrum(fname, i=0):
	"""
	PRE: Takes a list of list with the magnitude of the fft of the velocity autocorrelation function computed over different time frames 
	POST : Will print the spectrum corresponding to the ith element in the list 
	"""
	maglist=cPickle.load(open(fname, "rb"))
	freqs = np.linspace(-1e15/2.9979e10/2, 1e15/2.9979e10/2, len(maglist[i]))
	ekin = np.trapz(maglist[i])
	plt.close()
	plt.plot(freqs, [x/ekin for x in maglist[i]])
	plt.show()
	

def plot_three_stacked_spectra(fname, indices, fname_10):
	"""
	PRE   : Takes in a MAG file with three indices
	POST  : Will plot the three indicies stacked 
	"""
	maglist=cPickle.load(open(fname, "rb"))
	yticks = [10**x for x in range(-16,1)]
	xticks = range(0,6000, 500)
	nplots=len(indices)
	freqs = np.linspace(0, 1e15/2.9979e10/2, len(maglist[-1])/2)
	fig, ax = plt.subplots(nplots, 1, sharex='col', sharey='row')
	for i in range(nplots):
# =============================================================================
# 		if i==0: 
# 			maglistsmall=cPickle.load(open(fname_10, "rb"))
# 			freqsmall = np.linspace(0, 1e15/2.9979e10/2, len(maglistsmall[0])/2)
# 			ekin = np.max(maglistsmall[0])	
# 			ax[0].semilogy(freqsmall, [x/ekin for x in maglistsmall[0][len(maglistsmall[0])/2+1:]], label="{}-{} ps".format(0, 10), linewidth=0.4)
# 		else:
# =============================================================================
		ekin = np.max(maglist[indices[i]])	
		ax[i].semilogy(freqs, [x/ekin for x in maglist[indices[i]][len(maglist[indices[i]])/2+1:]], label="{}-{} ps".format(indices[i]*100, (indices[i]+1)*100), linewidth=0.4)

		ax[i].get_yaxis().set_minor_formatter(plt.NullFormatter())
		ax[i].set_ylim((10**-16,1.2e1))
		ax[i].set_xlim((0,6e3))
		ax[i].set_yticks(yticks[::6])
		ax[i].set_xticks(xticks[::2])
		ax[i].set_yticks(yticks[::4], minor=True)
		ax[i].grid(which='minor')
		ax[i].legend(loc='east', fontsize=18)
	
	fig.text(0.03, 0.5, 'Intensity [a.u.]', ha='center', va='center', rotation='vertical')
# =============================================================================
# 	plt.ylabel("Intensity [a.u.]")
# =============================================================================
	plt.xlabel(r"Frequency [cm$^{-1}$]")
	plt.autoscale(enable=False, tight=True)

	font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

	plt.rc('font', **font)
	plt.show()
	
def make_mtd_time_plot(time_plot_file):
	"""
	PRE  : Takes in a summary file for the isomerisation 
	POST : Gives out a double plot for vacuum and CBs
	"""
# =============================================================================
# 	GETS A USEABLE DIC FOR PLOTTING
# =============================================================================
	R=1.99E-03
	lnkkb_h=23.75951822
	with open(time_plot_file , "rb") as r:
		lines = r.readlines()
	split_res = [x.strip().split("KS") for x in lines]
	res_dic = {}
	for el in split_res:
		if el==['']: continue
		id_system = el[0][-13:-1]
		key = id_system[:6] 
		if key in res_dic:
			res_dic[key].append((id_system[-5:], el[1][1:]))
		else:
			res_dic[key] = [(id_system[-5:], el[1][1:])]
		
	for k in res_dic:
		res_dic[k] = list(zip(*res_dic[k]))
		
# =============================================================================
# 	PLOTS IT
# =============================================================================
	fig, ax = plt.subplots(2, 1, sharex='col', sharey='row')
	T=[300, 400, 500, 600, 700]
	for (i,k),s in zip(enumerate(res_dic), line_styles()):
		print i,k,s
		xtrend = np.linspace(10e3*0.001,10e3*0.004,100)
		line = res_dic[k]
		plottable = [10e3*1/float(x) for x in line[0]], [np.log(1/float(x[0])/float(x[1])) for x in zip(line[0], line[1])]
		if 'CB' in k:
			ax[1].scatter(plottable[0], plottable[1], marker=i, label=k, s=100)
			slope, intercept, r_value, p_value, std_err = scipy.stats.linregress([x[:] for x in plottable])
			ax[1].plot(xtrend, [x*slope+intercept for x in xtrend], '--', linewidth=1, linestyle=s)
			ax[1].set_xlim((10e3*0.001,10e3*0.004))
			ax[1].set_ylabel(r"ln($\frac{k}{T}$)")
# =============================================================================
# 			ax[1].legend(loc='lower left')
# =============================================================================
			ax[1].legend(loc='lower left', bbox_to_anchor=(-0.1, 1.05, 1.2, 0.1),
			          fancybox=True, shadow=False, ncol=4, mode='expand')
			ax[1].grid(True, alpha=0.2)
			print k, slope, intercept, -slope*R, (intercept-lnkkb_h)*R
			print '&'.join([str((x*np.exp(slope/x + intercept))**-1) for x in T])
		elif "vac" in k:
			ax[0].scatter(plottable[0], plottable[1], marker=i, label=k, s=100)
			slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(plottable)
			ax[0].plot(xtrend, [x*slope+intercept for x in xtrend], '--', linewidth=1, linestyle=s)
			ax[0].set_ylabel(r"ln($\frac{k}{T}$)")
# =============================================================================
# 			ax[0].legend(loc='lower left')
# =============================================================================
			ax[0].legend(loc='lower left', bbox_to_anchor=(-0.1, 1.05, 1.2, 0.1),
			          fancybox=True, shadow=False, ncol=4, mode='expand')
			ax[0].grid(True, alpha=0.2)
			print k, slope, intercept, -slope*R, (intercept-lnkkb_h)*R
			print '&'.join([str((x*np.exp(slope/x + intercept))**-1) for x in T])
	plt.xlabel(r"1/T 10$^{-3}$ [K$^{-1}$]")
	
	font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 17}

	plt.rc('font', **font)
	plt.tight_layout()
	plt.show()

def make_mtd_popping_rate_plot(time_plot_file):
	"""
	PRE  : Takes in a summary file for the popping 
	POST : Gives out a triple plot for o m p xylenes
	"""
# =============================================================================
# 	GETS A USEABLE DIC FOR PLOTTING
# =============================================================================
	R=1.99E-03
	lnkkb_h=23.75951822
	with open(time_plot_file , "rb") as r:
		lines = r.readlines()
	split_res = [x.strip().split("KS") for x in lines]
	res_dic = {}
	for el in split_res:
		if el==['']: continue
# =============================================================================
# 		id_system = el[0][-19:-1]
# =============================================================================
		id_system = el[0][-23:-1]
		key = '-'.join(id_system.split("-")[:-1])
		if key[0]=='-': key=key[1:]
		print el, key, id_system
		if key in res_dic:
			res_dic[key].append((id_system.split("-")[-1], el[1][1:]))
		else:
			res_dic[key] = [(id_system.split("-")[-1], el[1][1:])]
		
	for k in res_dic:
		res_dic[k] = list(zip(*res_dic[k]))
		
# =============================================================================
# 	PLOTS IT
# =============================================================================
	fig, ax = plt.subplots(3, 1, sharex='col', sharey='row')
	ref_dic={'o':(0, r'$o$-xylene'),'m':(1,r'$m$-xylene'),'p':(2,r'$p$-xylene')}
	lstyle=line_styles()
	for (i,k),s in zip(enumerate(sorted(res_dic)), lstyle):
		print i, k,s 
		j = ref_dic[k[0]][0]
		xtrend = np.linspace(0.001,0.004,100)
		line = res_dic[k]
		plottable = [1/float(x) for x in line[0]], [np.log(1/float(x[0])/float(x[1])) for x in zip(line[0], line[1])]
		ax[j].scatter(plottable[0], plottable[1], marker=6+i%2, label="{} in {}".format(ref_dic[k[0]][1], k[-3:].upper()))
		slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(plottable)
		ax[j].plot(xtrend, [x*slope+intercept for x in xtrend], '--', linewidth=1, linestyle=lstyle[i%2])
		ax[j].set_ylabel(r"ln($\frac{k}{T}$)")
		ax[j].legend(loc='lower left')
		ax[j].grid(True, alpha=0.2)
		ax[1].set_xlim((0.001,0.004))
		print k, slope, intercept, -slope*R, (intercept-lnkkb_h)*R
		T=700
		print (T*np.exp(slope/T + intercept))**-1

	plt.xlabel(r"1/T [K$^{-1}$]")
	plt.tight_layout()
	plt.show()

def make_KS_plot(KS_PLOTTING_LISTS, ax):
	"""
	PRE: Takes in a list of escape times and an axis, it will fill the axis with the data
	POST:Plots the KS test on the same pic, by fitting the cumulative distribution of the exponential to the curve then computing the KS test 
	"""
	def f_to_fit(x, scale):
		return 1-np.exp(-x/scale)
	
	ref_dic={'o':(0, r'$o$-xylene'),'m':(1,r'$m$-xylene'),'p':(2,r'$p$-xylene')}
	for f,ls in zip(sorted(KS_PLOTTING_LISTS)[:], line_styles()):
		name = f.split("/")[-1]
		print name
		with open(f, "rb") as r:
			time_dist=[float(x.strip()) for x in r.readlines()]
		sortedist = sorted(time_dist)
		p = 1. * np.arange(len(sortedist)) / (len(sortedist) - 1)
		#
		#
		xn = np.geomspace(sortedist[0], sortedist[-1], 200)
		p0 = np.median(sortedist)
		popt, pcov = scipy.optimize.curve_fit(f_to_fit, sortedist, p, p0=p0)
		
		ks_test_res=scipy.stats.kstest(rvs=time_dist, cdf='expon', args=(0,popt), N=len(time_dist))
		print ks_test_res[1]
	# FOR RATES OF ISOMERIZATION
		ax.semilogx(sortedist, p)
		ax.semilogx(xn, f_to_fit(xn, *popt), linestyle=ls, label=name.split('-')[-2][:-2]+" K")
	# sort both labels and handles by labels
	ax.set_xlim((1e-12,1e8))
	ax.text(0.02, 0.5, name.split("dump")[-1][1:7],
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes)
	ax.grid(True, alpha=0.2)
	# END RATES OF ISOMERIZATION

# =============================================================================
# 	# FOR RATES OF POPPING
# 		ax.semilogx(sortedist, p, label=name.split('-')[-2])
# 		#ax.semilogx(xn, f_to_fit(xn, *popt), linestyle=ls, label=name.split('-')[-2])
# 	ax.set_xlim((1e-12,1e8))
# 	ax.text(1, 0.5, "{0} in {1}".format(ref_dic[name.split('-')[-4][0]][1], name.split('-')[-3].upper()), 
# 							horizontalalignment='right',
# 							verticalalignment='top',
# 							transform=ax.transAxes)
# 	ax.grid(True, alpha=0.2)	
# 	# END RATES OF POPPING
# =============================================================================

# =============================================================================
# 	ax.set_xlabel('Time to reaction [s]')
# =============================================================================
	ax.set_ylabel('CDF')
	
def stack_ks_plots(list_of_list_of_escape_times):
	"""
	PRE: Takes in a list of lists of time distributions
	POST: Will read all the lists and create a stacked plot with it
	"""
	fig, ax = plt.subplots(6, 1, sharex="col")
	for i, l in enumerate(list_of_list_of_escape_times):
		make_KS_plot(l, ax[i])
	handles, labels = ax[i].get_legend_handles_labels()
	labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=True)) 
	plt.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5,8.9), ncol=5)
# =============================================================================
# 	labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: float(t[0]), reverse=False)) 
# 	plt.legend(handles, labels, loc="upper right", bbox_to_anchor=(1.3,7.5), ncol=1)
# =============================================================================
	plt.tight_layout()
	plt.xlabel('Time to escape [s]')
	
	
def line_styles():
	linestyles = dict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
	
	return sorted(linestyles.values())


def plot_instantaneous_binding_energy(complexfile, cbfile, guestfile):
	"""
	PRE  : Takes in three files containing the complex, CB and the guest's energy  
	POST : Will plot the instantaneous binding energy
	"""
	with open(complexfile, "rb") as r:
		complexlist = r.readlines()
		
	with open(cbfile, "rb") as r:
		cblist = r.readlines()
		
	with open(guestfile, 'rb') as r:
		guestlist = r.readlines()
	
	complexdata = list(zip(*[(x.strip().split(',')) for x in complexlist]))
	cbdata = list(zip(*[(x.strip().split(',')) for x in cblist]))
	guestdata = list(zip(*[x.strip().split(',') for x in guestlist]))
	
	complexdata = [float(x)/1000 for x in complexdata[0]],[float(x) for x in complexdata[1]]
	cbdata = [float(x.split('down')[1].split('-')[0])/1000+20 for x in cbdata[0]], [float(x) for x in cbdata[1]]
	guestdata = [float(x.split('down')[1].split('-')[0])/1000+20 for x in guestdata[0]], [float(x) for x in guestdata[1]]
	
	cbdatasorted = sorted(cbdata[0]), [k for _,k in sorted(zip(cbdata[0], cbdata[1]))]
	guestdatasorted = sorted(guestdata[0]), [k for _,k in sorted(zip(guestdata[0], guestdata[1]))]
	lw=0.2
	fig, axes = plt.subplots(2,2, sharex="col")
	for dtset, zeropoint, ax in zip([complexdata, cbdatasorted, guestdatasorted], [508.3467427, 467.013497, 41.209], axes.flat):
		print zeropoint
		ax.plot(dtset[0], [(x+zeropoint)*627.5 for x in dtset[1]], linewidth=lw)
		ax.set_ylabel(r"$E_{pot}$ - $E_0$ [kcal/mol] ")
		ax.set_xlabel(r"Time [ps]")
		ax.grid(alpha=0.2)
	
	ebindax=axes.flat[-1]
	ebindax.plot(guestdatasorted[0],[-(x[0]-x[1]-x[2])*627.5 for x in zip(complexdata[1], cbdatasorted[1], guestdatasorted[1])],linewidth=lw)
	ebindax.set_ylabel(r"$E_{bind}$ [kcal/mol]")
	ebindax.set_xlabel(r"Time [ps]")
	plt.tight_layout()
	ebindax.grid(alpha=0.2)
	
def plot_ener(num, enerfile):
	"""
	PRE : Takes in an ener file
	POST: Plots it
	"""
	with open(enerfile, "rb") as r:
		energ = r.readlines()
	ids = energ[0].strip().split("       ")
	energ = list(zip(*[(x.strip().split('      ')) for x in energ[1:]]))
	
	print len(energ)
	float_energ=[]
	for ls in energ:
		float_energ.append([float(x) for x in ls])
	
	print zip(range(len(energ)), ids)
	avgnum = sum(float_energ[num])/len(float_energ[num])
	plt.plot(float_energ[1], [(x-avgnum)*627.5 for x in float_energ[num]])
	

def plot_sp_not_stable_output(resfile="/home/macenrola/Documents/XYLENE/inputs/SP-DONT-WORK/resfile_with_convergence/resfile"):
	"""
	PRE: Takes in a resfile containing the name of the files of a MO TS in CB7
	POST: Plots its experimental cumulative distribution
	"""
	resdic={}
	with open(resfile,"rb") as r:
		for line in r:
			if "======"==line[:6]:
				cline=line
				resdic[cline]=[]
			else:
				resdic[cline].append(line)
	name_energies = []
	for el in resdic:
		if len(resdic[el])==3:
			name_energies.append((float(resdic[el][-1].strip().split(',')[-1]), el.strip()))
			
	name_energies_sorted = [(x,y) for x,y in sorted(name_energies)]
	
	print "BEST"
	for els in name_energies_sorted[:10]:
		print els
	print "WORST"
	for i, els in enumerate(name_energies_sorted[:]):
		print i, els
	
	with open("/home/macenrola/Documents/XYLENE/images/SP_DONT_WORK", "wb") as w:
		for i, els in enumerate(name_energies_sorted):
			w.write("{} {}\n".format((els[0]+0.4340774664)*627.5, i/float(len(name_energies_sorted))))
			if any([x in els[1] for x in ["_0.2_0.2_", "_2.0_0.8", "_2.0_0.5", "_1.0_0.8"]]):
				print "{} {} {} {}".format(i, els[1], (els[0]+0.4340774664)*627.5, i/float(len(name_energies_sorted)))
	
	print "ELS 500"
# =============================================================================
# 	for els in name_energies_sorted[501:510]:
# 		print "scp uccahcl@myriad.rc.ucl.ac.uk:/home/uccahcl/Scratch/XYLENE/pm6ts-stability-test/z-shift/{0} {0} &".format(els[1][6:-6])
# 	plt.plot([(x+0.4340774664)*627.5 for x,_ in name_energies_sorted[:]])   # TS HEIGHT
# =============================================================================
 	plt.xlabel("Index of MO@CB[7] starting point")
	plt.ylabel("$E_{TS}$ - $E_{m-xylene@CB[7]}$ [kcal/mol]")
# =============================================================================
# # =============================================================================
# # 		plt.plot(sorted([(float(x.split(",")[-1].strip())+0.6057797455-0.3393539069)*627.5 for x in r.readlines()])[:-58])   # BINDING OF TS
# # =============================================================================
# 		plt.plot(sorted([(float(x.split(",")[-1].strip())+0.4340774664)*627.5 for x in r.readlines()])[:-58])   # TS HEIGHT
# 		plt.xlabel("Index of MO@CB[7] starting point")
# 		plt.ylabel("$E_{bind}$ [kcal/mol]")
# 		
# 	with open(resfile,"rb") as r:
# 		lines = r.readlines()
# 		sp = [x.strip().split(',') for x in lines]
# 		names, energies = list(zip(*sp))
# 		names = [x for _,x in sorted(zip(energies, names))]
# 		
# 		print names, sorted(energies)
# =============================================================================
def plot_nice_FES_FOR_COVER(festxt, colvar_traj):
	"""
	PRE  : Takes in a txt file 
	POST : Will plot a nice FES plot in 3D for the cover picture
	"""
	def get_xyz(fname):
		"""
		pre: takes in a file name
		post: will plot the file
		"""
		with open(fname, 'rb') as r:
			lines = r.readlines()
	
		x = [float(k.strip().split()[0]) for k in lines]
		y = [float(k.strip().split()[1]) for k in lines]
		try:
			z = [float(k.strip().split()[2]) for k in lines]
			return x,y,z
		except:
			return x,y

	def get_z_from_interp(x, y, z_grid, x_traj, y_traj):
		"""
		PRE : Takes in the long vector forms x, y and z that correspond to the grid
		POST: Will use it to build a 2d interpolation function to generate the z_traj corresponding to the x_traj and y_traj
		"""

		f = interpolate.RectBivariateSpline(x[:,0], y[0,:], z_grid)
		z_traj = f.ev(x_traj, y_traj)
		return z_traj
	def plot_colormap(x_grid,y_grid,z_grid,name, xtraj, ytraj):
		"""
		pre: takes in x,y,z and will plot a colormap
		"""
		dim = 100
		x = np.asarray(x_grid).reshape((dim, dim))
		y = np.asarray(y_grid).reshape((dim, dim))
		z = np.asarray([k*627.5 for k in z_grid]).reshape((dim, dim))
		levels = plt.MaxNLocator(nbins=15).tick_values(z.min(), z.max())
		# get colormap
		ncolors = 256
		color_array = plt.get_cmap('gist_rainbow')(range(ncolors))
		
		# change alpha values
		color_array[:,-1] = np.linspace(1.0,0.0,ncolors)
		
		# create a colormap object
		map_object = mcolors.LinearSegmentedColormap.from_list(name='rainbow_alpha',colors=color_array)
		cmap = plt.get_cmap('Greys')
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		#cf = ax.contourf(x,y,z,levels=levels, cmap=cmap)
		
		ls = LightSource(0, 25)
		rgb = ls.shade(z, cmap=map_object, vert_exag=0.1, blend_mode='soft')
		cf = ax.plot_surface(x, y,z, rstride=1, cstride=1, facecolors=rgb, linewidth=0, antialiased=True, shade=True)
		start = 40000-15000
		dt = 2000
		print len(xtraj[start:start+dt])
		z_traj = get_z_from_interp(x, y, z, xtraj[start:start+dt], ytraj[start:start+dt])
		ax.plot(xtraj[start:start+dt], ytraj[start:start+dt], z_traj, linestyle=":", color="k")
		#plt.axhline(0.7,linestyle='--', linewidth=0.4, color='k', xmin=0, xmax=0.325)
		#plt.axvline(0.295,linestyle='--', linewidth=0.4, color='k', ymin=0.7)
		ax.annotate("", xy=(0.517, 0.6942), xytext=(0.517+0.001, 0.6942-0.0001), arrowprops=dict(arrowstyle="simple"))
		ax.set_xlabel(r"$C_1$")
		ax.set_ylabel(r"$C_2$")
		#cbar=fig.colorbar(cf, ax=ax)
		#cbar.ax.set_ylabel('Free energy (kcal/mol)')
		# ax.view_init(-126, -155) # for mountains
		ax.view_init(51, -112) # for valleys
		plt.axis('off')
		plt.grid(b=None)
		plt.show()
		
	x, y,z = get_xyz(festxt)
	_, co, ct = get_xyz(colvar_traj)
	
	plot_colormap(x, y, z, '-'.join(festxt.split('/')[-4:]), co, ct)

def plot_nice_FES(festxt, colvar_traj):
	"""
	PRE  : Takes in a txt file 
	POST : Will plot a nice FES plot
	"""
	def get_xyz(fname):
		"""
		pre: takes in a file name
		post: will plot the file
		"""
		with open(fname, 'rb') as r:
			lines = r.readlines()
	
		x = [float(k.strip().split()[0]) for k in lines]
		y = [float(k.strip().split()[1]) for k in lines]
		try:
			z = [float(k.strip().split()[2]) for k in lines]
			return x,y,z
		except:
			return x,y


	def plot_colormap(x,y,z,name, xtraj, ytraj):
		"""
		pre: takes in x,y,z and will plot a colormap
		"""
		dim = 100
		x = np.asarray(x).reshape((dim, dim))
		y = np.asarray(y).reshape((dim, dim))
		z = np.asarray([k*627.5 for k in z]).reshape((dim, dim))
		levels = plt.MaxNLocator(nbins=15).tick_values(z.min(), z.max())
		cmap = plt.get_cmap('PiYG')
		fig, ax = plt.subplots(constrained_layout=True)
		cf = ax.contourf(x,y,z,levels=levels, cmap=cmap)
		start = 40000-15000
		dt = 2000
		print len(xtraj[start:start+dt])
		ax.plot(xtraj[start:start+dt], ytraj[start:start+dt], linestyle=":", color="k")
		plt.axhline(0.7,linestyle='--', linewidth=0.4, color='k', xmin=0, xmax=0.325)
		plt.axvline(0.295,linestyle='--', linewidth=0.4, color='k', ymin=0.7)
		ax.annotate("", xy=(0.517, 0.6942), xytext=(0.517+0.001, 0.6942-0.0001), arrowprops=dict(arrowstyle="simple"))
		ax.set_xlabel(r"$C_1$")
		ax.set_ylabel(r"$C_2$")
		cbar=fig.colorbar(cf, ax=ax)
		cbar.ax.set_ylabel('Free energy (kcal/mol)')
		plt.show()
		
	x, y,z = get_xyz(festxt)
	_, co, ct = get_xyz(colvar_traj)
	
	plot_colormap(x, y, z, '-'.join(festxt.split('/')[-4:]), co, ct)
	
def plot_double_ring_e_pot(d=5.85, R=3.6):
	"""
	PRE : Takes in the radius of a double ring system with a charge located around the charged ring of radius R
	POST: Will plot the potential as a position of the charge around the two rings
	"""
	x = np.linspace(-5*d, 5*d, 200)
	epot = (1/((d-x)**2+R**2)**.5 + 1/(x**2+R**2)**.5)
	plt.plot(x, epot/max(epot))
	plt.xlabel("Distance between the charged ring and the point charge [Angstrom]")
	plt.ylabel("Electrostatic potential [a.u.]")
	
	
def get_energies_pm6_stability(summary_file):
	"""
	PRE: summary_file is obtained after running "for f in *out; do echo '==========' $f '==========='; grep RPM6D3 $f | tail -n5 | grep 'SCF Done:'; done > full_hessian_pm6stability "
	On a list of out files that contain a TS optimisation with calcall and a freq. between ==== === lines there will be only one line if the computation converged and several if it didn't huhu
	POST: Plots the distribution of the stability energies
	"""
	cline = None
	resdic = {}
	with open(summary_file, "rb") as r:
		for line in r:
			if line[0]=="=":
				cline=line
			elif cline not in resdic: #### didn't hit the === tag so it means we are in between lines
				resdic[cline] = line
			else: # means there are several lines between === lines
# =============================================================================
# 				print "{} already in the dic".format(cline)
# =============================================================================
				resdic[cline] = "NOT CONVERGED"
				continue
			
	energies = []
	with open(summary_file+"_tmp_check","wb") as w:
		with open(summary_file+"_failures_to_converge", "wb") as wf:
			for i, k in enumerate([y for _,y in sorted(zip(resdic.values(), resdic.keys()))]):
				if resdic[k]!="NOT CONVERGED":
					w.write("{}\t{}\n".format( k,  (float(resdic[k].strip().split()[-5])+0.4340774664)*627.5))
					print k,  (float(resdic[k].strip().split()[-5])+0.4340774664)*627.5
					energies.append(float( resdic[k].strip().split()[-5] ))
					print i
				else:
					wf.write(k.strip().strip("=").strip(" ") +"\n")


	with open(summary_file+"_numbers", "wb") as w:
		w.writelines(["{}\t{}\n".format(str(y[1]), str(y[0])) for y in zip(sorted([(x+0.4340774664)*627.5 for x in energies])[9:], np.linspace(0,1, len(energies)) )])
	plt.plot(sorted([(x+0.4340774664)*627.5 for x in energies])[9:])
	plt.ylabel("Energy of MO TS in CB7 vs protonated xylene in CB7 [kcal/mol]")
	plt.xlabel("index")
	
def compute_coordination_number_xyz(fname):
	"""
	PRE  : Takes in an xyz trajectory with an atom index and coordiation parameters 
	POST : Will return the coordination number as described in https://www.cp2k.org/exercises:2015_cecam_tutorial:mtd1
	"""
	r0, nn, nd = 1.5, 8., 14.
	atm_list, snapshot_array= return_list_of_snapshots_and_atom_list(fname)
	if len(atm_list) == 127:
		xylenelist = range(18)
		cblist = range(18, 127)
		atomindex = 18
	else:
		xylenelist = range(126, 144)
		atomindex = 144
		cblist = range(126)
# =============================================================================
# 	print [atm_list[x] for x in xylenelist]
# =============================================================================
	coords = []
	centroiddist = []
	for j, snap in enumerate(snapshot_array):
		coordnumber = 0
		for i in xylenelist:
			if i == atomindex:
				continue
			distsquared = sum([(x-y)**2 for x,y in zip(snap[i], snap[atomindex])])**.5
			coordnumber += (1-(distsquared/r0)**nn)/(1-(distsquared/r0)**nd)
		coords.append(coordnumber)
		xylcentroid, cbcentroid = [],[]
		for k, p in enumerate(snap):
			if k in xylenelist:
				xylcentroid.append(p)
			else:
				cbcentroid.append(p)
# =============================================================================
# 		print np.mean(np.array(xylcentroid), 0)
# 		print np.mean(cbcentroid, 0)
# =============================================================================
		centroiddist.append(sum([(x-y)**2 for x,y in zip(np.mean(xylcentroid, 0), np.mean(cbcentroid, 0))])**.5)
# =============================================================================
# 		print j, centroiddist[-1], coords[-1]
# =============================================================================

	try:
		out = [x>10 for x in centroiddist].index(True)
		n = [x>0.5 for x in coords[:out]].count(True)
		print out, coords[out], centroiddist[out], sum(coords[:out]), sum(coords), n, float(n)/out
		return sum(coords), float(n)/out
	except:
		return 0,0
	
def clearthedump(dumpfile, listofcrits, thresh):
	"""
	PRE  : Takes in a dump file and list of files with criteria
	POST : Will create a dump file without the files that don't pass the criteria
	"""
	critdic = {}
	for f in listofcrits:
		with open(f, "rb") as r:
			for line in r:
				parts = line.strip().split()
				if float(parts[2])>thresh:
					print "skipping", line
					continue
				else:
					fname = parts[0].split("/")
					fname = "{}-{}".format(fname[-2], fname[-1][:-14])
					critdic[fname] = parts
					
	with open(dumpfile+"-curated-dump", "wb"): pass
		
	with open(dumpfile, "rb") as r:
		for line in r:
			parts = line.split()[1]
			if parts[:32] in critdic:
# =============================================================================
# 			if float(line.strip().split()[-1])>1e-4:
# =============================================================================
				with open(dumpfile+"-curated-dump", "ab") as a: 
# =============================================================================
# 					a.write("{}\t{}\t{}".format( critdic[parts[:32]][1], critdic[parts[:32]][2], line))
# =============================================================================
					a.write(line)
			else: 
				continue
	
	
def make_average_energies(guest_list, ref):
	"""
	PRE  : Takes in a list of guest and a reference binding site, here cb8-mv+ 
	POST : Will print the average energy after the first 10 ps of equlibration 
	"""
	with open(ref, "rb") as r:
		lines = r.readlines()
	EREF = [float(x.strip().split()[4]) for x in lines[1000:]]
	
 	for f in guest_list:
		name = f.split("/")[-1].split(".")[0]
		try:
			comp_f = glob.glob(f.split(".")[0].replace(name, "MV1_CB8_{}*".format(name[1:]))+"*")[0]
		except:
			try:
				comp_f = glob.glob(f.split(".")[0].replace(name, "MV1_CB8_{}*".format(name[:]))+"*")[0]
			except:
				try:
					comp_f = glob.glob(f.split(".")[0].replace(name, "MV1_CB8_{}*".format(name[3:5]+name[:3]+name[-4:]))+"*")[0]
				except:
					print f, "not working"
					continue
		print f
		with open(f, "rb") as r:
			lines = r.readlines()
			EMOL = [float(x.strip().split()[4]) for x in lines[1000:]]

		print comp_f
		with open(comp_f, "rb") as r:
			lines = r.readlines()
			ECOMP = [float(x.strip().split()[4]) for x in lines[1000:]]		
		
		print "ENERGY DIFFERENCE:", (sum(ECOMP)/len(ECOMP)-sum(EMOL)/len(EMOL)-sum(EREF)/len(EREF))*627.5
		
def get_atom_number_in_traj_file(xyz_pos):
	with open(xyz_pos, "rb") as r:
		for line in r:
			return int(line)
		
def plot_summary_EKIN(fcomplex, gcomplex, smallfrag, bigfrag):
	"""
	PRE : Takes in EKIN files that contain a list of kinetic energies in kcal/mol
	POST: Plot them along with their variance
	"""
	for f in zip([fcomplex, gcomplex, smallfrag, bigfrag], ["E(G2@CB[7])","E(G2) vacuum","E(G2) @CB[7]","G2@E(CB[7])"], [.2,.2,0.2,0.2]):
		with open(f[0], "rb") as r:
			width =1000
			ekin = cPickle.load(r)
			moving_av = np.convolve(ekin, np.ones(width)/width)
			plt.plot(moving_av[10002:-width], label="{0} ({1:3.1f}/{2:3.1f})".format(f[1], np.mean(ekin[10002:]), np.var(ekin[10002:])), linewidth=1, alpha=1)
			plt.plot(ekin[10002:], linewidth=f[2], alpha=f[2])
	plt.title(gcomplex.split("/")[-1])
	plt.legend()
# =============================================================================
# 	plt.show()
# =============================================================================
	plt.savefig(gcomplex+"-curves.png")
	
def config_entropy(omega):
	"""
	PRE   : takes in a frequency 
	POST  : Returns the configurational entropy linked to that frequency
	"""
	R, k,T,hbar,h = 8.314, 1.3807e-23, 300, 1.054e-34,1.054e-34*2*np.pi  # those are the masses of  J / molK. kB T and hbar
	THETA = h*omega/k/T
	S = R*((THETA)/(np.exp(THETA)-1) -
				np.log(1-np.exp(-THETA)))
	
	return S
	
def get_entropy_from_xyz_file(fnameXYZ, label='CB8'):
	"""
	PRE  :  Using the method from Schlitter, Jurgen, and Matthias Massarczyk. "Estimating configurational entropy and energy of molecular systems from computed spectral density." arXiv preprint arXiv:1909.04726 (2019).
	        This method will extract the velocity power spectrum out of a xyz position trajectory and integrate its density along with an entropy expression
	POST :  Will return Spectrally Resolved Estimation (SRE) for entropy for the given trajectory
	"""
	with open(fnameXYZ, 'rb') as r:
		for line in r:
			num = int(line.strip())
			break
	freqs, magnitude, correlation = power_spectrum_from_velocityxyz(fnameXYZ, label )
	correlation_half = correlation[len(correlation)/2:]
	wavenum_half, cumcorrelation = (freqs[len(correlation)/2:], np.cumsum(correlation[len(correlation)/2:]))
	freqs_half = [2.99e10*x for x in wavenum_half] # in omega
	correlation_half, wavenum_half, freqs_half = correlation, freqs, [2.99e10*x for x in freqs]
	S = [config_entropy(x) for x in freqs_half]
	
	D = np.trapz(correlation_half, wavenum_half)
# =============================================================================
# 	plt.plot(wavenum_half, correlation[len(correlation)/2:]/max(correlation))
# 	plt.plot(wavenum_half, S/max(S))
# 	plt.plot(wavenum_half, cumcorrelation/max(cumcorrelation))
# 	plt.plot(wavenum_half, np.cumsum([(3.0*num-6)/D *x*y for x,y in zip(correlation_half, S)]))
# 	plt.plot(wavenum_half[1:], [(3*num-6)/D*x*y for x,y in zip(correlation_half[1:], S[1:])])
# 	plt.show()
# =============================================================================
	print  np.trapz([(3.0*num-6)/D *x*y for x,y in zip(correlation_half[1:], S[1:])], wavenum_half[1:])

def compare_CB8_ternary_energies(fname):
	"""
	PRE : Will take a cp2k energy file and assume that there is another one with a similar name
	POST: Will print the energy difference
	"""
	with open(fname, 'rb') as r:
		lines = r.readlines()
		MV2_eners = [float(x.strip().split()[4]) for x in lines[1:]]
		
	with open(fname.replace("MV2_", "/CHARGEONE/"), 'rb') as r:
		lines = r.readlines()
		MV1_eners = [float(x.strip().split()[4]) for x in lines[1:]]
		print "{}:\n{}".format(fname, (sum(MV2_eners)/len(MV2_eners) - sum(MV1_eners)/len(MV1_eners))*627.5)
		
def get_energy_from_XYZFile(xyzfile):
	"""
	PRE : takes an xyz file
	POST: Will extract the energy at every step and compute its average
	"""
	start = 450
	with open(xyzfile, "r") as r:
		lines = r.readlines()
	energies = [float(x.strip().split()[-1]) for x in lines if 'E =' in x ]
	plt.plot(energies)
	plt.savefig(xyzfile+"_ENERGIES.pdf")
	plt.close()
	average = sum(energies[start:])/len(energies[start:])	
	print average # average energy in hartree
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
	import scipy.signal 
	from scipy import interpolate
	import glob
	import cPickle
	from matplotlib import pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib.colors import LightSource
	import matplotlib.colors as mcolors
	import random
	import operator
	import math
	import sys
# =============================================================================
# 	# COMPUTE VOLUME
# 	mol = Chem.MolFromPDBFile('/home/macenrola/Desktop/xylene_polish/volumes_structures/MO-vac.pdb', removeHs=False)
# 	print AllChem.ComputeMolVolume(mol, boxMargin=1.0)
# =============================================================================
 	#get_entropy_from_xyz_file('/home/macenrola/Documents/CB8-electrochemistry/SPECTRA_CBs/CB5_neutral.xyz-pos-1.xyz', 'CB5')
	#get_entropy_from_xyz_file('/home/macenrola/Documents/CB8-electrochemistry/SPECTRA_CBs/CB6-neutral.com_OUT.out.xyz-pos-1.xyz', 'CB6')
 	#get_entropy_from_xyz_file('/home/macenrola/Documents/CB8-electrochemistry/SPECTRA_CBs/CB7-neutral.com_OUT.out.xyz-pos-1.xyz', 'CB7')
 	#get_entropy_from_xyz_file('/home/macenrola/Documents/CB8-electrochemistry/SPECTRA_CBs/CB8-neutral.xyz-pos-1.xyz', 'CB8')

# =============================================================================
# 	for f in sorted(glob.glob('/media/macenrola/41e645a2-504d-4468-adcc-106c8c763adb/CB8/archivedwl-53/trial_spectrum_tech/vacuum_1+/*xyz-pos-1.xyz')):
# 		print f 
# 		get_frequencies_quasiharmonic(f)
# 		get_entropy_from_xyz_file(f)
# =============================================================================
	for f in glob.glob('/media/macenrola/41e645a2-504d-4468-adcc-106c8c763adb/CB8/MV1_CB8*ALIGNED.xyz'):
		print f
		#get_entropy_from_xyz_file(f)
		get_energy_from_XYZFile(f)
		print f.replace('MV1', 'MV2_MV1')
		#get_entropy_from_xyz_file( f.replace('MV1', 'MV2_MV1'))
		get_energy_from_XYZFile( f.replace('MV1', 'MV2_MV1'))

# =============================================================================
# 	#RAMAN SHIFT LOOOOL
# 	power_spectrum_from_velocityxyz('/home/macenrola/Documents/CB8-electrochemistry/SPECTRA_CBs/CB5_neutral.xyz-pos-1.xyz', 'CB5')
# 	power_spectrum_from_velocityxyz('/home/macenrola/Documents/CB8-electrochemistry/SPECTRA_CBs/CB6-neutral.com_OUT.out.xyz-pos-1.xyz', 'CB6')
# 	power_spectrum_from_velocityxyz('/home/macenrola/Documents/CB8-electrochemistry/SPECTRA_CBs/CB7-neutral.com_OUT.out.xyz-pos-1.xyz', 'CB7')
# 	power_spectrum_from_velocityxyz('/home/macenrola/Documents/CB8-electrochemistry/SPECTRA_CBs/CB8-neutral.xyz-pos-1.xyz', 'CB8')
# 	#END RAMAN SHIFT
# =============================================================================
# =============================================================================
# 	#### DO THE EKIN FLOATING
# # =============================================================================
# # 	split_velocity_file("/home/macenrola/Documents/Thermal_energy_collector/OUTS/RANK_G46_PBNUM_91240295_COMPLEX.sdf.xyz-pos-1.xyz", 43)
# # =============================================================================
# 	if len(sys.argv) <3:
# 		print """Provide more argments (3) please. Use on GUEST FILES ONLY
# 		Use make_avg_coords.py slice fname or 
# 		make_avg_coords.py kin fname"""
# # =============================================================================
# # 		plot_E_kin("/home/macenrola/Documents/Thermal_energy_collector/OUTS/RANK_G46_PBNUM_91240295_COMPLEX.sdf.xyz-pos-1.xyz-large_frag.xyz-EKIN", "lol", "--")
# # =============================================================================
# 	elif len(sys.argv) == 3 and sys.argv[1] == "slice":
# 		print "will read the guest and provide the number to slice the complex then slice the complex"
# 		guest_name = sys.argv[2]
# 		n = get_atom_number_in_traj_file(guest_name)
# 		complex_name = guest_name.replace("GUEST", "COMPLEX")
# 		split_velocity_file(complex_name, n)
# 	elif len(sys.argv) == 3 and sys.argv[1] == "kin":
# 		print "will read the guest and compute the kinetic energy for the guest, the fragments and the complex"
# 		guest_name = sys.argv[2]
# 		complex_name = guest_name.replace("GUEST", "COMPLEX")
# 		small_frag_name = complex_name + "-small_frag.xyz" 
# 		big_frag_name = complex_name + "-large_frag.xyz"
# 		for k in [guest_name,complex_name,small_frag_name,big_frag_name]:
# 			get_kinetic_energy_from_position_file(k)
# 	elif len(sys.ar.gv) == 3 and sys.argv[1] == "plot":
# 		print "will read the guest and compute the kinetic energy for the guest, the fragments and the complex"
# 		guest_name = sys.argv[2]
# 		complex_name = guest_name.replace("GUEST", "COMPLEX")
# 		small_frag_name = complex_name.replace("-EKIN", "-small_frag.xyz-EKIN" )
# 		big_frag_name = complex_name.replace("-EKIN","-large_frag.xyz-EKIN")
# 		plot_summary_EKIN(complex_name, guest_name, small_frag_name, big_frag_name)
# 		
# 	### END OF EKIN FLOATING
# =============================================================================
	
	
# =============================================================================
# 	get_frequencies_quasiharmonic("/home/macenrola/Documents/CB8-electrochemistry/get_entropy_quasiharmonic/Benzene_sol.out.xyz-pos-1.xyz")
# =============================================================================
# =============================================================================
# 	plt.plot(*power_spectrum_from_velocityxyz("/home/macenrola/Documents/CB8-electrochemistry/get_entropy_quasiharmonic/trial_spectrum_tech/vacuum_2+/MV1_CB8_sol.com_OUT.out.xyz-pos-1.xyz"))
# =============================================================================
# =============================================================================
# ASSIGNING THE PAIRS OF COMPLEXES
# 	print "GETTING THE G2"
# 	flist = glob.glob("/home/macenrola/Documents/CB8-electrochemistry/ENERS/tolaunch_vac*/*.ener")
# 	for f in sorted(flist):
# 		if "_CB8_" in f:
# 			print "removing", f
# 			flist.remove(f)
# 
# 	print "PRINTING G2"
# 	for f in sorted(flist):
# 		print f
# 	
# 	print "PAIRING G2 and COMPLEX G2"
# 	make_average_energies(flist,
# 					   "/home/macenrola/Documents/CB8-electrochemistry/ENERS/tolaunch_vacuum/MV1_CB8_sol.com_OUT.out.xyz-1.ener")
# 	
# =============================================================================
	
# =============================================================================
# 	compute_coordination_number_xyz("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/XYZ_sanity/300/101-mxylene-prot-cb7-H-single.inp-pos-1.xyz")
# =============================================================================
# =============================================================================
# ## COMPUTE THE COORDINATIONS
# 	record = "/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/600-CB-RES"
# 	with open(record, "wb"): pass
# 	for f in sorted(glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/XYZ_sanity/600/*-cb*-H-single.inp-pos-1.xyz")):
# 		print f
# 		res = compute_coordination_number_xyz(f)
# 		with open(record, 'ab') as a:
# 			a.write("{}\t{}\t{}\n".format(f, res[0], res[1]))
# =============================================================================
# =============================================================================
# # CLEARING THE DUMP
# 	clearthedump("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/DUMP_MTS-dump", glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/*00-CB-RES"), 0.2)
# 
# =============================================================================
# =============================================================================
# 	with open("/home/macenrola/Documents/H-S-compensation/PUBCHEM_ONLY_C_only_19_atoms", "rb") as r:
# 		lines = r.readlines()
# 		elems = random.sample(lines, k =230)
# 	with open("/home/macenrola/Documents/H-S-compensation/C19_sample_200.can", "wb") as w:
# 		w.writelines(elems) 
# =============================================================================
	# process_z_matrix_trajectory('cb6.inp-pos-1-aligned.gzmat')
	# for f in ['/home/macenrola/Documents/XYLENE/inputs/for-reaction-frozen-cb/MO-CB6.inp-pos-1-aligned-just-CB6.xyz',
	# 			'/home/macenrola/Documents/XYLENE/inputs/for-reaction-frozen-cb/MO-CB7.inp-pos-1-aligned-just-CB6.xyz',
	# 			'/home/macenrola/Documents/XYLENE/inputs/for-reaction-frozen-cb/MP-CB6.inp-pos-1-aligned-just-CB6.xyz',
	# 			'/home/macenrola/Documents/XYLENE/inputs/for-reaction-frozen-cb/MP-CB7.inp-pos-1-aligned-just-CB6.xyz']:
	# process_xyz_coordinates('/home/macenrola/Documents/XYLENE/base-systems-equilibrated/starting-trajectories-for-frozen/CB6-aligned-from-MP-centered.xyz')
	# for f in glob.glob("/home/macenrola/Documents/XYLENE/base-systems-equilibrated/equilibrated+NVE-long/just*aligned-centered.xyz"):
	# 	make_covariance_matrix(f)
# =============================================================================
# 	power_spectrum_from_velocityxyz("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/sample_vel-large-frag.xyz")
# =============================================================================
	# f="/home/macenrola/Documents/heptylamine/TRAJS/300-heptylamine.inp-pos-1.xyz.sdf-w-charge.sdf-aligned.sdf"
	#for f in glob.glob("/home/macenrola/Documents/heptylamine/TRAJS/SDFs/*00-heptylamine.inp-pos-1.xyz.sdf-w-charge.sdf-aligned.sdf-rzlists"):
	#	print f
		# edit_xyz_file_to_add_format_charge(f)
		# align_trajectory_with_axes_according_to_cb(f+"-w-charge.sdf")
	#	compute_cylindrical_transform_for_sdf_trajectory(f)
# =============================================================================
# 	make_covariance_matrix("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/sample_pos_just_cb.xyz")
# =============================================================================
# =============================================================================
# 	for f in glob.glob('/home/macenrola/Documents/XYLENE/base-systems-equilibrated/equilibrated+NVE-long/M*.inp-pos-1.xyz-aligned-centered.xyz'):
# 		make_covariance_matrix(f)
# =============================================================================
# =============================================================================
# 	for f in glob.glob('/home/macenrola/Documents/XYLENE/base-systems-equilibrated/equilibrated+NVE-long/just-*-prot.inp-pos-1.xyz-aligned-centered.xyz'):
# 		make_covariance_matrix(f)
# =============================================================================
# =============================================================================
# 	convert_position_snapshots_to_velocity_snapshots('/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/guest-host-sep/just_xylene.xyz')
# 	convert_position_snapshots_to_velocity_snapshots('/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/guest-host-sep/sample_pos.xyz')
# 	convert_position_snapshots_to_velocity_snapshots('/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/guest-host-sep/just_cb.xyz')
# 
# =============================================================================
# =============================================================================
# 	convert_position_snapshots_to_velocity_snapshots('/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/sample_pos_coupling.xyz')
# =============================================================================
# =============================================================================
# 	split_velocity_file("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/split_by_chunks/sample_vel_coupling.xyz")
# =============================================================================
# =============================================================================
# 	split_xyz_file_in_chuncks("/home/macenrola/Documents/Thesis/XYLENE/coupling/frequency_move/sample_vel_coupling.xyz")
# =============================================================================
	
	
# =============================================================================
# 	magnitude=[]
# 	for f in sorted(glob.glob("/home/macenrola/Documents/Thesis/XYLENE/coupling/frequency_move/sample_vel_coupling.xyz-part-*.xyz"))[:]:
# 		print f 
# 		mag = power_spectrum_from_velocityxyz(f)
# 		magnitude.append(mag)
# 	cPickle.dump(magnitude, open("/home/macenrola/Documents/Thesis/XYLENE/coupling/frequency_move/sample_vel_coupling.xyz-MAG", "wb"))
# 	
# 		
# =============================================================================
# =============================================================================
# 	plot_sequence_of_ffts("/home/macenrola/Documents/Thesis/XYLENE/coupling/frequency_move/sample_vel_coupling.xyz-MAG")
# =============================================================================
# =============================================================================
# 	get_kinetic_energy_from_velocity_file("/home/macenrola/Documents/Thesis/XYLENE/coupling/sample_vel_coupling.xyz")
# =============================================================================
# =============================================================================
# 	split_velocity_file_bis("/home/macenrola/Documents/Thesis/XYLENE/coupling/sample_vel_coupling.xyz")
# =============================================================================
# =============================================================================
# 	for i in range(10):
# =============================================================================
# =============================================================================
# 	plot_single_spectrum("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/ALL_100k/sample_vel_coupling.xyz-MAG", i=5)
# =============================================================================
# =============================================================================
# 	plot_double_ring_e_pot()
# =============================================================================
	
# =============================================================================
# 	
# 	## PLOT NICE FES
# 	plot_nice_FES("/home/macenrola/Documents/Thesis/XYLENE/manuscript_update/plotting_new/1-MP-CB6.inp-1_45000.restart.out", "/home/macenrola/Documents/Thesis/XYLENE/manuscript_update/plotting_new/1-MP-CB6.inp-COLVAR.metadynLog")
# =============================================================================
	
# =============================================================================
# 	## PLOT NICE FES
# 	plot_nice_FES_FOR_COVER("/home/macenrola/Documents/Thesis/XYLENE/manuscript_update/plotting_new/1-MP-CB6.inp-1_45000.restart.out", "/home/macenrola/Documents/Thesis/XYLENE/manuscript_update/plotting_new/1-MP-CB6.inp-COLVAR.metadynLog")
# 
# =============================================================================
# =============================================================================
# 	# ECDF FOR SP DONT WORK
# 	plot_sp_not_stable_output()
# =============================================================================
	
# =============================================================================
# 	plot_ener(4, "/home/macenrola/Documents/XYLENE/base-systems-equilibrated/equilibrated+NVE-long/just-mxylene-prot.inp-1.ener")
# =============================================================================
# =============================================================================
# 	
# 	# PLOT THE REACTIVE TRAJECTORY
# 	
# 	plot_instantaneous_binding_energy("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/breaking-down-cb-and-guest/complexres",
# 								   "/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/breaking-down-cb-and-guest/res/resfilegh2",
# 								   "/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/breaking-down-cb-and-guest/res/resfilegh1")
# 
# # =============================================================================
# =============================================================================
# 	# KS STACKED POPPING OLD
# 	stack_ks_plots([
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/popping/DUMP_MTS-dump-oxylene-cb6-*KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/popping/DUMP_MTS-dump-oxylene-cb7-*KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/popping/DUMP_MTS-dump-mxylene-cb6-*KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/popping/DUMP_MTS-dump-mxylene-cb7-*KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/popping/DUMP_MTS-dump-pxylene-cb6-*KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/popping/DUMP_MTS-dump-pxylene-cb7-*KS"),
# 			])
# 	
# =============================================================================
# =============================================================================
# 
# 	# KS STACKED POPPING NEW BUT UGLY
# 	stack_ks_plots([import matplotlib.colors as mcolors

# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/DUMP_MTS-dump-curated-dump-oxylene-prot-cb6-*.0-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/DUMP_MTS-dump-curated-dump-oxylene-prot-cb7-*.0-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/DUMP_MTS-dump-curated-dump-mxylene-prot-cb6-*.0-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/DUMP_MTS-dump-curated-dump-mxylene-prot-cb7-*.0-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/DUMP_MTS-dump-curated-dump-pxylene-prot-cb6-*.0-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/DUMP_MTS-dump-curated-dump-pxylene-prot-cb7-*.0-KS"),
# 			])
# 	
# =============================================================================
	
# =============================================================================
# 	# POPPING RATE
# 	make_mtd_popping_rate_plot(glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_2_TARGETS/popping_times")[0])
# 	#make_mtd_popping_rate_plot(glob.glob("/home/macenrola/Documents/XYLENE/inputs/popping-from-cb/SLOW_POPPING_W_TARGET/OUTS/popping_times")[0])
# =============================================================================
# =============================================================================
# 	# KS STACKED _ISOMERISATION
# 	stack_ks_plots([
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/DUMP_MTS_CLEANED_FES-dump-MO-vac-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/DUMP_MTS_CLEANED_FES-dump-MP-vac-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/DUMP_MTS_CLEANED_FES-dump-MO-CB6-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/DUMP_MTS_CLEANED_FES-dump-MP-CB6-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/DUMP_MTS_CLEANED_FES-dump-MO-CB7-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/outputs/DUMP_MTS_CLEANED_FES-dump-MP-CB7-*-KS")
# 			])
# 	
# 	
# =============================================================================
# =============================================================================
# 	# SLOW KS STACKED _ISOMERISATION
# 	stack_ks_plots([
# 			glob.glob("/home/macenrola/Documents/XYLENE/images/DUMP_ISO/DUMP_ALL-dump-MP-vac-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/images/DUMP_ISO/DUMP_ALL-dump-MP-CB6-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/images/DUMP_ISO/DUMP_ALL-dump-MP-CB7-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/images/DUMP_ISO/DUMP_ALL-dump-MO-vac-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/images/DUMP_ISO/DUMP_ALL-dump-MO-CB6-*-KS"),
# 			glob.glob("/home/macenrola/Documents/XYLENE/images/DUMP_ISO/DUMP_ALL-dump-MO-CB7-*-KS")
# 			])
# 	
# 	
# =============================================================================
# =============================================================================
# 	
# #MAKE MTD TIME PLOT
# 	make_mtd_time_plot("/home/macenrola/Desktop/xylene_polish\/slow_summary_file_without_low_p")
# 	
# =============================================================================
	
# =============================================================================
	# PLOT STACKED SPECTRA
# =============================================================================
#  	#plot_three_stacked_spectra("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/ALL_100k/sample_vel_coupling.xyz-MAG", [0,1,2,5,9], "/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/ALL_FREQS/sample_vel_coupling.xyz-MAG")
#  	plot_three_stacked_spectra("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/ALL_100k/sample_vel_coupling.xyz-MAG", [0,1,2,5], "/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/ALL_FREQS/sample_vel_coupling.xyz-MAG")
# 
# =============================================================================
# =============================================================================
# =============================================================================
# # PLOT EKIN
# 	plot_E_kin("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/sample_vel_coupling.xyz-small_frag.xyz-EKIN",  r'$m$-xylene', "-.")
# 	plot_E_kin("/home/macenrola/Documents/XYLENE/correlation_for_reaction/slow-reaction-MP-CB6/vibrational_analysis/traj_from_mode_368/sample_vel_coupling.xyz-large_frag.xyz-EKIN", r'CB[6]', ':')
# 
# 
# =============================================================================


# =============================================================================
# 	#ENERGY PM6 STABILITY
# 	get_energies_pm6_stability("/home/macenrola/Documents/XYLENE/inputs/SP-DONT-WORK/z_shift_full_hessian/full_hessian_pm6stabilityfull_hessian_pm6stability_MO_in_CB7")
# 
# 
# 
# =============================================================================
