
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

def parse_xyz_block(xyz_block):
	"""
	PRE: Takes in an xyz block as a list of lines
	POST: Will return an array of array as np.array(np.array(), np.array(),...) each sub array contains [x,y,z] coord of a specific atom
	"""
	# let go of the atom names, the order is supposed to keep that information
	coords = np.array([np.array([float(y) for y in x.strip().split()[1:]]) for x in xyz_block[2:]], dtype=np.float64) # skip the two first lines that hold atom number and energy info
	return coords

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

def power_spectrum_from_velocityxyz(xyzfile,l):
	"""
	PRE : Takes in an xyz file containing velocity (most likely in bohr/au_time)
	A lot of confusion: the FT of the autocorrelation function of a signal X is the power spectrum of X 
	The power spectrum of X is the *square* of the amplitude of the FT of X
	Then everything makes sense
	IDEALLY ONE SHOULD GET THE VELOCITY DIRECTLY FROM THE TRAJECTORY TO AVOID JUMPS IF COMPUTED FROM THE CARTESIAN COORDINATES. RAW VELOCITIES WOULD BE GREAT
	POST: Will produce the power spectrum
	"""
	#print "Getting the snaps"
	atm_list, snapshot_array= return_list_of_snapshots_and_atom_list(xyzfile)
	#
	
	#print "Got the snaps, now reshaping"
	snapshot_array = [np.reshape(x, len(atm_list)*3) for x in list(snapshot_array)][450:] #-25000:
	snapshot_array = [x*1e-10 for x in snapshot_array]
	
	#### CONVERT SNAPS OF POSITION TO SPEED
	snapshot_speed = np.gradient(snapshot_array, 1/ts, edge_order=2, axis=0)
	#snapshot_array = snapshot_speed
	mass_matrix, sqrt_mass_matrix= make_mass_matrix_from_atom_list(atm_list)
	#### END CONVERSION
	
	#print "reshaping level 1"
	snapshot_array=np.array(snapshot_array).T
	zero_mean_array = []
	for el in snapshot_array:
		el = el-np.mean(el)
		zero_mean_array.append(np.array(el))
	snapshot_array = np.array(zero_mean_array)
	print snapshot_array.shape
	#print "Reshaping level 2"
	total_correlation = None
	#print "Computing the correlation"
#################  DOING THE FT of the autocorrelation function of X
# =============================================================================
# 
# 	for i, component in enumerate(snapshot_array): #
# 		print component.shape, np.diag(sqrt_mass_matrix)[i], i, ''.join(([''.join([x]*3) for x in atm_list]))[i]
# 		component = [x*np.diag(sqrt_mass_matrix)[i] for x in component] 
# 		n=len(component) #n=10000
# 		temp = scipy.signal.correlate(component, component, mode='same')
# 		#temp = np.array([np.abs(A)**2 for A in numpy.fft.fftshift(numpy.fft.fft(np.array(component)))])
# 		if i==0: total_correlation = temp
# 		else: total_correlation = total_correlation + temp
# 	#print "Summing the correlation"	
# 	magnitude  =[np.abs(A) for A in numpy.fft.fft(total_correlation)][:n/2]
# 	angle      =[np.angle(A) for A in numpy.fft.fftshift(numpy.fft.fft(total_correlation))]	
# 	freqs = np.linspace(0, 1/ts/2.99e10/2, n/2)
# =============================================================================
###################
############### DOING THE square of the FT of X
	for i, component in enumerate(snapshot_array): # 
		#print component.shape, np.diag(sqrt_mass_matrix)[i], i, ''.join(([''.join([x]*3) for x in atm_list]))[i]
		component = [x*np.diag(sqrt_mass_matrix)[i] for x in component] 
		#stdev = statistics.stdev(component)
		#component = [x for x in component if np.abs(x)<stdev]
# =============================================================================
# 		if i==0:
# 			print statistics.stdev(component)
# 			plt.plot(component)
# 			plt.show()
# =============================================================================
		temp = np.array([np.abs(A)**2 for A in (numpy.fft.fft(np.array(component)))][:])
		if i==0: magnitude = temp
		else: magnitude = magnitude + temp
	n=len(magnitude)	
	magnitude = magnitude[:n/2]
	freqs = np.linspace(0, 1/ts/2.99e10/2, n/2)	
################
############### DOING THE square of the FT of X through the welch method directly to smoothe outliers
# =============================================================================
# 	nw = 15000
# 	for i, component in enumerate(snapshot_array): # 
# 		#print component.shape, np.diag(sqrt_mass_matrix)[i], i, ''.join(([''.join([x]*3) for x in atm_list]))[i]
# 		component = [x*np.diag(sqrt_mass_matrix)[i] for x in component]
# 		n=len(component) #n=10000
# 		freqs, temp = scipy.signal.welch(component, fs=1/ts, nperseg=nw, return_onesided=True, noverlap=3000)
# 		if i==0: magnitude = temp
# 		else: magnitude = magnitude + temp
# 	freqs = [x/2.99e10 for x in freqs]
# =============================================================================
################	
	
# =============================================================================
# 	plt.plot(freqs, magnitude/max(magnitude), label=l) # gives the spectrum with x axis in cm-1
# 	plt.legend()
# 	plt.show()
# =============================================================================
	return  freqs, [x**2*y for x,y in zip(freqs, magnitude)], total_correlation 


	
	
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
	wavenum, magnitude, correlation = power_spectrum_from_velocityxyz(fnameXYZ, label )
# =============================================================================
# 	tosave = wavenum, magnitude, correlation
# 	with open(fnameXYZ+"_temp", "wb") as w: 
# 		cPickle.dump(tosave, w)
# 	with open(fnameXYZ+"_temp", "rb") as r: 
# 		wavenum, magnitude, correlation = cPickle.load(r)
# =============================================================================
	####### Here we apply a ad-hoc rescaling to avoid noise:
	toclean = [index for index, els in enumerate(wavenum) if els > 4000]
#	print toclean, wavenum[toclean[0]]
	
	magnitude[toclean[0]:] = [0]*len(toclean)
	magnitude[:toclean[0]] = [x-magnitude[toclean[0]-1] for x in magnitude[:toclean[0]]]
	magnitude = [x if x>0 else 0 for x in magnitude]
	
# =============================================================================
# 	plt.plot(wavenum, magnitude/max(magnitude), label='cb8 complex') # gives the spectrum with x axis in cm-1
# 	plt.legend()
# 	plt.show()
# =============================================================================
	
# =============================================================================
# 	correlation_half = correlation[len(correlation)/2:]
# =============================================================================
	#magnitude_half = magnitude[len(magnitude)/2:]
	#wavenum_half = freqs[len(magnitude)/2:]
# =============================================================================
# 	wavenum_half, cumcorrelation = (freqs[:len(correlation)/2], np.cumsum(correlation_half))
# =============================================================================
	cummag = np.cumsum(magnitude)

	freqs = [2.99e10*x for x in wavenum] # in omega
	#correlation_half, wavenum_half, freqs_half = correlation, freqs, [2.99e10*x for x in freqs]
	S = [config_entropy(x) for x in freqs]
	
	start=1
# =============================================================================
# 	D = np.trapz(correlation_half, wavenum_half)
# =============================================================================
	D = np.trapz(magnitude[start:], wavenum[start:])

# =============================================================================
# 	plt.plot(wavenum_half,correlation_half/max(correlation))
# 	#plt.plot(wavenum_half, S/max(S))
# 	#plt.plot(wavenum_half, cumcorrelation)
# 	#plt.plot(wavenum_half, np.cumsum([(3.0*num-6)/D *x*y for x,y in zip(correlation_half, S)]))
# 	#plt.plot(wavenum_half[1:], [(3*num-6)/D*x*y for x,y in zip(correlation_half[1:], S[1:])])
# =============================================================================
	plt.plot(wavenum,magnitude/max(magnitude))
	plt.plot(wavenum, S/max(S))
	plt.plot(wavenum, cummag/max(cummag))
	plt.title("Configurational Entropy: {} (J/mol K)".format(np.trapz([(3.0*num-6)/D *x*y for x,y in zip(magnitude[start:], S[start:])], wavenum[start:])))
	plt.savefig(fnameXYZ+'-IMAGE.pdf')
	plt.close()
	#plt.plot(wavenum_half, np.cumsum([(3.0*num-6)/D *x*y for x,y in zip(magnitude_half, S)]))
	#plt.plot(wavenum_half[1:], [(3*num-6)/D*x*y for x,y in zip(magnitude_half[1:], S[1:])])
	#plt.show()
	print  np.trapz([(3.0*num-6)/D *x*y for x,y in zip(magnitude[start:], S[start:])], wavenum[start:])


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
# 	ts = .5e-15	
# 	get_entropy_from_xyz_file("/home/macenrola/Documents/cb8/ANILINE-pos-1.xyz")
# =============================================================================
# =============================================================================
# 	ts =1e-15
# 	#get_entropy_from_xyz_file('/media/macenrola/41e645a2-504d-4468-adcc-106c8c763adb/CB8/archivedwl-53/trial_spectrum_tech/vacuum_1+/MV1_CB8_sesamol_sol.com_OUT.out.xyz-pos-1.xyz')
# 	get_entropy_from_xyz_file('/media/macenrola/41e645a2-504d-4468-adcc-106c8c763adb/CB8/MV1_CB8_sesamol_sol.com_OUT.out.xyz-pos-1_ALIGNED.xyz')
# 	
# =============================================================================
	ts =1e-15
	for f in glob.glob('/media/macenrola/41e645a2-504d-4468-adcc-106c8c763adb/CB8/MV1_CB8_*-ALIGNED.xyz')[:]:
		print f
		get_entropy_from_xyz_file(f)
		print f.replace('MV1', 'MV2_MV1')
		get_entropy_from_xyz_file( f.replace('MV1', 'MV2_MV1'))

