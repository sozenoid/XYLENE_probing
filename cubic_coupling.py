#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:51:59 2019

@author: macenrola
"""

def get_block_between_tags(fname, begin_tag, end_tag):
	"""
	PRE: Takes in a file and two tags, 
	POST: Will return the text between the first occurence of these two tags, between the first and the second occurence if the begin and end tag is the same
	"""
	with open(fname, "rt") as r: 
		lines=r.readlines()
	
	begin = lines.index(begin_tag)
	end = lines[begin+1:].index(end_tag)
	return lines[begin: begin+end+1]

def get_data_from_out(outfile):
	"""
	PRE  : The output comes from an input file with the route `#n PM6D3 Opt=(Restart) Freq=(Anharmonic, PrintDerivatives, InternalModes, cubic, HPModes)`
		The starting geometry should probably be minimised to at least Tight or very tight in Gaussian
	POST : Will return the cubic coefficients, the cartesian modes 
		Double check that the modes coefficients are given to 5 decimals
	"""
	tag_modes = " Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering\n"
	mode_block = get_block_between_tags(outfile, tag_modes, tag_modes)
# =============================================================================
# 	print "".join(mode_block)
# =============================================================================
	# Then gets the mode numbers and positions in the block based on the "A" sequencies in the block
	AA_indices = [i-1 for i, j in enumerate(mode_block) if "Frequencies ---" in j]
	spacing = AA_indices[1]-AA_indices[0]
	mode_dic={}
	for i in AA_indices:
# =============================================================================
# 		print "".join(mode_block[i-1:i-1+spacing])
# =============================================================================
		mode_num = mode_block[i-1].strip().split()
		meta = [x.strip().split() for x in mode_block[i+1:i+5]] # gets the frequency, reduced masses, force constants and IR intensities
		meta_vals = map(list, zip(*[x[x.index("---")+1:] for x in meta])) # gets the values only, by mode
		mode_vals = map(list, zip(*[x.strip().split() for x in mode_block[i+6:i+spacing-1]])) #gets the mode data along with the kind of atom
		for k, i in enumerate(mode_num):
			mode_dic[i] = (meta_vals[k], mode_vals[3+int(k)]) # adds the meta values and the mode values to the dictionary
		mode_dic["order"] = mode_vals[:3]
		
	cubic_dic = {}
	begin_cubic = " :        CUBIC FORCE CONSTANTS IN NORMAL MODES         :\n"
	end_cubic = " :       QUARTIC FORCE CONSTANTS IN NORMAL MODES        :\n"
	cubic_block = get_block_between_tags(outfile, begin_cubic, end_cubic)
# =============================================================================
# 	print "".join(cubic_block)
# =============================================================================
	for i, line in enumerate(cubic_block[9:-5]):
		parts = line.strip().split()
		cubic_dic["{}-{}-{}".format(*parts[:3])] = parts[3:6]
	#print "".join(cubic_block)

	return mode_dic, cubic_dic

def make_permutations(cub):
	"""
	PRE:Takes in a dictionary made from a gaussian output
	POST: Will fill up the dictionary with the permutation ijk > kij, kji, jik, jki, ikj, ijk
	"""
	for w in cub.keys():
		i,j,k = w.split("-")
		for els in itertools.permutations((i,j,k)):
			cub["{}-{}-{}".format(*els)] = cub[w]
	return cub

def procrustean_fit_for_xylene_part_of_eigenvectors(mod_xyl, mod_comp):
	"""
	PRE  : Provides the eigenvector as a dictionary as returned by get_data_from_out with the xylene part of the eigenvectors similarly ordered 
	POST : will return the correlation coefficient and the aligned version of the modes 
	"""
	
	# BUILD THE MATRICES FOR THE PROCUSTEAN TEST
	mat_xyl = []
	mat_comp = []
	
	for k in mod_xyl:
		print k, mod_xyl[k]

def get_mode_localisation_on_xylene(mode, mode_xyl_atoms):
	"""
	PRE: Takes in modes formatted as by get_data_from_out
	POST: Will print the sum of amplitudes corresponding to the xylene moiety, they are assumed to be located first in the list of atoms
	Observation: The xylene guest tends to have more movement than expected from its isolated form. It gets a larger share of kinetic energy from the NM distribution than when in the vacuum
	Conversely the CB gets less of that energy. The sum of amplitude of normed normal modes is 70 in the mxylene-h-single-CB6 system while 48 in vacuum
	Here are the results for energy localisation, fascinating
	mxylene-h-single-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	69.9628961769
	mxylene-h-single-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	66.8494180851
	mxylene-protonated-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	71.2384169403
	mxylene-protonated-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	68.9233001733
	mxylene-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	64.3905017355
	mxylene-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out
	62.4629103182
	"""
	xyl_localised_modes=[]
	for k in sorted(mode)[:-1]: # to avoid the order 
		LOC_ON_XYL_COEFF=sum([float(x)**2 for x in mod[k][1][:mode_xyl_atoms*3-1]])
		xyl_localised_modes.append((LOC_ON_XYL_COEFF, int(k), mod[k][0][0]))
		
	return sorted(xyl_localised_modes, key=lambda x: x[1])
	



def C_FUJISAKI(omega, h, kb, T):
	"""
	PRE   : Takes in frequency, h, kb and T
	POST  : Will return the parameter used in the FUJISAKI paper page
	Fujisaki, Hiroshi, Lintao Bu, and John E. Straub. "Vibrational energy relaxation (VER) of a CD stretching mode in cytochrome c." arXiv preprint q-bio/0403019 (2004).
	"""
	return 1/h/omega*(1-math.exp(-h*omega/kb/T))/(1+math.exp(-h*omega/kb/T))

def delta(omega, z, l=10):
	"""
	PRE: Takes in frequency in Hz with the lifetime given in cm-1
	POST: Returns the value of the thawed delta according to Fujisaki, Hiroshi, Lintao Bu, and John E. Straub. "Vibrational energy relaxation (VER) of a CD stretching mode in cytochrome c." arXiv preprint q-bio/0403019 (2004).
	"""
# =============================================================================
# 	if abs(omega/29979245800.)<1:
# 		print "resonance: diff is {} delta is {}".format(omega/29979245800., (1/math.pi * l/(l**2+(omega)**2)))
# 		return 0
# =============================================================================
	return 1/math.pi * l*z/((l*z)**2/4.+(omega)**2)
	

def compute_decay_rates_fujisaki(cubic_dic, mode_dic, mode_nbr, xyl_mode_list, CB_mode_list):
	"""
	PRE  : Takes a cubic dictionary for a complex, the mode dictionary for frequency information, the complex mode number to examine, the list of xylene assigned mode numbers, and the CB assigned mode numbers
	POST : Returns the decay rate for the mode to the xylene assigned group (2/2), the CB assigned group (2/2) or one of each.
	The formula is given by hbar pi/8/omega_alpha Sum_beta_gamma |Phi_alpha_beta_gamma|/omega_beta/omega_gamma (1+n_beta+n_gamma) *delta(omega_alpha - omega_beta - omega_gamma)
	The occupation number is given by (exp(hbar omega_alpha/k/T) -1) **(-1) and Phi_alpha_beta_gamma represents the cubic coupling coefficients
	The delta are "thawed" using the lorentzian approximation delta(x) = 1/pi * gamma/(gamma**2+x**2) with gamma being the lifetime and gamma= 3 cm-1
	the force constant present like this in the gaussian output
 ........................................................
 :        CUBIC FORCE CONSTANTS IN NORMAL MODES         :
 :                                                      :
 : FI =  Reduced values [cm-1]  (default input)         :
 : k  =  Cubic Force Const.[AttoJ*amu(-3/2)*Ang(-3)]    :
 : K  =  Cubic Force Const.[Hartree*amu(-3/2)*Bohr(-3)] :
 :......................................................:
	"""
	WXX, WCC, WCX, WXM, WCM, WMM = 0,0,0,0,0,0 # cumulative decay rate for xylene xylene denominated frequencies (WXX), CB and CB (WCC), CB and Xylene (WCX), Mixed and Mixed (WMM), Xylene and Mixed (WXM), CB and Mixed (WCM)
	cWXX, cWCC, cWCX, cWXM, cWCM, cWMM = 0,0,0,0,0,0 # associated counters
	mode_nbr = str(mode_nbr)
	for k in sorted(cubic_dic)[:]:

		cwdecay = 0
		nums = k.split('-')
		freqs = [float(mode_dic[x][0][0]) for x in nums] # Gets for mode num x the metadata (0) containing the frequency (0) in cm-1
		masses = [float(mode_dic[x][0][1]) for x in nums] #masses in amu
		if mode_nbr == nums[0]:
			fcmode = freqs[nums.index(mode_nbr)]
			nums.remove(mode_nbr)
			freqs.remove(fcmode)
			tfr = abs(float(cubic_dic[k][0])/(fcmode-freqs[0]-freqs[1]))
# =============================================================================
# 			if abs(fcmode-freqs[0]-freqs[1])/29979245800.<1.:
# 				print "resonance: diff is {} modes are {} {} {}, coupling coeff is :{}".format((fcmode-freqs[0]-freqs[1])/29979245800., nums[0], nums[1], mode_nbr, float(cubic_dic[k][1])*1e-18*6e23/4.18e3)
# =============================================================================
			#cwdecay = h/16./fcmode*(float(cubic_dic[k][1])*1e-18*amu**(-3/2.0)*1e-10**(-3))**2*masses[0]*masses[1]*masses[2]/freqs[0]/freqs[1]*(1+occ_num(freqs[0], h, kb,T)+ occ_num(freqs[1],h, kb,T))*delta(fcmode-freqs[0]-freqs[1])
			#cwdecay = h/16/fcmode.*(float(cubic_dic[k][1])*1e-18*amu**(-3./2.0)*(1e-10)**(-3))**2 /freqs[0]/freqs[1]*(1+occ_num(freqs[0], h, kb,T) + occ_num(freqs[1],h, kb,T))*delta(fcmode-freqs[0]-freqs[1],c)
			cwdecay = (float(cubic_dic[k][0]))**2 *(1+occ_num(freqs[0]) + occ_num(freqs[1]))*delta(fcmode-freqs[0]-freqs[1],1)
# =============================================================================
# 			cwdecay = C_FUJISAKI(fcmode, h, kb, T) * h**2/2 * (float(cubic_dic[k][1])*1e-18*amu**(-3/2.0)*1e-10**(-3))**2/freqs[0]/freqs[1] * (
# 									(1+occ_num(freqs[0], h, kb,T) + occ_num(freqs[1],h, kb,T) + 2*occ_num(freqs[1],h, kb,T)*occ_num(freqs[0],h, kb,T))*(delta(fcmode+freqs[0]+freqs[1])+delta(-fcmode+freqs[0]+freqs[1])) + (
# 											(occ_num(freqs[0], h, kb,T) + occ_num(freqs[1],h, kb,T) + 2*occ_num(freqs[1],h, kb,T)*occ_num(freqs[0],h, kb,T))*(delta(fcmode+freqs[0]-freqs[1])+delta(-fcmode+freqs[0]-freqs[1]))
# 											)
# 									)
# =============================================================================

			if all([x in xyl_mode_list for x in nums]):
				WXX+=cwdecay
				cWXX+=1
			elif all([x in CB_mode_list for x in nums]):
				WCC+=cwdecay
				cWCC+=1
			elif any([x in xyl_mode_list for x in nums]) and any([x in CB_mode_list for x in nums]): 
				WCX+=cwdecay
				cWCX+=1
			elif any([x in xyl_mode_list for x in nums]):
				WXM+=cwdecay
				cWXM+=1
			elif any([x in CB_mode_list for x in nums]):
				WCM+=cwdecay
				cWCM+=1
			else:
				WMM+=cwdecay
				cWMM+=1
		else: 
			continue
	print "Mode freq: ", fcmode/c, "overall: ", sum([WXX, WCC, WCX, WXM, WCM, WMM]) 
	print "counts: ", cWXX, cWCC, cWCX, cWXM, cWCM, cWMM
	#print "avg rate contribution: ",W1/cw1, W2/cw2, W3/cw3
	print "Overall rates: ", zip(['WXX', 'WCC', 'WCX', 'WXM', 'WCM', 'WMM'], [WXX, WCC, WCX, WXM, WCM, WMM])
	return [x for x in [WXX, WCC, WCX, WXM, WCM, WMM]]
	

def compute_decay_rates_fujisaki(cubic_dic, mode_dic, mode_nbr, xyl_mode_list, CB_mode_list):
	"""
	PRE  : 
	POST :
	"""

def get_lambda(mode_nbrs, mod, cub):
	"""
	PRE: Takes in a mode number, a cubic dictionary and mode dictionary and K and T
	POST: Returns the lambda_ijk value accordking to buttrin. It is assumed that the cm-1 format can yield an energy through *c*h which then can give lamdaijk as energy once divided by kb*T
	"""	
	cubic_fluctuation = 0
	mode_nbrs = [str(x) for x in mode_nbrs]
	for it in [str(x+1) for x in range(len(mod)-1)]:
		if it not in mode_nbrs: 
			try:
				cubic_fluctuation += occ_num(float(mod[it][0][0]))*(float(cub["{0}-{1}-{1}".format(it, mode_nbrs[0])][0]) - float(cub["{0}-{1}-{1}".format(it, mode_nbrs[1])][0]) -float(cub["{0}-{1}-{1}".format(it, mode_nbrs[2])][0]))**2
			except:
# =============================================================================
# 				print "{0}-{1}-{1}".format(it, mode_nbrs[0]), "{0}-{1}-{1}".format(it, mode_nbrs[1]), "{0}-{1}-{1}".format(it, mode_nbrs[2]), "not going through"
# =============================================================================
				pass
	return cubic_fluctuation/4./2/kbtcm # returns the value in J (/2/kb/T*(h*c)**2) or cm-1 (/2/kb/T*(h*c))

def get_lambda_array(mod, cub):
	"""
	PRE : takes the frequency and cubic information
	POST: Will build a three dimentional array of lambdas for future use, the numbering starts at 0 not 1 like the mode numbers
	"""
	dim = len(mod)-1
	lambda_array = np.zeros([dim, dim, dim])
	for i in range(dim):
		print i 
		for j in range(dim):
			for k in range(dim):
				lambda_array[i][j][k] = get_lambda([i+1,j+1,k+1], mod, cub)
	return lambda_array

def kijk(i,j,k, vijk, tauijk, lambdaijk, mod):
	"""
	PRE: Given the Tauijk 
	POST: Returns the kijk
	"""
	Vijktilde = vijk[i][j][k]/(1+(j==k))
	hbar = h
# =============================================================================
# 	res = (2*math.pi) /hbar * (Vijktilde*h*c)**2 *1./8/(4*math.pi *kb*T*h*c* lambdaijk)**.5 *1/(
# 			1+ (4*math.pi*(h*c*Vijktilde)**2*tauijk)/(8*hbar*h*c*lambdaijk))
# =============================================================================
# =============================================================================
# 	res = c* (Vijktilde)**2 *1./8/(4*math.pi * kbtcm* lambdaijk)**.5 *1/(
# 			1+ (c*4*math.pi*Vijktilde**2*tauijk)/(8*lambdaijk))
# =============================================================================
# =============================================================================
# 	res = c* (Vijktilde)**2 *1/(4* kbtcm* lambdaijk)**.5 *1/(
# 			1+ (c*4*Vijktilde**2*tauijk)/(lambdaijk))
# =============================================================================
	acc = 0
	for z in range(len(mod)-1):
		fr = float(mod[str(z+1)][0][0])
		acc +=fr**2 * (vijk[i][i][z] - vijk[j][j][z] -vijk[k][k][z])**2/kbtcm/2/math.sinh(fr/2/kbtcm)
	kna = c*Vijktilde**2*(math.pi/lambdaijk/kbtcm)**.5
	ka = c/2/math.pi * (1/lambdaijk * acc)**.5
	res = kna*ka/(kna+ka)

	return res
			
				
def wijk(freqs,lambdaijk,kijk):
	"""
	PRE: Given the kijk
	POST: Returns the wijk
	"""
	deltaijk = (-float(freqs[0])+float(freqs[1])+float(freqs[2]))
	return kijk*math.exp(-(deltaijk+lambdaijk)**2/4/lambdaijk/kbtcm)
	

def wjki(freqs,lambdaijk,kijk):
	"""
	PRE: Given the kijk
	POST: Returns the wjki
	"""
	deltaijk = (-float(freqs[0])+float(freqs[1])+float(freqs[2]))
	return kijk*math.exp(-(-deltaijk+lambdaijk)**2/4/lambdaijk/kbtcm)

def rlij(wi_jl, wjl_i, wij_l, wl_ij, l, i, j):
	"""
	PRE : Given the wijk and wjki
	POST: Returns the rlij
	"""
	res = (1+l==j)**3*(
			wi_jl * (1+occ_num(float(mod[str(i+1)][0][0]))) * occ_num(float(mod[str(j+1)][0][0])) - (
					wjl_i * (occ_num(float(mod[str(i+1)][0][0])))*(1+occ_num(float(mod[str(j+1)][0][0])))
					)
	) + (1+i==j)*(
			wij_l*(1+occ_num(float(mod[str(i+1)][0][0])))*(1+occ_num(float(mod[str(j+1)][0][0]))) - (
					wl_ij*occ_num(float(mod[str(i+1)][0][0]))*occ_num(float(mod[str(j+1)][0][0]))
					)
			)
	return res

def taulint(l, rijkarray):
	"""
	PRE  : Given rlij, l is a mode number starting a 0
	POST : Returns taulint
	"""
	acc=0
	for i in range(len(mod)-1):
		for j in range(len(mod)-1):
			acc+= rijkarray[l][i][j]

	return 1./acc

def taul(tlint):
	"""
	PRE: Given taulint
	POST: Returns taul
	"""
# =============================================================================
# 	return tlint
# =============================================================================
	return (1./50e-12 + 1/tlint)**-1


def tauijk(mode_nbrs, mod, lambdaijk, taull, vijk):
	"""
	PRE  : Given taul array and lambda ijk and cubic coefficient array the mode numbers will start at 0
	POST : Returns tauijk
	"""
	acc=0 

	i,j,k = mode_nbrs[0], mode_nbrs[1], mode_nbrs[2]
	if lambdaijk[i][j][k] < 0 :
		print "lambda is negative"
	for l in range(len(taull)):
		if l not in mode_nbrs:
			acc+= (vijk[l][i][i]- vijk[l][j][j] - vijk[l][k][k])**2/8/kbtcm/taull[l] * occ_num(float(mod[str(l+1)][0][0]))
			if taull[l]<0:
				pass
# =============================================================================
# 				print "a-cc is negative"
# =============================================================================
	res= (1/lambdaijk[i][j][k]*acc)**-1
	return res



def compute_decay_rates_burin(cubic_dic, mode_dic, mode_nbr, xyl_mode_list, CB_mode_list):
	"""
	PRE  : Takes a cubic dictionary for a complex, the mode dictionary for frequency information, the complex mode number to examine, the list of xylene assigned mode numbers, and the CB assigned mode numbers
	POST :  Returns the decay rate using the self consistent approach used by 
	Burin, Alexander L., et al. "Semiclassical model for vibrational dynamics in polyatomic molecules: investigation of internal vibrational relaxation." The Journal of Physical Chemistry C 114.48 (2010): 20510-20517.
 ........................................................
 :        CUBIC FORCE CONSTANTS IN NORMAL MODES         :
 :                                                      :
 : FI =  Reduced values [cm-1]  (default input)         :
 : k  =  Cubic Force Const.[AttoJ*amu(-3/2)*Ang(-3)]    :
 : K  =  Cubic Force Const.[Hartree*amu(-3/2)*Bohr(-3)] :
 :......................................................:
	"""
	h = 6.62607015e-34 # Js
	kb = 1.380649e-23 # J/K
	T = 20 # K
	c = 29979245800# cm/s

def make_vijk_array(cub):
	"""
	PRE :  Takes in the cub 
	POST: Returns an array
	"""
	dim = len(mod)-1
	cubarray = np.zeros([dim, dim, dim])
	for keys in cub:
		parts = [int(x) for x in keys.split('-')]
		cubarray[parts[0]-1][parts[1]-1][parts[2]-1] = float(cub[keys][0])
	return cubarray

def rescale_vijk_array(vijk, mod):
	"""
	PRE : Takes in the vijk array and rescale it according to the paper of tesar
	"""
	dim = len(mod)-1
	vijk_new = np.zeros([dim, dim, dim])
	for i in range(dim):
		for j in range(dim):
			for k in range(dim):
				acc = 0 		
				for z in range(dim):
					fr = float(mod[str(z+1)][0][0])
					acc += (vijk[i][i][z] - vijk[j][j][z] - vijk[k][k][z])**2 / 2 /fr**2 * math.tanh(fr/4/kbtcm) * (1-((fr)/(4*kbtcm*math.sinh(fr/kbtcm)))**2)
					
					
				tmp = vijk[i][j][k]*math.exp(-acc)
				vijk_new[i][j][k] = tmp
	return vijk_new
def occ_num(omega):
	"""
	PRE: Takes in an angular frequency in cm-1 and a temperature in K
	POST: Returns the occupation number
	"""
	return (math.exp(h*omega*c/T/kb)-1)**(-1)

def llike(oga, ogb, ogc, oi, oj, ok, gi, gj, gk):
	return (oga+ogb+ogc)/((oi+oj+ok)**2+(gi+gj+gk)**2/4)

def get_omega_k(omegak, k, cubdic, omega0, omegas, gammas):
	"""
	PRE: Starting from the initial frequencies and the corrected frequencies and the decay rate
	POST: Will return the correction to omega_k
	"""
	acc = 0
	for i in range(len(omegas)):
		for j in range(len(omegas)):
			if i == k or j == k: continue
			omegai, omegaj, omegak, gammai, gammaj, gammak = omegas[i],omegas[j],omegak,gammas[i],gammas[j],gammas[k]
			ni = occ_num(omegai)
			nj = occ_num(omegaj)
			key = "{}-{}-{}".format(i+1, j+1, k+1)
			
			try:
				if float(cubdic[key][0])**2<1e15: # If the anharmonicities are too big, the computation just crashes
					acc += float(cubdic[key][0])**2 * (
							(ni+nj+1)*(llike(omegai, omegaj, -omegak, omegai, omegaj, -omegak, gammai, gammaj, gammak)) + 
							(ni+nj+1)*(llike(omegai, omegaj, omegak, omegai, omegaj, omegak, gammai, gammaj, gammak)) +
							(nj-ni)*(llike(omegai, -omegaj, -omegak, omegai, -omegaj, -omegak, gammai, gammaj, gammak)) +
							(ni-nj)*(llike(omegai, -omegaj, -omegak, omegai, -omegaj, -omegak, gammai, gammaj, gammak))
							)
			except:
				continue
			
	return omega0[k] - acc/16  - omegak


def get_decay_k(gammak, k, cubdic, omega0, omegas, gammas):
	"""
	PRE: Starting from the initial frequencies and the corrected frequencies and the decay rate
	POST: Will return the correction to gamma_k
	"""
	acc = 0
	print "step {}".format(k)
	for i in range(len(omegas)):
		for j in range(len(omegas)):
			if i == k or j == k: continue
			omegai, omegaj, omegak, gammai, gammaj, gammak = omegas[i],omegas[j],omegas[k],gammas[i],gammas[j],gammak
			ni = occ_num(omegai)
			nj = occ_num(omegaj)
			key = "{}-{}-{}".format(i+1, j+1, k+1)
			
			try:
				if float(cubdic[key][0])**2<1e15: # If the anharmonicities are too big, the computation just crashes
					temp = float(cubdic[key][0])**2 * (
							(ni+nj+1)*(llike(gammai, gammaj, gammak, omegai, omegaj, -omegak, gammai, gammaj, gammak)) + 
							(ni+nj+1)*(llike(gammai, gammaj, gammak, omegai, omegaj, omegak, gammai, gammaj, gammak)) +
							(nj-ni)*(llike(gammai, gammaj, gammak, omegai, -omegaj, -omegak, gammai, gammaj, gammak)) +
							(ni-nj)*(llike(gammai, gammaj, gammak, omegai, -omegaj, -omegak, gammai, gammaj, gammak))
							)
					acc+=temp
			except:
				continue
	return acc/16 - gammak

def breakdown_decay_by_mode(k, ARG_LIST):
	"""
	PRE: Takes in pre optimized gammas as ARG_LIST = cubdic, omegas, gammas, xyl_mode_list, CB_mode_list
	POST: Returns a list with index i correponds to the ith mode contribution to decay (probs times 2)
	"""
	print k
	biggest_list = []
	cubdic, omegas, gammas, xyl_mode_list, CB_mode_list = ARG_LIST
	kthdecay_breakdown = [0]*len(omegas)
	WXX, WCC, WCX, WXM, WCM, WMM = 0,0,0,0,0,0 # cumulative decay rate for xylene xylene denominated frequencies (WXX), CB and CB (WCC), CB and Xylene (WCX), Mixed and Mixed (WMM), Xylene and Mixed (WXM), CB and Mixed (WCM)
	cWXX, cWCC, cWCX, cWXM, cWCM, cWMM = 0,0,0,0,0,0 # associated counters
	acc = 0
	for i in range(len(omegas)):
		for j in range(len(omegas)):
			if i == k or j == k: continue
			nums = [i,j]
			omegai, omegaj, omegak, gammai, gammaj, gammak = omegas[i],omegas[j],omegas[k],gammas[i],gammas[j],gammas[k]
			ni = occ_num(omegai)
			nj = occ_num(omegaj)
			key = "{}-{}-{}".format(i+1, j+1, k+1)
			
			try:
				if float(cubdic[key][0])**2<1e15: # If the anharmonicities are too big, the computation just crashes
					temp = 1./16*float(cubdic[key][0])**2 * (
							(ni+nj+1)*(llike(gammai, gammaj, gammak, omegai, omegaj, -omegak, gammai, gammaj, gammak)) + 
							(ni+nj+1)*(llike(gammai, gammaj, gammak, omegai, omegaj, omegak, gammai, gammaj, gammak)) +
							(nj-ni)*(llike(gammai, gammaj, gammak, omegai, -omegaj, -omegak, gammai, gammaj, gammak)) +
							(ni-nj)*(llike(gammai, gammaj, gammak, omegai, -omegaj, -omegak, gammai, gammaj, gammak))
							)
					kthdecay_breakdown[i]+=temp
					kthdecay_breakdown[j]+=temp
				
					if all([x in xyl_mode_list for x in nums]):
						WXX+=temp
						cWXX+=1
					elif all([x in CB_mode_list for x in nums]):
						WCC+=temp
						cWCC+=1
					elif any([x in xyl_mode_list for x in nums]) and any([x in CB_mode_list for x in nums]): 
						WCX+=temp
						cWCX+=1
					elif any([x in xyl_mode_list for x in nums]):
						WXM+=temp
						cWXM+=1
					elif any([x in CB_mode_list for x in nums]):
						WCM+=temp
						cWCM+=1
					else:
						WMM+=temp
						cWMM+=1
					acc+=temp
#						print "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(k, i, j, omegak, omegai, omegaj, temp)
					biggest_list.append((k, i, j, omegak, omegai, omegaj, temp))
			except:
# =============================================================================
# 				print(traceback.format_exc())
# =============================================================================
				continue
	for w in sorted(biggest_list, key = lambda x: x[-1], reverse=True)[:20]:
		print w
	print "Mode freq: ", omegas[k], "overall: ", sum([WXX, WCC, WCX, WXM, WCM, WMM]), acc
	print "counts: ", cWXX, cWCC, cWCX, cWXM, cWCM, cWMM
#	print "avg rate contribution: ",W1/cw1, W2/cw2, W3cw3
	print "Overall rates: ", zip(['WXX', 'WCC', 'WCX', 'WXM', 'WCM', 'WMM'], [WXX, WCC, WCX, WXM, WCM, WMM])
	return zip(['WXX', 'WCC', 'WCX', 'WXM', 'WCM', 'WMM'], [WXX, WCC, WCX, WXM, WCM, WMM])

if __name__ == "__main__":
	import numpy as np
	import glob
	import math
	import itertools
	import cPickle
	import scipy
	import scipy.optimize
	import traceback
	from multiprocessing import Pool
	import matplotlib.pyplot as plt
	h = 6.62607015e-34 # Js
	kb = 1.380649e-23 # J/K
	T = 300 # K
	c = 29979245800# cm/s
	amu = 1.66054e-27 # kg
	kbtcm = kb*T/h/c

	flist = ["/Users/hugueslambert/Desktop/xylene/cubic_coupling/{}".format(x) for x in [
			"mxylene.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out",
			"mxylene-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out",
			"mxylene-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out",
			"mxylene-protonated-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out",
			"mxylene-protonated-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out"]]
#	flist = sorted(glob.glob("/Users/hugueslambert/Desktop/xylene/cubic_coupling/mxylene-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out"))
	print flist
#	flist = ["/Users/hugueslambert/Desktop/xylene/cubic_coupling/mxylene.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out"]
# =============================================================================
# 	flist = ["/home/macenrola/Documents/XYLENE/inputs/cubic_coupling/OUTS/cubic_coupling/mxylene.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out"]
# =============================================================================
## GAMMA OPTI WITH SCF
#	def get_resg(k0, CUB_O0_O_G):
#		cub, freqs0, temp_omega, temp_gamma = CUB_O0_O_G
#		return scipy.optimize.brentq(get_decay_k, -10, 20000, args=(k0,cub,freqs0,temp_omega,temp_gamma))
#	
#	def star_get_resg(k0_other_args):
#		return get_resg(*k0_other_args)
#	p = Pool(None)
#
#	for f in flist:
#		print f
#		mod, cub = get_data_from_out(f)
#		freqs0 = sorted([float(mod[x][0][0]) for x in sorted(mod.keys())[:-1]])
#		try:
#			with open(f+"_GAMMAS_seqo", "rb") as r:
#				gammas0 = [float(x.strip()) for x in r.readlines()]
#		except:
#			gammas0 = [150.]*len(freqs0)
#
#		print gammas0
#		cub = make_permutations(cub)
#		dim = len(mod)-1
#		
#		# ALL SCF 
#		temp_omega = freqs0
#		temp_gamma = gammas0
#		convo, convg = 0e10, 1e10
#		for z in range(100):
#			if (convo+convg) < 0.001: break
#			print "iteration {}, total residues {}".format(z, convo+convg)
#			cgamma =[]
#			print "POOL IS STARTING"
#			gamma_num = range(len(freqs0))
#			second_arg = [cub, freqs0, freqs0, temp_gamma]
#			cgamma = p.map(star_get_resg,  itertools.izip(gamma_num, itertools.repeat(second_arg)))
#			print "POOL IS OVER"
#			convg = sum([(x-y)**2 for x,y in zip(cgamma, temp_gamma)])
#			print "it:{}, difference in gamma {}".format(z, convg), ["{0:3.3f}".format(x-y) for x,y in zip(cgamma, temp_gamma)]
#			temp_gamma = cgamma
#
#		print freqs0	
#		print "final omega", temp_omega
#		print "final gamma", [x/1e12*c for x in temp_gamma]
#		with open(f+"_GAMMAS", "wb") as w:
#			w.writelines(["{}\n".format(x) for x in temp_gamma])
	
	
	
# =============================================================================
# ## ASSIGNAMTIONS OF THE GAMMAS	
	def star_breakdown_decay_by_mode(ALL_ARGS):
		return breakdown_decay_by_mode(*ALL_ARGS)
	for f in flist[:3]:
		mod, cub = get_data_from_out(f)
	 	cub = make_permutations(cub)
	 	xyl_mod = get_mode_localisation_on_xylene(mod, 18)
	 	for w in sorted(xyl_mod):
			 print w
	 	with open(f + "_GAMMAS", "rb") as r:
	 		gammas0 = [float(x.strip()) for x in r.readlines()]
			print len(gammas0)
	 	freqs0 = sorted([float(mod[x][0][0]) for x in sorted(mod.keys())[:-1]])
		print len(freqs0)
		plottable_xyl_loc = [] # as tuples (xylene localisation number, decay to xylene, total decay, mode number, mode frequency, mode number as defined by xylene localistaion )
#		p = Pool(4)
#		second_arg = [cub, freqs0, gammas0, [y-1 for x,y,_ in xyl_mod if x>0.90], [y-1 for x,y,_ in xyl_mod if x<0.1]]
#		breakdowns = p.map_async(star_breakdown_decay_by_mode,  itertools.izip(range(len(freqs0)), itertools.repeat(second_arg)))
#		p.close()
	 	for k in range(len(freqs0)):#[302]:#range(359,375):#[317,318]:##[365]:
#		k = 361
			breakdown = breakdown_decay_by_mode(k, [cub, freqs0, gammas0, [y-1 for x,y,_ in xyl_mod if x>0.90], [y-1 for x,y,_ in xyl_mod if x<0.1] ])
			plottable_xyl_loc.append((xyl_mod[k][0], breakdown[0][1]+breakdown[2][1]+breakdown[3][1], sum([x[1] for x in breakdown]), k, float(xyl_mod[k][2]), sorted(xyl_mod, key= lambda x: x[0]).index(xyl_mod[k])))
	 	print plottable_xyl_loc
	# # =============================================================================
	# # 	l = sorted(plottable_xyl_loc)[::-1]
	# # =============================================================================
	 	l = sorted(plottable_xyl_loc)[::-1]
#	 	plt.plot([x[0] for x in l])
#	 	plt.plot([x[1]/max([x[2] for x in l]) for x in l])
#	 	plt.plot([x[2]/max([x[2] for x in l]) for x in l])
#	 	
#	 	plt.show()
		cPickle.dump(l, open(f+'_graph_dump', 'wb'))
# =============================================================================
#	for ls in glob.glob("/Users/hugueslambert/Desktop/xylene/cubic_coupling/*graph_dump")[-1:]:
#		l= cPickle.load(open(ls, 'rb'))
#	 	plt.plot([x[0] for x in l])
#	 	plt.plot([x[1]/max([x[2] for x in l]) for x in l])
#	 	plt.plot([x[2]/max([x[2] for x in l]) for x in l])
# 	
#	 	plt.show()
# =============================================================================
# 	mod, cub = get_data_from_out('/home/macenrola/AcPhCN.com_OUT.out')
# =============================================================================
# =============================================================================
# 	flist = glob.glob("/home/macenrola/Documents/XYLENE/inputs/cubic_coupling/OUTS/cubic_coupling/*.out")
# =============================================================================
# =============================================================================
# 	flist = ["/home/macenrola/Documents/XYLENE/inputs/cubic_coupling/OUTS/cubic_coupling/mxylene-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out"]
# =============================================================================

# =============================================================================
# 	flist = glob.glob("/home/macenrola/Documents/XYLENE/inputs/cubic_coupling/OUTS/cubic_coupling/xylenes/*.out")
# =============================================================================
# =============================================================================
# 	mod, cub = get_data_from_out(flist[0])
# 	print mod.keys()
# 	for k, i in enumerate(["{}\t{}\t{}".format(x, y, z) for x, y, z in get_mode_localisation_on_xylene(mod, 18)]):
# 		print k,i
# =============================================================================

		
	

