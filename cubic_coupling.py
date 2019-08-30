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
	for k in sorted(mode)[:-1]:
		LOC_ON_XYL_COEFF=sum([float(x)**2 for x in mod[k][1][:mode_xyl_atoms]])
		xyl_localised_modes.append((LOC_ON_XYL_COEFF, k))

	return sorted(xyl_localised_modes)[::-1]
	

def occ_num(omega, h, kb, T):
	"""
	PRE: Takes in an angular frequency in cm-1 and a temperature in K
	POST: Returns the occupation number
	"""
	return (math.exp(h*omega/T/kb)-1)**(-1)

def C_FUJISAKI(omega, h, kb, T):
	"""
	PRE   : Takes in frequency, h, kb and T
	POST  : Will return the parameter used in the FUJISAKI paper page
	Fujisaki, Hiroshi, Lintao Bu, and John E. Straub. "Vibrational energy relaxation (VER) of a CD stretching mode in cytochrome c." arXiv preprint q-bio/0403019 (2004).
	"""
	return 1/h/omega*(1-math.exp(-h*omega/kb/T))/(1+math.exp(-h*omega/kb/T))

def delta(omega, c, l=3):
	"""
	PRE: Takes in frequency in Hz with the lifetime given in cm-1
	POST: Returns the value of the thawed delta according to Fujisaki, Hiroshi, Lintao Bu, and John E. Straub. "Vibrational energy relaxation (VER) of a CD stretching mode in cytochrome c." arXiv preprint q-bio/0403019 (2004).
	"""
# =============================================================================
# 	if abs(omega/29979245800.)<1:
# 		print "resonance: diff is {} delta is {}".format(omega/29979245800., (1/math.pi * l/(l**2+(omega)**2)))
# 		return 0
# =============================================================================
	return 1/math.pi * l*c/((l*c)**2+(omega)**2)
	

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
	h = 6.62607015e-34 # Js
	kb = 1.380649e-23 # J/K
	T = 20 # K
	amu = 1.66054e-27 # kg
	c = 29979245800# cm/s
	WXX, WCC, WCX, WXM, WCM, WMM = 0,0,0,0,0,0 # cumulative decay rate for xylene xylene denominated frequencies (WXX), CB and CB (WCC), CB and Xylene (WCX), Mixed and Mixed (WMM), Xylene and Mixed (WXM), CB and Mixed (WCM)
	cWXX, cWCC, cWCX, cWXM, cWCM, cWMM = 0,0,0,0,0,0 # associated counters
	mode_nbr = str(mode_nbr)
	for k in cubic_dic:
		cwdecay = 0
		nums = k.split('-')
		freqs = [float(mode_dic[x][0][0])*c for x in nums] # Gets for mode num x the metadata (0) containing the frequency (0) in cm-1
		masses = [float(mode_dic[x][0][1]) for x in nums] #masses in amu
		if mode_nbr == nums[0]:
			fcmode = freqs[nums.index(mode_nbr)]
			nums.remove(mode_nbr)
			freqs.remove(fcmode)
# =============================================================================
# 			if abs(fcmode-freqs[0]-freqs[1])/29979245800.<1.:
# 				print "resonance: diff is {} modes are {} {} {}, coupling coeff is :{}".format((fcmode-freqs[0]-freqs[1])/29979245800., nums[0], nums[1], mode_nbr, float(cubic_dic[k][1])*1e-18*6e23/4.18e3)
# =============================================================================
			#cwdecay = h/16./fcmode*(float(cubic_dic[k][1])*1e-18*amu**(-3/2.0)*1e-10**(-3))**2*masses[0]*masses[1]*masses[2]/freqs[0]/freqs[1]*(1+occ_num(freqs[0], h, kb,T)+ occ_num(freqs[1],h, kb,T))*delta(fcmode-freqs[0]-freqs[1])
			cwdecay = h/16.*(float(cubic_dic[k][1])*1e-18*amu**(-3./2.0)*1e-10**(-3))**2 /freqs[0]/freqs[1]/fcmode*(1+occ_num(freqs[0], h, kb,T) + occ_num(freqs[1],h, kb,T))*delta(fcmode-freqs[0]-freqs[1],c)
			#cwdecay = h/16.*(float(cubic_dic[k][0])*c*2*math.pi)**2 /freqs[0]/freqs[1]/fcmode*(1+occ_num(freqs[0], h, kb,T) + occ_num(freqs[1],h, kb,T))*delta(fcmode-freqs[0]-freqs[1],c)
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
	


def get_lambda(mode_nbrs, mod, cub):
	"""
	PRE: Takes in a mode number, a cubic dictionary and mode dictionary and K and T
	POST: Returns the lambda_ijk value accordking to buttrin. It is assumed that the cm-1 format can yield an energy through *c*h which then can give lamdaijk as energy once divided by kb*T
	"""	
	cubic_fluctuation = 0
	temp_cub_dic = {}
	mode_nbrs = [str(x) for x in mode_nbrs]
	for i in sorted(cub):
		nums = i.split('-')
		for els in mode_nbrs:
			if nums.count(els) == 2:
				nums.remove(els)
				nums.remove(els)
				if int(nums[0]) not in temp_cub_dic:
					temp_cub_dic[int(nums[0])] = [(i, cub[i])]
				else:
					temp_cub_dic[int(nums[0])].append((i, cub[i]))
	for k in sorted(temp_cub_dic):
		fi = -1
		if len(temp_cub_dic[k])==3:
			for i, els in enumerate(temp_cub_dic[k]):
				parts=els[0].split('-')
				if parts.count(mode_nbrs[0]) == 2:
					parts.remove(mode_nbrs[0])
					parts.remove(mode_nbrs[0])
					fi=i
					break
			index = [0,1,2]
			index.remove(fi)
			cubic_fluctuation += ((float(temp_cub_dic[k][fi][1][0])-float(temp_cub_dic[k][index[0]][1][0])-float(temp_cub_dic[k][index[1]][1][0]))*c*h)**2*occ_num(float(mod[parts[0]][0][0])*c,h, kb,T)
	return .25*cubic_fluctuation/T/kb
# =============================================================================
# 		if i not in mode_nbrs:
# 			try:
# 				temp_flux = []
# 				for els in mode_nbrs:
# 					temp_flux.append(float(cub['{}-{}-{}'.format(i, els,els)][0]))
# 				for els in mode_nbrs:
# 					print '{}-{}-{}'.format(i, els,els)
# 				cubic_fluctuation+= ((temp_flux[0]-temp_flux[1]-temp_flux[2])*c*h)**2*occ_num(float(mod[i][0][1])*c,h, kb,T)
# 			except:
# 				continue
# 		else: 
# 			continue
# 	print .25*cubic_fluctuation/kb/T/h/c
# =============================================================================
		
def kijk(i,j,k, Vijk, lambdaijk, tauijk):
	"""
	PRE: Given the Tauijk 
	POST: Returns the kijk
	"""
	Vijktilde = Vijk/(1+(j==k))
	hbar = h/math.pi/2
	return (2*math.pi) /hbar * (Vijktilde)**2 *1/8/(4*math.pi * kb * T* lambdaijk)**.5 *1/(
			1+ (4*math.pi*Vijktilde**2*tauijk)/(8*hbar*lambdaijk))

def wijk(mod, i,j,k,lambdaijk,kijk):
	"""
	PRE: Given the kijk
	POST: Returns the wijk
	"""
	deltaijk = h/math.pi/2*(-mod[i][0][0]+mod[j][0][0]+mod[k][0][0])*c
	return kijk*math.exp(-(deltaijk+lambdaijk)**2/4/lambdaijk/kb/T)
	

def wjki(mod,i,j,k,lambdaijk,kijk):
	"""
	PRE: Given the kijk
	POST: Returns the wjki
	"""
	deltaijk = h/math.pi/2*(-mod[i][0][0]+mod[j][0][0]+mod[k][0][0])*c
	return kijk*math.exp(-(-deltaijk+lambdaijk)**2/4/lambdaijk/kb/T)

def rlij(mod, l, i, j, Vijk):
	"""
	PRE : Given the wijk and wjki
	POST: Returns the rlij
	"""
	lambdaijl = get_lambda([i,j,l], mod, cub)
	tauijl = tauijk([i,j,l], lambdaijl)
	kijl = kijk(i,j,l, Vijk, lambdaijl, tauijl)
	lambdajli = get_lambda([j,l,i], mod, cub)
	taujli = tauijk([j,l,i], lambdajli)
	kjli = kijk(j,l,i, lambdajli, taujli)
	lambdaijl = get_lambda([i,j,l], mod, cub)
	tauijl = tauijk([i,j,l], lambdaijl)
	kijl = kijk(i,j,l, lambdaijl, tauijl)
	lambdalij = get_lambda([l,i,j], mod, cub)
	taulij = tauijk([l,i,j], lambdalij)
	klij = kijk(l,i,j, lambdalij, taulij)
	
	return (1+l==j)**3*(
			wijk(mod, i, j, l, lambdaijl, kijl) * (1+occ_num(float(mod[i][0][0])*c,h, kb,T)) * occ_num(float(mod[j][0][0])*c,h, kb,T) - (
					wjki(mod, j, l, i, lambdajli, kjli) * (occ_num(float(mod[i][0][0])*c,h, kb,T))*(1+occ_num(float(mod[j][0][0])*c,h, kb,T))
					)
	) + (1+i==j)*(
			wjki(mod, i, j, l, lambdaijl, kijl)*(1+occ_num(float(mod[i][0][0])*c,h, kb,T))*(1+occ_num(float(mod[j][0][0])*c,h, kb,T)) - (
					wijk(mod, l,i,j, lambdalij, klij)*occ_num(float(mod[i][0][0])*c,h, kb,T)*occ_num(float(mod[j][0][0])*c,h, kb,T)
					)
			)
			

def taulint(l):
	"""
	PRE  : Given rlij
	POST : Returns taulint
	"""
	acc=0
	for i in range(len(mod)-1):
		for j in range(len(mod)-1):
# =============================================================================
# 			try:
# =============================================================================
			vijk = float(cub["{}-{}-{}".format(*sorted([l,i+1,j+1])[::-1])][0])*c*h
			temp = rlij(mod, l, i+1, j+1, vijk )
			print temp 
			acc+=temp
# =============================================================================
# 			except:
# 				pass
# =============================================================================
	return acc

def taul(l):
	"""
	PRE: Given taulint
	POST: Returns taul
	"""
	return (1/10e-12 + 1/taulint(l))**-1
# =============================================================================
# 	return 50e-12
# =============================================================================

def tauijk(mode_nbrs, lambdaijk):
	"""
	PRE  : Given taul
	POST : Returns tauijk
	"""
	acc=0 
	temp_cub_dic = {}
	for i in range(len(mod)-1):
		print mode_nbrs
		if i+1 in mode_nbrs:
			continue
		else:
			temp_cub_dic[i+1] = []
			for els in mode_nbrs:
				try:
					key = "{}-{}-{}".format(*sorted([i+1, els, els])[::-1])
					temp_cub_dic[i+1].append((key, cub[key]))
				except:
					print "not here"

# =============================================================================
# 	for i in sorted(cub):
# 		nums = i.split('-')
# 		for els in mode_nbrs:
# 			if nums.count(els) == 2:
# 				nums.remove(els)
# 				nums.remove(els)
# 				if int(nums[0]) not in temp_cub_dic:
# 					temp_cub_dic[int(nums[0])] = [(i, cub[i])]
# 				else:
# 					temp_cub_dic[int(nums[0])].append((i, cub[i]))
# =============================================================================
	print temp_cub_dic
	mode_nbrs = [str(x) for x in mode_nbrs]
	for k in sorted(temp_cub_dic):
		fi = -1
		print k, temp_cub_dic[k]
		if len(temp_cub_dic[k])==3:
			for i, els in enumerate(temp_cub_dic[k]):
				parts=els[0].split('-')
				if parts.count(mode_nbrs[0]) == 2:
					parts.remove(mode_nbrs[0])
					parts.remove(mode_nbrs[0])
					fi=i
					break
			index = [0,1,2]
			index.remove(fi)
			acc += ((float(temp_cub_dic[k][fi][1][0])-float(temp_cub_dic[k][index[0]][1][0])-float(temp_cub_dic[k][index[1]][1][0]))*c*h)**2*occ_num(float(mod[parts[0]][0][0])*c,h, kb,T)/taul(k)
	return (1/lambdaijk*acc/8/T/kb)**-1
	
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
	
if __name__ == "__main__":
	import numpy as np
	import glob
	import math
	h = 6.62607015e-34 # Js
	kb = 1.380649e-23 # J/K
	T = 20 # K
	c = 29979245800# cm/s
	amu = 1.66054e-27 # kg

# =============================================================================
# 	#for f in glob.glob('/Users/hugueslambert/Desktop/xylene/cubic_coupling/*out'):
# 	#	print len(get_data_from_out(f)[0].keys())
# 	path = '/home/macenrola/Documents/XYLENE/inputs/cubic_coupling/OUTS/freqs_from_tight_geom/'
# 	#mod_xyl, _ = get_data_from_out("/Users/hugueslambert/Desktop/xylene/cubic_coupling/mxylene.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out")
# 	couples = [
# 			'mxylene-h-single-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
# 			'mxylene-h-single-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
# 			'mxylene-protonated-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
# 			'mxylene-protonated-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
# 			'mxylene-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
# 			'mxylene-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
# 			]
# 	for f in couples[0:1]:
# 		print f
# 		mod, cub= get_data_from_out(path+f)
# 		mod_xyl, _ = get_data_from_out("{}just_xylenes/just_xyl-{}.xyz.com_OUT.out".format(path,f))
# 		i=0
# 
# 		xyl_localised_modes=get_mode_localisation_on_xylene(mod, len(mod_xyl['1'][1]))
# 		for i,k in enumerate(xyl_localised_modes):
# 			print k[0], k[1], i
# 		lim_xyl=121 # 121 limit for 10% on xylene, 90% on xylene is 28 and 99 % is 16
# 		start_cb=210 # 210 limit for 99% on CB, 122 limit of 90% CB
# 		WXX, WCC, WCX, WXM, WCM, WMM = 0,0,0,0,0,0
# 		for els in  [x[1] for x in xyl_localised_modes[-start_cb:]]:
# 			print els
# 			tWXX, tWCC, tWCX, tWXM, tWCM, tWMM = compute_decay_rates(cub, mod, els, [x[1] for x in xyl_localised_modes[:lim_xyl]],  [x[1] for x in xyl_localised_modes[start_cb:]])
# 			WXX+=tWXX
# 			WCC+=tWCC
# 			WCX+=tWCX
# 			WXM+=tWXM
# 			WCM+=tWCM
# 			WMM+=tWMM
# 			
# 			
# 		print "cumulated values are: ",WXX, WCC, WCX, WXM, WCM, WMM
# =============================================================================
# =============================================================================
# 	mod, cub = get_data_from_out('/home/macenrola/AcPhCN.com_OUT.out')
# 	WXX, WCC, WCX, WXM, WCM, WMM = 0,0,0,0,0,0
# 	for i in sorted([int(x) for x in mod.keys()[:-1]]):
# 		print i, "Freq is: {}".format(mod[str(i)][0][0])
# 		tWXX, tWCC, tWCX, tWXM, tWCM, tWMM = compute_decay_rates_fujisaki(cub, mod, i, sorted(mod.keys())[:-1],  sorted(mod.keys())[:-1])
# 		WXX+=tWXX
# 		WCC+=tWCC
# 		WCX+=tWCX
# 		WXM+=tWXM
# 		WCM+=tWCM
# 		WMM+=tWMM
# 	print "cumulated values are: ", WXX, WCC, WCX, WXM, WCM, WMM
# =============================================================================
	mod, cub = get_data_from_out('/home/macenrola/AcPhCN.com_OUT.out')

	l= get_lambda([49-x for x in [8,45,47]], mod, cub)
	print l/h/c
# =============================================================================
# 	print tauijk([8,21,20], l)
# =============================================================================
	