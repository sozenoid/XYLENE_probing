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

def get_mode_localisation_on_xylene(mode, mode_xyl):
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
		LOC_ON_XYL_COEFF=sum([float(x)**2 for x in mod[k][1][:len(mode_xyl['1'][1])]])
		xyl_localised_modes.append((LOC_ON_XYL_COEFF, k))

	return sorted(xyl_localised_modes)[::-1]
	

def occ_num(omega, h, k, T):
	"""
	PRE: Takes in an angular frequency in cm-1 and a temperature in K
	POST: Returns the occupation number
	"""
	return (math.exp(h*omega/T/k)-1)**(-1)

def delta(omega, l=3*29979245800):
	"""
	PRE: Takes in frequency in cm-1 with the lifetime given in cm-1
	POST: Returns the value of the thawed delta according to Fujisaki, Hiroshi, Lintao Bu, and John E. Straub. "Vibrational energy relaxation (VER) of a CD stretching mode in cytochrome c." arXiv preprint q-bio/0403019 (2004).
	"""
	return 1/math.pi * l/(l**2+omega**2)
	

def compute_decay_rates(cubic_dic, mode_dic, mode_nbr, xyl_mode_list, CB_mode_list):
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
	T = 300 #K
	amu = 1.66054e-27 # kg
	WXX, WCC, WCX, WXM, WCM, WMM = 0,0,0,0,0,0 # cumulative decay rate for xylene xylene denominated frequencies (WXX), CB and CB (WCC), CB and Xylene (WCX), Mixed and Mixed (WMM), Xylene and Mixed (WXM), CB and Mixed (WCM)
	cWXX, cWCC, cWCX, cWXM, cWCM, cWMM = 0,0,0,0,0,0 # associated counters
	mode_nbr = str(mode_nbr)
	for k in cubic_dic:
		cwdecay = 0
		nums = k.split('-')
		freqs = [float(mode_dic[x][0][0])*29979245800 for x in nums] # Gets for mode num x the metadata (0) containing the frequency (0) in cm-1
		masses = [float(mode_dic[x][0][1]) for x in nums]
		if mode_nbr in nums:
			fcmode = freqs[nums.index(mode_nbr)]
			nums.remove(mode_nbr)
			freqs.remove(fcmode)
			
			cwdecay = (float(cubic_dic[k][1])*1e-18*1e-10**(-3))**2*masses[0]*masses[1]*masses[2]/freqs[0]/freqs[1]*(1+occ_num(freqs[0], h,kb,T)+ occ_num(freqs[1],h,kb,T))*delta(fcmode-freqs[0]-freqs[1])
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
			pass
	print "counts: ", cWXX, cWCC, cWCX, cWXM, cWCM, cWMM
	#print "avg rate contribution: ",W1/cw1, W2/cw2, W3/cw3
	print "Overall rates: ", WXX, WCC, WCX, WXM, WCM, WMM
	return [x*h/16/fcmode for x in [WXX, WCC, WCX, WXM, WCM, WMM]]
	
if __name__ == "__main__":
	import numpy as np
	import glob
	import math
	#for f in glob.glob('/Users/hugueslambert/Desktop/xylene/cubic_coupling/*out'):
	#	print len(get_data_from_out(f)[0].keys())
	path = '/home/macenrola/Documents/XYLENE/inputs/cubic_coupling/OUTS/freqs_from_tight_geom/'
	#mod_xyl, _ = get_data_from_out("/Users/hugueslambert/Desktop/xylene/cubic_coupling/mxylene.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out")
	couples = [
			'mxylene-h-single-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-h-single-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-protonated-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-protonated-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			]
	for f in couples[0:1]:
		print f
		mod, cub= get_data_from_out(path+f)
		mod_xyl, _ = get_data_from_out("{}just_xylenes/just_xyl-{}.xyz.com_OUT.out".format(path,f))
		i=0

		xyl_localised_modes=get_mode_localisation_on_xylene(mod, mod_xyl)
		for i,k in enumerate(xyl_localised_modes):
			print k[0], k[1], i
		lim_xyl=121 # 121 limit for 10% on xylene, 90% on xylene is 28 and 99 % is 16
		start_cb=210 # 210 limit for 99% on CB, 122 limit of 90% CB
		WXX, WCC, WCX, WXM, WCM, WMM = 0,0,0,0,0,0
		for els in  [x[1] for x in xyl_localised_modes[lim_xyl:start_cb]]:
			print els
			tWXX, tWCC, tWCX, tWXM, tWCM, tWMM = compute_decay_rates(cub, mod, els, [x[1] for x in xyl_localised_modes[:lim_xyl]],  [x[1] for x in xyl_localised_modes[start_cb:]])
			WXX+=tWXX
			WCC+=tWCC
			WCX+=tWCX
			WXM+=tWXM
			WCM+=tWCM
			WMM+=tWMM
			
		print "cumulated values are: ",WXX, WCC, WCX, WXM, WCM, WMM