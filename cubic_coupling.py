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
		meta = [x.strip().split() for x in mode_block[i+1:i+5]] # gets the frequency, reduced masses, force constans and IR intensities
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

if __name__ == "__main__":
	import numpy as np
	import glob
	#for f in glob.glob('/Users/hugueslambert/Desktop/xylene/cubic_coupling/*out'):
	#	print len(get_data_from_out(f)[0].keys())
	path = '/home/macenrola/Documents/XYLENE/inputs/cubic_coupling/OUTS/cubic_coupling/'
	#mod_xyl, _ = get_data_from_out("/Users/hugueslambert/Desktop/xylene/cubic_coupling/mxylene.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out")
	couples = [
			'mxylene-h-single-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-h-single-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-protonated-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-protonated-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-CB6.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			'mxylene-CB7.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out.xyz.com_OUT.out',
			]
	for f in couples:
		print f
		mod,_ = get_data_from_out(path+f)
		mod_xyl, _ = get_data_from_out("{}xylenes/just_xyl-{}.xyz.com_OUT.out".format(path,f))
		procrustean_fit_for_xylene_part_of_eigenvectors(mod_xyl, mod)
# =============================================================================
# 		for k in sorted(mod)[:-1]: # avoid the order 
# 			for j in sorted(mod_xyl)[:-1]: # leave out the last one to avoid the ORDER that contains the order of the atoms defining the modes, it's assumed they line up
# 				corr= sum([float(x)*float(y) for x, y in zip(mod_xyl[j][1], mod[k][1][:len(mod_xyl[j][1])])])
# 				if corr>0.5:
# 					print k, j , corr
# =============================================================================
