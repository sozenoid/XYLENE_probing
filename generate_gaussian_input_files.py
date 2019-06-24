#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:51:04 2019

@author: macenrola
"""

def format_gaussian_input_from_xyz(xyz_file):
	"""
	PRE: Takes in a valid xyz file generated by openbabel or a gaussian output file
	POST: will produce a gaussian input file
	"""
	name=xyz_file.split("/")[-1]
	route="#N PM6D3 opt=(ts, noeigentest, modredundant, calcall, maxcyc=999) freq"
	freeze=" D       2       3       9      10 F"
	checkpoint="%Chk={}.chk".format(name)
	mem="%mem=120gb"
	procs="%NProcShared=24"

# =============================================================================
# 	Fine tuning the charge and multiplicity
# =============================================================================
	if 'diradical' in name or 'N2' in name:
		charge_mult="1 3"
		print "{} is diradical or n2".format(name)
	else:
		charge_mult = "2 1"
		
		
# =============================================================================
# 	Fine tuning the route section
# =============================================================================
	if 'TS' in name:
		route="#n wB97XD/6-31G(d) opt=(ts, noeigentest, modredundant, calcfc, maxcyc=999) maxdisk=100GB freq"
	else:
		route="#n wB97XD/6-31G(d) opt=(modredundant, maxcyc=999) maxdisk=100GB freq"
		
	
	if name[-4:]==".xyz":
		with open(xyz_file, 'rb') as r:
			coords = r.readlines()[2:]
	elif name[-4:]==".out":
		coords = extract_xyz_from_gaussian_out(xyz_file)

	with open(xyz_file+".com", "wb") as w:
		w.write(procs+"\n")
		w.write(checkpoint+"\n")
		w.write(mem+"\n")
		w.write(route+"\n\n")
		w.write(name+"\n\n")
		w.write(charge_mult+"\n")
		w.writelines(coords)
		w.write("\n")
		# w.write("notatoms=1-{}\n".format(len(coords)-126))
		w.write(freeze+"\n")
		w.write("\n")


def extract_xyz_from_gaussian_out(outfile):
	"""
	PRE  : Takes in a valid Gaussian File 
	POST : Returns the last xyz coordinates
	"""
	atom_dic={'1':'H', '6':'C', '7':'N', '8':'O'}
	with open(outfile, "rb") as r:
		lines = r.readlines()
	
	last_rot = -1
	for i, line in enumerate(lines):
		if "Rotational constants" in line:
			last_rot=i
	
	last_atom=-1
	for i, line in enumerate(lines):
		if "Coordinates (Angstroms)" in line:
			last_atom=i
			
	if last_rot==-1 or last_atom==-1:
		print "FILE INCOMPLETE"
		return None
	
	return ["{0}    {1:+2.8f}    {2:+2.8f}    {3:+2.8f}\n".format(atom_dic[x.strip().split()[1]], *[float(y) for y in x.strip().split()[-3:]]) for x in lines[last_atom+3: last_rot-1]]
	
if __name__ == "__main__":
	import sys
	
	if len(sys.argv)==1:
		print """
		Please provide input files
		You need to run the following lines before having the command to work
			echo "alias gengauss='python `pwd`/generate_gaussian_input_files.py'" >> ~/.bashrc
			source ~/.bashrc
		Then you can use the script in any directory using:
			gengauss *xyz
		and edit the script by going to the generate_gaussian_input_files.py file and doing
			nano generate_gaussian_input_files.py
		
		"""
	if len(sys.argv)==2:
		format_gaussian_input_from_xyz(sys.argv[1])
	if len(sys.argv)>2:
		for f in sys.argv[2:]:
			print extract_xyz_from_gaussian_out(f)
			format_gaussian_input_from_xyz(f)
