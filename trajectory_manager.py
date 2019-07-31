#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 09:08:16 2019

@author: macenrola
"""

import multiprocessing
import subprocess
import sys
import glob

def launch_trajectory(traj_name):
	cmd = "/home/uccahcl/cp2k/exe/local/cp2k.sopt -i {0} -o {0}.out"
	return subprocess.call(cmd.format(traj_name).split(), shell=False)


def print_trajectory(traj_name):
	cmd = "/home/uccahcl/cp2k/exe/local/cp2k.sopt -i {0} -o {0}.out"
	print cmd.format(traj_name)
	
def make_inputs(suffix_to_target, node_size=24):
	"""
	PRE: Will take a suffix, collect all the matching files and generate inputs files for this script
	POSE: Will print a nfile/node_size files
	"""
	# This is the Thomas pattern
	pattern="""#!/bin/bash -l 
#$ -S /bin/bash
#$ -l h_rt=10:0:00
#$ -l mem=5G
#$ -N mtd
#$ -pe smp 24
#$ -cwd
#$ -P Free
#$ -A UCL_chemM_Lee
/home/uccahcl/XYLENE_probing/trajectory_manager.py {}
"""
	flist = glob.glob("*{}".format(suffix_to_target))
	ftowrite = []
	for i, f in enumerate(sorted(flist)):
		ftowrite.append(f)
		if (i%node_size==0 or i==len(flist)-1) and i!=0:
			with open("traj_launcher_{}.sh".format(i/node_size), "wt") as w:
				w.write(pattern.format(" ".join(ftowrite)))
				ftowrite=[]

if __name__ == '__main__':
	pool = multiprocessing.Pool(None)
	if len(sys.argv)==1:
		print """You need to provide arguments
		Use `trajectory_manager.py suffix inp` to generate input files 
		"""
	elif sys.argv[1]=="suffix" and len(sys.argv)==3:
		make_inputs(sys.argv[2])
	else:
		tasks = sys.argv[1:]
		results = []
		r = pool.map_async(print_trajectory, tasks, callback=results.append)
		r.wait() # Wait on the results
		print results
		r = pool.map(launch_trajectory, tasks)
		print r