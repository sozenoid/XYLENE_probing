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
import os
def launch_trajectory(traj_name):
	cmd = "cp2k.sopt -i {0} -o {0}.out"
	return subprocess.call(cmd.format(traj_name).split(), shell=False)


def print_trajectory(traj_name):
	cmd = "cp2k.sopt -i {0} -o {0}.out"
	print cmd.format(traj_name)
	
def make_inputs(pattern_to_target, node_size=40):
	"""
	PRE: Will take a suffix, collect all the matching files and generate inputs files for this script
	POSE: Will print a nfile/node_size files
	"""
	# This is the Thomas pattern
	pattern="""#!/bin/bash

#PBS -N 350-adam
#PBS -l select=1:ncpus=40:mpiprocs=40:ompthreads=1
#PBS -l walltime=48:0:00
#PBS -o 350-adam.out
#PBS -e 350-adam.err
#PBS -l software="CP2K"

cd {}
export OMP_NUM_THREADS=1
# Now run the program
export CP2K_DATA_DIR=/home/users/app/common_apps/chemistry/cp2k/dev/8.0_Jan2020_gcc9_mkl/data/
# module load knl-software intel-mpi cp2k/8.0-dev
module load .legacy-stratus-software  knl-software  intel-mkl cp2k/8.0-dev
/home/users/astar/ihpc/stuhlamb/XYLENE_probing/trajectory_manager.py {}
"""
	flist = glob.glob("*{}".format(pattern_to_target))
	ftowrite = []
	for i, f in enumerate(sorted(flist)):
		ftowrite.append(f)
		if (i+1)%node_size==0 and i!=0:
			with open("traj_launcher_{}.pbs".format(i/node_size), "wt") as w:
				w.write(pattern.format(os.getcwd()," ".join(ftowrite)))
				ftowrite=[]
	if ftowrite!=[]:
		with open("traj_launcher_{}.pbs".format(i/node_size), "wt") as w:
			w.write(pattern.format(os.getcwd()," ".join(ftowrite)))

if __name__ == '__main__':
	pool = multiprocessing.Pool(None)
	if len(sys.argv)==1:
		print """You need to provide arguments
		Use `trajectory_manager.py suffix inp` to generate input files for all files in the current directory that have the suffix "inp"
		"""
	elif sys.argv[1]=="suffix" and len(sys.argv)==3:
		make_inputs(sys.argv[2])
	else:
		# pool = multiprocessing.Pool(len(sys.argv[1:]))
		tasks = sys.argv[1:]
		results = []
		r = pool.map_async(print_trajectory, tasks, callback=results.append)
		r.wait() # Wait on the results
		print results
		r = pool.map(launch_trajectory, tasks)
		print r
