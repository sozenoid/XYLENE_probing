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
	return subprocess.call(cmd.format(traj_name), shell=False)


def print_trajectory(traj_name):
	cmd = "/home/uccahcl/cp2k/exe/local/cp2k.sopt -i {0} -o {0}.out"
	print cmd.format(traj_name)

if __name__ == '__main__':
	pool = multiprocessing.Pool(None)
	tasks = sys.argv[1:]
	results = []
	r = pool.map_async(print_trajectory, tasks, callback=results.append)
	r.wait() # Wait on the results
	print results
	r = pool.map_async(launch_trajectory, tasks, callback=results.append)
	r.wait() # Wait on the results
	print results