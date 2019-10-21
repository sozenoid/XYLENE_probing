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

def plot_xy(x,y, name,f):
	"""
	pre: takes in x and y as lists of float
	post: prints the graph
	"""
	plt.close()
	plt.plot(x,y)
	plt.xlabel("COLVAR 1 (au)")
	plt.ylabel("energy (kcal/mol)")
	plt.title(name)
	plt.savefig(f+'.png')
	print '{} is done'.format(name)
	plt.show()

def plot_xyz(x,y,z, name):
	"""
	pre: takes in x, y, z as three vectors of float
	post: print the 3d graph
	"""
	dim = 100
	plt.close()
	x = np.asarray(x).reshape((dim, dim))
        y = np.asarray(y).reshape((dim, dim))
        z = np.asarray(z).reshape((dim, dim))
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(x,y,z,cmap=cm.coolwarm)
	ax.set_xlabel("coord nbr (final carbon)")
	ax.set_ylabel("coord nbr (initial carbon)")
	ax.set_zlabel("energy (kcal/mol)")
	plt.title(name+' min={0:2.2f} kcal/mol'.format(z.min()))
	plt.show()
	# print('{} has energy of {} kcal/mol'.format(name, z.min()))
	print(z.min())

def plot_colormap(x,y,z,name):
	"""
	pre: takes in x,y,z and will plot a colormap
	"""
	dim = 100
	x = np.asarray(x).reshape((dim, dim))
        y = np.asarray(y).reshape((dim, dim))
        z = np.asarray(z).reshape((dim, dim))
	levels = MaxNLocator(nbins=15).tick_values(z.min(), z.max())
	cmap = plt.get_cmap('PiYG')
	fig, ax = plt.subplots(constrained_layout=True)
	cf = ax.contourf(x,y,z,levels=levels, cmap=cmap)
	ax.set_xlabel("coord nbr (initial carbon)")
	ax.set_ylabel("coord nbr (final carbon)")
	cbar=fig.colorbar(cf, ax=ax)
	cbar.ax.set_ylabel('Free energy (kcal/mol)')
	plt.savefig(f+'.png')
	# plt.show()
	plt.close()


def compute_time(fname, thres=[0.3,0.7], T=300):
	"""
	PRE: Takes in the COLVAR output file of a metadynamic simulation
	POST: Will compute the accelerated time it took for the system's colvar
	      to pass the thres point
	      it is assumed that the number of thres elements fits the number of colvars, actually here 2 colvars are hardcoded
	      The energy is given in hartree and the temperature in Kelvin, here 300 by default, one hartree is 627.5 kcal/mol and
	      R=1.987E-3 kcal/K/mol
	      The timestep is expected to be 1 fs
	"""
	Ha2kcal = 627.5
	R=1.987E-3
	beta=(R*T)**-1
	with open(f, 'rb') as r:
		lines = [x.strip().split() for x in r.readlines()[:-1]]
#	for i in vars:
#		print i
	accelerated_time = 0
	current_time = 0
	for var in lines:
		try:
			var = [float(x) for x in var]
		except:
			print "file {} probably corrupted".format(f)
			break
		current_time=var[0]

		if len(thres)==2:
			crit=[var[1]<=thres[0], var[2]>=thres[1]]
			badcrit=[var[1]<=thres[0], var[2]<=thres[0]] # if both carbons have very low coordinations
			accelerated_time+=np.exp(beta*Ha2kcal*var[7])
		elif len(thres)==1:
			crit=[var[1]>=thres[0]]
			badcrit=[False]				     # it is a bad sign ghh
			accelerated_time+=np.exp(beta*Ha2kcal*var[4])
		if all(badcrit): 
			print "{} is bad at line {}".format(fname, var)
			break
		if all(crit):
			print "in the zone: MTD_time: {}; Real_equilvalent: {}".format(var[0], accelerated_time*1e-15)
			return current_time, accelerated_time*1e-15 # factor 1e-15 added to account for fs units
	return 0.0,0.0


def compute_avg_time_from_fit(timedist):
	"""
	PRE: Takes in a range of times assumed to belong to an exponential distribution
	POST: fits the experimental cdf of the timedist to 1-exp(-t/tau) to find tau
	"""
	sortedist = sorted(timedist)
	p = 1. * np.arange(len(sortedist)) / (len(sortedist) - 1)
	#
	def f_to_fit(x, scale):
		return 1-np.exp(-x/scale)
	#
	xn = np.linspace(sortedist[0], sortedist[-1], 200)
	p0 = np.median(sortedist)
	popt, pcov = optimize.curve_fit(f_to_fit, sortedist, p, p0=p0)
	#
	plt.plot(sortedist, p, 'or')
	plt.plot(xn, f_to_fit(xn, *popt))
	plt.show()

	return popt, sortedist, p

def KS_test(timedist, scale):
	"""
	PRE: a sequence of escape times from the botton energy well (computed from the MTD time and the bias potential)
	POST: Returns wether this time distribution follows a Poisson distribution AKA the law of rare events
	"""
	print stats.kstest(rvs=timedist, cdf='expon', args=(0,scale), N=len(timedist))

def reformat_all_dump(time_dist_dump_file):
	"""
	PRE:  A file formatted as follows 1000.0  28-prot-mxylene-MO-TS.inp-COLVAR.metadynLog     20840.0 13071211.9322
	      where the first colum is the temperature, the second a file name containing an indication of the system here MO
	      the third column is the MTD time and the last the real time
	POST: will produce  different files named after the temperature and system and containing only the real time
	"""
# =============================================================================
# 	system_markers=['MO-CB6', 'MO-CB7', 'MP-CB6', 'MP-CB7', 'MO-vac', 'MP-vac',
# 	'oxylene-cb6', 'pxylene-cb6', 'mxylene-cb6', 'oxylene-cb7', 'pxylene-cb7', 'mxylene-cb7',
# 	'adamantanol_cb7']
# =============================================================================
	system_markers=['oxylene-prot-cb6', 'mxylene-prot-cb6', 'pxylene-prot-cb6', 'oxylene-prot-cb7', 'mxylene-prot-cb7', 'pxylene-prot-cb7']
	systems = dict()
	with open(time_dist_dump_file, 'rb') as r: lines = [x.strip().split() for x in r.readlines()]
	#
	for els in lines:
		for mkr in system_markers:
			if mkr in els[1]:
				label = "{}-{}".format(mkr, str(els[0]))
				if label not in systems:
					systems[label] = [els]
				else:
					systems[label].append(els)
				break
	print systems.keys()

	for k in systems.keys():
		with open(time_dist_dump_file+"-{}-KS".format(k), "wb") as w:
			w.writelines([x[-1]+"\n" for x in systems[k]])

def print_average_energy_for_xyz_file(xyzfile):
	"""
	PRE  : Takes in a xyz trajectory
	POST : Computes the average energy and prints it, along with the number of snapshot
	"""
	E_list=[]
	with open(xyzfile, 'rb') as r:
		for line in r:
			if "E =" in line:
				E_list.append(float(line.split()[-1]))

	print "The energy is {} ({} snapshots)".format(sum(E_list)/len(E_list), len(E_list))

def plot_big_array_of_fes(flist):
	"""
	PRE : Takes a big array of FES txt files
	POST: Will plot a grid of plots with all the images
	"""
	nplots = len(flist)
	print nplots, nplots%10, nplots/10
	nrows=20
	flist = sorted(flist)
# =============================================================================
# 	fig, ax = plt.subplots(nrows, nplots/nrows, sharex='col', sharey='row')
# =============================================================================
	fig, axes = plt.subplots(nrows+1,nplots/nrows)

	cmap = plt.get_cmap('PiYG')
	dim = 100
	for i, ax in enumerate(axes.reshape(-1)):
		if i==nplots:
			break
		name = '-'.join(flist[i].split('/')[-4:])
		print i, name
		x,y,z=get_xyz(flist[i])
		title = [name.split('-')[0]] + name.split('-')[3:6]
		x = np.asarray(x).reshape((dim, dim))
		y = np.asarray(y).reshape((dim, dim))
		z = np.asarray(z).reshape((dim, dim))
		levels = MaxNLocator(nbins=15).tick_values(z.min(), z.max())
		cf = ax.contourf(x,y,z, cmap=cmap, levels=levels)
		ax.set(aspect='equal')
		ax.set_title('-'.join(title))
		ax.get_xaxis().set_major_formatter(plt.NullFormatter())
		ax.get_yaxis().set_major_formatter(plt.NullFormatter())
		ax.title.set_fontsize(3)
		ax.title.set_position([.5, 0.8])
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.tight_layout()
	plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False,
	labeltop=False) # labels along the bottom edge are off
	plt.show()
	
	
	
if __name__ == "__main__":
	import glob
	from matplotlib import cm
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib.colors import BoundaryNorm
	from matplotlib.ticker import MaxNLocator
	import sys
	from scipy import stats, optimize
	import numpy as np
	import os
	cwd = os.getcwd()
	if len(sys.argv)==1:
		print """
==========NO ARGUMENT===========
Please use 1D or 2D fes data produced by graph.popt and ending in '.txt'
of a 2 colvar output file from a metadynamic run ending in 'Log'
pwd = {}""".format(cwd)
		plot_big_array_of_fes(glob.glob("/home/macenrola/Documents/XYLENE/inputs/for-reaction-flexible-cb/DUMP_SLOW_REACT/fes_validation/*00/*-*/MO-CB6/*txt"))
	elif len(sys.argv)==2:
		flist=[sys.argv[1]]
	else:
		flist=sys.argv[1:]
	MP=False
	dumpfile='{}/DUMP_MTS'.format(cwd)
        with open(dumpfile, 'wb'):
		print("writing to {}".format(dumpfile))
	timedist =[]
	for f in sorted(flist):
		name = '-'.join(f.split('/')[-2:])
		# print name
		if name[-4:]==".txt":
			X =get_xyz(f)
			if 'MP' in name and not MP:
				print 'MP now'
				MP=True
			if len(X) ==2:
				x,y=X
				plot_xy([i for i in x], [k*627.5 for k in y], name,f)
			else:
				x,y,z=X
				# print f
				plot_colormap(x,y, [k*627.5 for k in z],name)

		if name[-3:]=='Log':
			print name
			T=float(f.split('/')[0])
			with open(dumpfile,'ab') as a:
				ctime,atime=compute_time(f,thres=[10], T=T)
				if ctime > 300:
					a.write("{}\t{}\t{}\t{}\n".format(T,name,ctime, atime))
					timedist.append(atime)
				else:
					print "accelerated time out of bounds ({})".format(ctime)
		if name[-2:]=="KS":
			with open(f, 'rb') as r: lines=sorted([float(x.strip()) for x in r.readlines()])[:]
			# print lines
			#scale=float(name.split('-')[-2])
			try:
				fit_scale, sortedist, p =compute_avg_time_from_fit(lines)
                	        print "{}-{}".format(name, fit_scale[0])
				with open(f+"_for_plot_{}".format(int(fit_scale[0])), 'wb') as w:
					w.writelines("\n".join(["{} {}".format(x[0], x[1]) for x in zip(sortedist,p)]))
        	                #print lines
			except : print "impossible to fit and or plot for {}; lines are {}".format(name, lines)

			#scale=np.median(lines)/np.log(2)
			KS_test(lines, fit_scale)

		if name[-4:]=="dump":
			reformat_all_dump(f)

		if name[-4:]==".xyz":
			print name
			print_average_energy_for_xyz_file(f)
