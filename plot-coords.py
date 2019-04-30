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
	plt.xlabel("Angle (degree)")
	plt.ylabel("energy (kcal/mol)")
	plt.title(name)
	plt.savefig(f+'.png')
	print '{} is done'.format(name)
	# plt.show()

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


def compute_time(fname, thres=[0.2,0.8], T=300):
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
		vars = [[float(y) for y in x.strip().split()] for x in r.readlines()]
#	for i in vars:
#		print i
	accelerated_time = 0
	current_time = 0
	for var in vars:
		current_time=var[0]
		crit=[var[1]<=thres[0], var[2]>=thres[1]]
		accelerated_time+=np.exp(beta*Ha2kcal*var[7])
		if all(crit):
			print "in the zone: MTD_time: {}; Real_equilvalent: {}".format(var[0], accelerated_time)
			break
	return current_time, accelerated_time

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

	return popt

def KS_test(timedist, scale):
	"""
	PRE: a sequence of escape times from the botton energy well (computed from the MTD time and the bias potential)
	POST: Returns wether this time distribution follows a Poisson distribution AKA the law of rare events
	"""
	print stats.kstest(rvs=timedist, cdf='expon', args=(0,scale), N=1000)

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
	if len(sys.argv)==1:
		print """
==========NO ARGUMENT===========
Please use 1D or 2D fes data produced by graph.popt and ending in '.txt'
of a 2 colvar output file from a metadynamic run ending in 'Log'
"""
	elif len(sys.argv)==2:
		flist=[sys.argv[1]]
	else:
		flist=sys.argv[1:]
	MP=False
	dumpfile='/home/macenrola/Documents/XYLENE/pm6-mtd/double-coords-vdw/prod/300k-logs/MOCB7-DUMP'
        with open(dumpfile, 'wb'): pass
	timedist =[]
	for f in sorted(flist):
		name = '-'.join(f.split('/')[-2:])
		print name
		if name[-4:]==".txt":
			X =get_xyz(f)
			if 'MP' in name and not MP:
				print 'MP now'
				MP=True
			if len(X) ==2:
				x,y=X
				plot_xy([i*180/3.1316 for i in x], [k*627.5 for k in y], name,f)
			else:
				x,y,z=X
				# print f
				plot_colormap(x,y, [k*627.5 for k in z],name)

		if name[-3:]=='Log':
			print name
			with open(dumpfile,'ab') as a:
				ctime,atime=compute_time(f)
				if atime<=1E35:
					a.write("{}\t{}\t{}\n".format(name,ctime, atime))
					timedist.append(atime)
				else:
					print "accelerated time is over 1E35 ({})".format(atime)
		if name[-2:]=="KS":
			with open(f, 'rb') as r: lines=[float(x.strip()) for x in r.readlines()]
			#scale=float(name.split('-')[-2])
			fit_scale=compute_avg_time_from_fit(lines)
			print "fit scale is:{}".format(fit_scale)
			scale=np.median(lines)/np.log(2)
			KS_test(lines, fit_scale)
