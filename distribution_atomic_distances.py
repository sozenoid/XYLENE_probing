##########
#
#
# Find distribution of atomic distances from a xyz file 
#
#########

def read_conformations_from_xyz_file(xyz_file):
	"""
	PRE: takes an xyz file that represents a complex
	POST: will return the distribution of non bonded atom-atom distances between host and guest
		as a dictionary
	"""
	with open(xyz_file, "rb") as r:
		tempconf = []
		dictconf = {}
		indexconf= -1
		for line in r:
			if "i =" in line:
				dictconf[int(indexconf)]=tempconf[:-1]
				tempconf=[]
				indexconf=line.split()[2][:-1]
				continue
			tempconf.append(line.strip())
	return dictconf

def compute_distance_distribution_for_set_of_configurations(conf_dict):
	"""
	PRE  : Takes in dictionary of conformations 
	POST : will return the cdf and pdf of the distances
	"""
	indexesMol1=range(26)
	indexesMol2=range(27,152)
	dist_dic = {}
	for k in sorted(conf_dict):
		if k>1001: break
		mol1atmpos=[]
		mol2atmpos=[]
		#print k
		for pos in range(len(conf_dict[k])):
			parts=conf_dict[k][pos].strip().split()
			if pos in indexesMol1:
				mol1atmpos.append(parts)
			else:
				mol2atmpos.append(parts)

		for xyz1 in mol1atmpos:
			for xyz2 in mol2atmpos:
				dist=sum([(float(x[0])-float(x[1]))**2 for x in zip(xyz1[1:], xyz2[1:])])**.5
				atm_type="{}-{}".format(xyz1[0], xyz2[0])
				#print atm_type, dist
				if atm_type in dist_dic:
					dist_dic[atm_type].append(dist)
				else:
					dist_dic[atm_type] = [dist]
	return dist_dic

def extract_dists_from_many_files(flist):
	"""
	PRE: Takes in a list of files
	POST: Will extract the atom-atom distances for each of them and compile a dictionary of distances by atom types
	"""
	dist_dic={}
	for f in flist:
		temp_conf_dict=read_conformations_from_xyz_file(f)
		temp_dist_dict=compute_distance_distribution_for_set_of_configurations(temp_conf_dict)
		for k in temp_dist_dict:
			if k in dist_dic:
				dist_dic[k].extend(temp_dist_dict[k])
			else:
				dist_dic[k]=temp_dist_dict[k]
	return dist_dic

def plot_dist_dic(dist_dic):
	"""
	PRE: Takes in a dictionary of distances for different atom types
	POST: Will print the pdf of the distances for each atom types and then together
	"""
	for k in dist_dic:
		tempdist = dist_dic[k]
		x=[x/float(len(tempdist)) for x in range(len(tempdist))]
		sx = sorted(tempdist)
		bins = np.linspace(sx[0], sx[-1], 50)
		#bins=[1,2,3,4,5,6,7,8,9]
		hist, hist_bins = np.histogram(sx, bins, density=True)
		centered_hist_bins=[a+(hist_bins[1]-hist_bins[0])/2 for a in hist_bins[:-1]]
		#sxcentered = [a+(sx[1]-sx[0]) for a in sx[:-1]]
		#f=interpolate.interp1d(sx, x, 'cubic')
		#sx_low=np.linspace(sx[0], sx[-1], 100)
		#dxdd = np.diff(np.array(f(sx_low)))
		plt.plot(centered_hist_bins, hist)
		#plt.plot(sx,x)
		plt.title(k)
		plt.show()


def merge_dics(old_dic, new_dic):
	"""
	PRE: Takes an old dic and a new dic
	POST: Will copy the new keys to the old dic along with the values if the key does not exist
	      If the key is already present in the old dic, the old value will be extended with the new value
	      The function returns nothing and updates the old dic
	"""
	for nk in new_dic:
		if nk not in old_dic:
			old_dic[nk] = new_dic[nk]
		else:
			old_dic[nk].extend(new_dic[nk])
	return

if __name__ == "__main__":
	import matplotlib.pyplot as plt
	import numpy as np
	from scipy import interpolate
	import cPickle
	import glob
	testfile="/home/macenrola/Documents/XYLENE/pm6-mtd/double-coords-vdw/popping/sanity-check-adamantane/test/8-adamantanol_cb7.inp-pos-1.xyz"
	print "haha"
	flist=glob.glob("/home/macenrola/Documents/XYLENE/pm6-mtd/double-coords-vdw/popping/sanity-check-adamantane/300/*.inp-pos-1.xyz")
	all_confs_dict={}
	flist=['/home/macenrola/Documents/XYLENE/pm6-mtd/double-coords-vdw/popping/sanity-check-adamantane/opt_adamantanol/last_iter_adamantanol@cb7.xyz']
	for i, f in enumerate(flist):
		print "{}/{}".format(i, len(flist)), f, "  ".join(["key:{} ({})".format(k, len(all_confs_dict[k])) for k in all_confs_dict])
		dict=read_conformations_from_xyz_file(testfile)
		dist_dic=compute_distance_distribution_for_set_of_configurations(dict)
		merge_dics(all_confs_dict, dist_dic)

	plot_dist_dic(all_confs_dict)
