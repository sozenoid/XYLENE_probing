def parse_z_matrix(z_matrix):
	"""
	PRE: Takes in a z matrix as produced by obabel as a text block
	POST: returns a dictionary representation of the matrix
		one key holds the structure ("struct") the other keys are named after the variable names that are expected to be consistent

	"""
	if len(z_matrix.split())<3:
		return None

	zmatdic={}
	lines=z_matrix.split('\n')
	pastvariables=False
	for l in lines:
		if "Variables:" in l:
			pastvariables=True
		if pastvariables and '=' in l:
			temp_line=l.strip().split('=')
			zmatdic[temp_line[0]]=temp_line[1].strip()
		if 'structure' not in zmatdic and not pastvariables:
			zmatdic['structure']=[l]
		elif not pastvariables:
			zmatdic['structure'].append(l)
	zmatdic["structure"]="\n".join(zmatdic["structure"])
	return zmatdic

def average_z_matrix(list_of_z_matrix_as_dic):
	"""
	PRE: Takes a list of dics representing z matrices
	POST: Returns a single dic with all values averaged and another dic with the deviations
	"""
	example_dic=dict.fromkeys(list_of_z_matrix_as_dic[0].keys())
	clustered_dic={}
	avg_clustered_dics=[]
	##### SPLIT THE ORIGINAL LIST BY CLUSTER NUMBER
	for dic in list_of_z_matrix_as_dic:
		if dic["cluster"] not in clustered_dic:
			clustered_dic[dic["cluster"]] = [dic]
		else:
			clustered_dic[dic["cluster"]].append(dic)

	#### COMPUTE THE AVERAGE Z MATRIX CLUSTER WISE
	print [len(x) for x in clustered_dic.values()]
	for c in clustered_dic:
		res_dic = dict.fromkeys(clustered_dic[c][0].keys())
		res_dic["structure"] = clustered_dic[c][0]["structure"]

		for k in res_dic:
			temp_value=[]
			if k=='structure': continue
			for dic in clustered_dic[c]:
				temp_value.append(float(dic[k]))

			# print k,temp_value, statistics.stdev(temp_value)
			# temp_value=[x if x<180 else abs(x-360) for x in temp_value]
			# print k,temp_value, statistics.stdev(temp_value)

			res_dic[k]=1.0*sum(temp_value)/len(temp_value)
			# if res_dic[k]<0: res_dic[k]=res_dic[k]+360
		avg_clustered_dics.append(res_dic)
	return avg_clustered_dics

def write_z_matrix_from_dic(z_matrix_dic, fname):
	"""
	PRE: takes in a z matrix written as dic
	POST: will produce a file formatted as obabel z matrix
	"""
	with open(fname, 'wb') as w:
		w.write(z_matrix_dic["structure"]+"\n")
		w.write("Variables:\n")
		for k in sorted(z_matrix_dic.keys()):
			if k =="structure":
				continue
			w.write("{}={}\n".format(k, z_matrix_dic[k]))

def cluster_list_of_zdics(list_of_z_matrix_as_dic):
	"""
	PRE: Takes in a list of z matrices as cluster_list_of_zdics
	POST: will use PCA or kmeans to cluster the coordinates.
	"""
	full_list=[]
	for dic in list_of_z_matrix_as_dic:
		tmplist=[]
		for k in sorted(dic.keys()):
			if k=="structure":
				continue
			tmplist.append(dic[k])
		full_list.append(tmplist)
		# print tmplist
	full_list_array= np.array([np.array(x) for x in full_list])
	# print(full_list_array.shape)

	######## COMPUTE THE PCA
	# pca=PCA(n_components=2)
	# transformed_coords =  pca.fit_transform(full_list_array)
	# transformed_coords_pd = (pd.DataFrame(transformed_coords))
	# transformed_coords_pd.columns = ["PC1" ,"PC2"]

	####### PLOT AS A SANITY CHECK
	# ax = transformed_coords_pd.plot(kind='scatter', x='PC2', y='PC1', figsize=(16,8))
	# plt.show()
	# print transformed_coords
	# print(pca.explained_variance_ratio_)
	# plt.scatter(transformed_coords.T[0], transformed_coords.T[1])
	# plt.show()

	#### COMPUTE THE KMEANS AND BUILD THE ESTIMATOR
	kmeans = KMeans(n_clusters=100)
	clusters = kmeans.fit_predict(full_list_array)

	#### PLOT AS A SANITY CHECK
	# transformed_coords_pd['cluster'] = pd.Series(clusters.labels_, index=transformed_coords_pd.index)
	# axk =transformed_coords_pd.plot(
    # kind='scatter',
    # x='PC2',y='PC1',
    # c=transformed_coords_pd.cluster.astype(np.float),
    # figsize=(16,8))
	# plt.show()

	####### ADD A KEY TO EACH DICTIONARY INDICATING WHICH CLUSTER
	for i, dic in enumerate(list_of_z_matrix_as_dic):
		dic["cluster"] = clusters[i]
	return list_of_z_matrix_as_dic

def process_z_matrix_trajectory(zmatrix_traj_file):
	"""
	PRE : takes in a zmatrix_traj_file
	POST: will extract all the indivudual blocks and process them using parse_z_matrix
	"""
	temp_z=[]
	zmatlist=[]
	current_snapshot=0
	with open(zmatrix_traj_file, 'rb') as r:
		for i,lines in enumerate(r):
			# if i==100000: break
			if '#' in lines:
				current_snapshot+=1
				if current_snapshot%100==0:
					print "current snapshot is {}".format(current_snapshot)
				if temp_z!=[]:
					temp_z_mat=parse_z_matrix(''.join(temp_z))
					if temp_z_mat is not None:
						zmatlist.append(temp_z_mat)
				temp_z=[]
				temp_z.append(lines)
			else:
				temp_z.append(lines)

		if temp_z!=[]:
			# temp_z.append(lines)
			temp_z_mat=parse_z_matrix(''.join(temp_z))
			if temp_z_mat != []:
				if zmatlist!=[]:
					nkeys = len(zmatlist[-1].keys())
					if nkeys!=len(temp_z_mat.keys()):
						print "corrupted file, expected {} keys but found {}.\n Total indivudual snapshots found {}".format(nkeys, len(temp_z_mat.keys()), len(zmatlist))
					else:
						zmatlist.append(temp_z_mat)
			temp_z=[]

	clustered_zmatlist = cluster_list_of_zdics(zmatlist)
	########## PRINTING THOSE WINNERS
	avg_z_matlist = average_z_matrix(clustered_zmatlist) # will contain as many matrices as there are clusters
	print len(avg_z_matlist)
	for i, avg_z_mat in enumerate(avg_z_matlist):
		write_z_matrix_from_dic(avg_z_mat, "avg_resdic-{}.gzmat".format(i))

	return

if __name__ == "__main__":
	import statistics
	import numpy as np
	import numpy.core.multiarray
	from sklearn.decomposition import PCA
	from sklearn.cluster import KMeans
	import matplotlib.pyplot as plt
	import pandas as pd

	with open('huhhu.gzmat') as r:
		lines=r.readlines()
	process_z_matrix_trajectory('cb6.inp-pos-1-aligned.gzmat')
