def interpolate_poscars_keep_selective_flags(poscar1, poscar2, nimages, path=None):
	"""
	:param poscar1 and poscar2: Take poscar1 and poscar2 that correspond to the reagent and product of a chemical reaction
					Both poscar have:
										The same cell parameters
										The same number, type and atom orderning
										The same keywords cart selective
										Consistent flags T T T of F F F or T F T for respective atoms, only the flags, vectors and keywords from poscar1 will be used though (overwriting poscar2 in the last file even)
	:nimages: The number of images, poscar1 and poscar2 are not included
	:return: a POSCAR file in each of the files; for example for 3 images; there will be 00 01 02 03 04 files with poscar1 in 00 and poscar2 in 04; the files 01 02 03 containt the POSCARs of the images
	"""

	def build_poscar_dic_representation(poscar):
		"""
		:param poscar: take in a poscar file
		:return: a dictionary the keys are name, scale, v1, v2, v3, atmnum, keywords, and numbers for the atomic positions
		"""
		posdic = {}
		with open(poscar, 'rb') as r:
			lines = [x.strip().split() for x in r.readlines()]
		print lines
		posdic['name'] = lines[0][0]
		posdic['scale'] = float(lines[1][0])
		posdic['v1'] = [float(x) for x in lines[2]]
		posdic['v2'] = [float(x) for x in lines[3]]
		posdic['v3'] = [float(x) for x in lines[4]]
		posdic['atmnum'] = [int(x) for x in lines[5]]
		posdic['keywords'] = lines[6]
		for i, line in enumerate(lines[7:]):
			posdic[i] = [float(x) for x in line[:3]] + line[3:]
		recenter_dic_coords(posdic)
		return posdic

	def build_poscar_from_dic_repr(dic):
		"""
		:param dic: from a dic  produced by build_poscar_dic_representation
		:return: a string that correspond to a poscar
		"""
		print dic
		name = "{0!s}\n"
		lattice_constant = "{0:.5f}\n"
		vector = '     {0:.6f}    {1:.6f}    {2:.6f}\n'
		number_of_atoms = "     {0:d}"
		typeprocess = "{0!s} "
		positionline = "  {0: .6f}  {1: .6f}  {2: .6f}   {3!s} {4!s} {5!s}\n"
		postring = ''.join([
			name.format(dic['name']),
			lattice_constant.format(dic['scale']),
			vector.format(*tuple(dic['v1'])),
			vector.format(*tuple(dic['v2'])),
			vector.format(*tuple(dic['v3'])),
			''.join([number_of_atoms.format(x) for x in dic['atmnum']])+'\n',
			''.join([typeprocess.format(x) for x in dic['keywords']]) + '\n'
		])
		for i in range(len(dic.keys())-7):
			postring = postring + positionline.format(*tuple(dic[i]))
		return postring

	def recenter_dic_coords(dic):
		"""
		:param dic: Takes in a dic as provided by the function  build_poscar_dic_representation
		:return: will modify its cartesian coordinates so that its center appears at the middle of the cell (otherwise visualization becomes very confusing)
		"""
		x_center_ref, y_center_ref, z_center_ref = dic['v1'][0]*dic['scale']/2., dic['v2'][1]*dic['scale']/2., dic['v3'][2]*dic['scale']/2.
		x_center, y_center, z_center = [],[],[]
		for i in range(len(dic.keys())-7):
			cx,cy,cz = dic[i][:3]
			x_center.append(cx)
			y_center.append(cy)
			z_center.append(cz)

		xc = sum(x_center)*1.0/len(x_center)
		yc = sum(y_center) * 1.0 / len(y_center)
		zc = sum(z_center) * 1.0 / len(z_center)

		print xc, yc, zc

		for i in range(len(dic.keys())-7):
			tc = dic[i][:3]
			dic[i][:3] = [x[0]+x[1]-x[2] for x in zip(tc, [x_center_ref, y_center_ref, z_center_ref],[xc, yc, zc])]

	def build_interpolated_dics(posdic1, posdic2, nimages):
		"""
		:param posdic1 posdic2: two dics corresponding to poscar1 and poscar2 as returned by build_poscar_dic_representation
		:param nimages: the number of images
		:return: a list of dics corresponding to the images poscars, there will be nimages+2 dics with the first one corresponding to poscar1 and the last one corresponding to poscar2
		"""
		import numpy as np
		diclist= []
		for k in range(nimages+2):# Creates the diclist
			if k==0:
				diclist.append(posdic1)
			# elif k==nimages+1:
			# 	diclist.append(posdic2)
			# 	pass
			else:
				diclist.append({'name':posdic1['name']+'->'+posdic2['name'],
								'scale':posdic1['scale'],
								'v1':posdic1['v1'],
								'v2': posdic1['v2'],
								'v3': posdic1['v3'],
								'atmnum':posdic1['atmnum'],
								'keywords':posdic1['keywords']
								})

		for atom_line in range(len(posdic1.keys())-7): # will iterate over the atom numbers (there are 6 fields that are not line numbers)
			i = atom_line
			tmp_pos1 = posdic1[i]
			tmp_pos2 = posdic2[i]
			flags1 = tmp_pos1[3:]
			x_range = np.linspace(tmp_pos1[0], tmp_pos2[0], nimages+2)
			y_range = np.linspace(tmp_pos1[1], tmp_pos2[1], nimages+2)
			z_range = np.linspace(tmp_pos1[2], tmp_pos2[2], nimages+2)
			for j, intermediate_point_coords in enumerate(zip(x_range, y_range, z_range)):
				if j==0:
					continue
				elif j==nimages+1:
					pass
				diclist[j][i] = list(intermediate_point_coords) + flags1
		return diclist


	posdic1 = build_poscar_dic_representation(poscar1)
	posdic2 = build_poscar_dic_representation(poscar2)

	interpolated_dics = build_interpolated_dics(posdic1, posdic2, nimages)

	# if path==None:
	# 	import os
	# 	path = os.getcwd()
	#
	for i, dic in enumerate(interpolated_dics):
		s = build_poscar_from_dic_repr(dic)
		with open(path+'0{}/POSCAR'.format(i), 'wb') as w:
			w.write(s)


def recenter_direct_poscar(poscar):
	"""
	:param poscar: Takes in a direct poscar and recenters it
	:return: another file named poscar_recentered that will contain the recentered data
	"""
	with open(poscar, 'rb') as r:
		with open(poscar+'_recentered', 'wb') as w:
			for i, line in enumerate(r):
				if i<8 or line==' \n':
					w.write(line)
				else:
					coords = [float(x) for x in line.strip().split()]
					coords_center = [x if abs(x)<0.5 else x-1.0 for x in coords]
					w.write('  {0: .16f}  {1: .16f}  {2: .16f}\n'.format(*coords_center))


if __name__ == "__main__":
	import sys
	# interpolate_poscars_keep_selective_flags(sys.argv[1], sys.argv[2], sys.argv[3])
	# interpolate_poscars_keep_selective_flags('/home/macenrola/Documents/vasp/NEB_mpXylene/mXylene-Protonated', '/home/macenrola/Documents/vasp/NEB_mpXylene/pXylene-Protonated',
	# 										 8, '/home/macenrola/Documents/vasp/NEB_mpXylene/')
	recenter_direct_poscar('/home/macenrola/Documents/vasp/xylene/contcar')