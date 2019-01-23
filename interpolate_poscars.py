def interpolate_poscars_keep_selective_flags(poscar1, poscar2, nimages, path=None):
	"""
	:param poscar1 and poscar2: Take poscar1 and poscar2 that correspond to the reagent and product of a chemical reaction
					Both poscar have:
										The same cell parameters
										The same number, type and atom orderning
										The same keywords cart selective
										Consistent flags T T T of F F F or T F T for respective atoms
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
		posdic['keywords'] = lines[6][0]
		for i, line in enumerate(lines[7:]):
			posdic[i] = [float(x) for x in line[:3]] + line[3:]
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
		typeprocess = "{0!s}\n"
		positionline = "  {0: .6f}  {1: .6f}  {2: .6f}   {3!s} {4!s} {5!s}\n"
		postring = ''.join([
			name.format(dic['name']),
			lattice_constant.format(dic['scale']),
			vector.format(*tuple(dic['v1'])),
			vector.format(*tuple(dic['v2'])),
			vector.format(*tuple(dic['v3'])),
			''.join([number_of_atoms.format(x) for x in dic['atmnum']])+'\n',
			typeprocess.format(dic['keywords'])
		])
		for i in range(len(dic.keys())-7):
			postring = postring + positionline.format(*tuple(dic[i]))
		return postring

	def build_interpolated_dics(posdic1, posdic2, nimages):
		"""
		:param posdic1 posdic2: two dics corresponding to poscar1 and poscar2 as returned by build_poscar_dic_representation
		:param nimages: the number of images
		:return: a list of dics corresponding to the images poscars, there will be nimages+2 dics with the first one corresponding to poscar1 and the last one corresponding to poscar2
		"""

		diclist= 'hahahahaha'
		return diclist


	posdic1 = build_poscar_dic_representation(poscar1)
	print build_poscar_from_dic_repr(posdic1)
	# posdic2 = build_poscar_dic_representation(poscar2)
	# diclist = build_interpolated_dics(posdic1, posdic2, nimages)
	# if path==None:
	# 	import os
	# 	path = os.getcwd()
	#
	# for i, dic in enumerate(diclist):
	# 	print i
	# 	s = build_poscar_from_dic_repr(dic)
	# 	with open(path+'0{}/POSCAR'.format(i), 'wb') as w:
	# 		w.write(s)



if __name__ == "__main__":
	import sys
	# interpolate_poscars_keep_selective_flags(sys.argv[1], sys.argv[2], sys.argv[3])
	interpolate_poscars_keep_selective_flags('/home/macenrola/Documents/vasp/NEB_mpXylene/mXylene-Protonated', '/home/macenrola/Documents/vasp/NEB_mpXylene/oXylene-Protonated',
											 8, '/home/macenrola/Documents/vasp/NEB_mpXylene/')
