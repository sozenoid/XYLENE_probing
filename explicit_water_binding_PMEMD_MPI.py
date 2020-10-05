import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
from subprocess import call
import os
# sys.path.append('/home/macenrola/Thesis/Molecule_operations/')
working_directory = '/home/uccahcl/Scratch/AMBER/'
bin_path_mpi = '/shared/ucl/apps/amber/amber-16-mpi/bin/'
bin_path_serial = '/shared/ucl/apps/amber/amber-16/bin/'
#from mol_ops_amber import make_pdb_complex_with_named_residues, align_mol

def create_topology_file(RDKIT_BLOCK, MPI=True):
	"""
	PRE: Takes in a RDKIT_BLOCK (sdf file in text) with valid 3D coordinates. A minimum energy conformation is needed.
		 ...The Coordinates main axis (defined by PCA) are expected to be aligned with the z axis and centered around 0 (for fitting with the CB7). 
		 The proper modules are loaded (the AMBER binaries are in the path and amber.sh has been sourced)
		 This script needs to be called AFTER sourcing of the /home/macenrola/Thesis/AMBER/amber16/amber.sh file, like this for example:
		 (my-rdkit-env) macenrola@macenrola-MacBookPro:~/Thesis/AMBER/converge_with_amber$ source /home/macenrola/Thesis/AMBER/amber16/amber.sh; python AMBER_CONVERGE.py
		 IF MPI==True, the topology files will NOT be created. The parallel version of amber doesn't have all the fancy antechamber and all
		 IF MPI==False, amber can't run in parallel but all the file preparation can be performed. We suggest running in MPI==False mode to prepare the files by loading a serial version of the amber module then going to MPI for the production
	POST: The prmtop and rst7 files are produced
	NOTE: A brief annealing could be great to estimate the best conformation of the host inside CB7
	"""

    ###### CREATION OF TOPOLOGY FILES
	
	RDKIT_BLOCK = align_mol(RDKIT_BLOCK)
	mol = Chem.MolFromMolBlock(RDKIT_BLOCK, removeHs = False)

	path = working_directory
	guest_name =  mol.GetProp('_Name')# Diadamantyl diammonium
	host_name = 'CB'
	fguest = path+guest_name+'.pdb'
	fCB = path+'{}.pdb'.format(host_name)
	fcomplex = path+guest_name+'-{}.pdb'.format(host_name)
	script_name_solv1 = 'make_prmtop_rst_first_solvation.sh'
	script_name_solv2 = 'make_prmtop_rst_second_solvation.sh'
	script_name_min_params = 'make_min_params.sh'
	script_name_min_run = 'run_min_algorithm.sh'

	
	if MPI != True:
		make_pdb_complex_with_named_residues(RDKIT_BLOCK, fguest, fCB, fcomplex) 
		CBmol2_file = 'CB.mol2'
		water_file = 'water.pdb'
		with open(path+script_name_solv1, 'wb') as w:
			w.write(get_prmtop_rst_script_solvation(charge=Chem.GetFormalCharge(mol)))
			call('chmod +x {}'.format(path+script_name_solv1), shell=True)
		with open(path+CBmol2_file, 'wb') as w:
			w.write(get_CB_mol2())
		with open(path+water_file, 'wb') as w:
			w.write(get_water_pdb())	

		call('{} {} {}'.format(path+script_name_solv1, guest_name, host_name), shell=True) # First attempt at solvation to know how many water molecules to remove
																				 
		watermol_to_remove = get_water_number_to_remove(path+'water-tleap.log', path+'{}-tleap.log'.format(guest_name), path+'{}-tleap.log'.format(host_name), path+'{}-{}-tleap.log'.format(guest_name, host_name))
		print watermol_to_remove
		with open(path+script_name_solv2, 'wb') as w:
			w.write(get_prmtop_rst_script_solvation(charge=Chem.GetFormalCharge(mol), watermol_to_remove=watermol_to_remove))
			call('chmod +x {}'.format(path+script_name_solv2), shell=True)

		call('{} {} {}'.format(path+script_name_solv2, guest_name, host_name), shell=True) # Second attempt at solvation with the proper number of water (1500 total)

		for rstfiles in [path+'{}.rst7'.format(x) for x in ['water', guest_name, host_name, '{}-{}'.format(guest_name, host_name)]]: # Correct the edges of the cubes
			correct_rst_files_for_cuboid(rstfiles)

		with open(path+script_name_min_params, 'wb') as w:
			w.write(get_minimization_route_script())
			call('chmod +x {}'.format(path+script_name_min_params), shell=True)
		call(path+script_name_min_params, shell=True)

		with open(path+script_name_min_run, 'wb') as w:
			w.write(get_minimization_run_script('water', guest_name, host_name, '{}-{}'.format(guest_name, host_name), bin_path_mpi))
			call('chmod +x {}'.format(path+script_name_min_run), shell=True)
	
	#### RUN THE ACTUAL SCRIPT FOR AVERAGING PURPOSES
	call(path+script_name_min_run, shell=True)



	###### PREPARE THE SYSTEM AND ADDS THE WATER
	### MINIMIZE TO RELAX BAD CONTACTS

def get_water_number_to_remove(waterLog, guestLog, CbLog, complexLog):
	"""
	PRE: Given the logs of previous tleap runs, 
	POST: It returns the number of water molecules that need to be removed in order to reach 1500 units for each of the 4 files, water, guest, cb, complex
	"""
	logs = [waterLog, guestLog, CbLog, complexLog]
	waters_to_remove = []

	for log in logs:
		with open(log, 'rb') as r:
			for line in r:
				if 'WAT' in line:
					waters_to_remove.append(line.split()[-1])
	return tuple([int(x)-1500 for x in waters_to_remove])

def correct_rst_files_for_cuboid(rstfile, maxdim =-1, suffix='_iso.rst7'):
	"""
	PRE: A RST7 file is provided
	POST: The last line of the RST7 file will be modified so that the box now is a perfect cube. Expends the smallest edges of the cube so that they match the largest edge. According to amber tutorial21
		  http://ambermd.org/tutorials/advanced/tutorial21/section1.htm 
		  /home/macenrola/Thesis/AMBER/explicit_water_binding/BENZENE-CB.rst7 -> /home/macenrola/Thesis/AMBER/explicit_water_binding/BENZENE-CB_iso.rst7
	"""
	with open(rstfile, 'rb') as r: # reads the file
		inlines = r.readlines()
	if maxdim < 0:
		maxdim = max([float(x) for x in inlines[-1].split()[:3]]) # locate the largest edge
	else: pass
	nline = '  {0:.7f}  {0:.7f}  {0:.7f}  {1:.7f}  {1:.7f}  {1:.7f}\n'.format(maxdim, 90) # corrects the last line
	inlines[-1] = nline # replaces it in the original file
	with open(rstfile[:-5]+suffix, 'wb') as w:
		w.writelines(inlines)


def get_prmtop_rst_script_solvation(bin_path = bin_path_serial, charge=0, watermol_to_remove=(0,0,0,0)):
	"""
	PRE: bin_path defines the path to the executable for tleap, antechamber, pmemd.MPI
		 charge specifies the total charge of the molecule, ideally 0 since the am1-bcc tends to spread the charge over the whole molecule which is not physical in most cases
		 watermol_to_remove indicates the number of watermolecules to remove to attain 1500 water molecule in all solvation model. The ORDER IS WATER, GUEST, CB7, COMPLEX 
	"""
	script = '\n'.join([
	'guest=${1:-test_guest}',
	'CB=${2:-CB}',
	'',
	'rm water.in',
	'echo "source $AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" >> water.in',
	'echo "structure = loadpdb water.pdb"  >> water.in',
	'echo "solvatebox structure TIP3PBOX 16.50 iso" >> water.in',
	'\n'.join(['echo "remove structure structure.{}" >> water.in'.format(x) for x in range(1501, 1501+watermol_to_remove[0])[::-1]]),
	'echo "saveamberparm structure water.prmtop water.rst7" >> water.in',
	'echo "quit" >> water.in',
	'',
	'{}tleap -s -f water.in > water-tleap.log'.format(bin_path),
	'# Default atom type is GAFF -at gaff',
	'{}antechamber -i $guest.pdb -fi pdb -o $guest.mol2 -fo mol2 -c bcc -nc {} -s 2 -j 4'.format(bin_path, charge),
	'# Unit names in the text are not variables but actual residue/unit names ><',
	'{}parmchk2 -i $guest.mol2 -f mol2 -o $guest.frcmod -f frcmod'.format(bin_path),
	'rm $guest.in',
	'echo "source leaprc.gaff" >> $guest.in',
	'echo "source $AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" >> $guest.in',
	'echo "loadAmberParams $guest.frcmod" >> $guest.in',
	'echo "GST = loadMol2 $guest.mol2" >> $guest.in',
	'echo "check GST" >> $guest.in',
	'echo "loadamberparams gaff.dat" >> $guest.in',
	'echo "solvatebox GST TIP3PBOX 14.5 iso" >> $guest.in',
	'\n'.join(['echo "remove GST GST.{}" >> $guest.in'.format(x) for x in range(1501, 1501+watermol_to_remove[1])[::-1]]),
	'echo "saveamberparm GST $guest.prmtop $guest.rst7" >> $guest.in',
	'echo "saveoff GST $guest.lib" >> $guest.in',
	'echo "quit" >> $guest.in',
	'{}tleap -s -f $guest.in > $guest-tleap.log'.format(bin_path),
	#'{}pmemd.MPI -O -i min.in -o $guest.out -p $guest.prmtop -c $guest.rst7 -r $guest.rst -x $guest.mdcrd -inf $guest.mdinfo'.format(bin_path),
	'',
	#'{}antechamber -i $CB.pdb -fi pdb -o $CB.mol2 -fo mol2 -s 2 -c gas -pf y -j 5 -nc 0'.format(bin_path), # Uses a CB.mol2 file that contain RESP charges, computed once and for all
	'{}parmchk2 -i $CB.mol2 -f mol2 -o $CB.frcmod -f frcmod'.format(bin_path),
	'rm $CB.in',
	'echo "source leaprc.gaff" >> $CB.in',
	'echo "source $AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" >> $CB.in',
	'echo "loadAmberParams $CB.frcmod" >> $CB.in',
	'echo "CB7 = loadMol2 $CB.mol2" >> $CB.in',
	'echo "check CB7" >> $CB.in',
	'# echo "source leaprc.gaff" >> $CB.in',
	'echo "loadamberparams gaff.dat" >> $CB.in',
	'echo "solvatebox CB7 TIP3PBOX 10.3 iso" >> $CB.in',
	'\n'.join(['echo "remove CB7 CB7.{}" >> $CB.in'.format(x) for x in range(1501, 1501+watermol_to_remove[2])[::-1]]),
	'echo "saveamberparm CB7 $CB.prmtop $CB.rst7" >> $CB.in',
	'echo "saveoff CB7 $CB.lib" >> $CB.in',
	'echo "quit" >> $CB.in',
	'{}tleap -s -f $CB.in > $CB-tleap.log'.format(bin_path_serial),
	#'{}pmemd.MPI -O -i min.in -o $CB.out -p $CB.prmtop -c $CB.rst7 -r pmemd.MPI_restart.rst -x $CB.mdcrd -inf $CB.mdinfo'.format(bin_path),
	'',
	'rm $guest-$CB.in',
	"echo '$guest-$CB.in'",
	'echo "source leaprc.gaff" >> $guest-$CB.in',
	'echo "source $AMBERHOME/dat/leap/cmd/oldff/leaprc.ff14SB" >> $guest-$CB.in',
	'echo "loadoff $guest.lib" >> $guest-$CB.in',
	'echo "loadoff $CB.lib" >> $guest-$CB.in',
	'echo "loadamberparams $CB.frcmod" >> $guest-$CB.in',
	'echo "loadamberparams $guest.frcmod" >> $guest-$CB.in',
	'echo "CB7 = loadmol2 $CB.mol2" >> $guest-$CB.in',
	'echo "GST = loadmol2 $guest.mol2" >> $guest-$CB.in',
	'echo "COMPLEX = loadPDB $guest-$CB.pdb" >> $guest-$CB.in',
	'echo "savemol2 COMPLEX $guest-$CB.mol2 1" >> $guest-$CB.in',
	'echo "solvatebox COMPLEX TIP3PBOX 10.3 iso" >> $guest-$CB.in',
	'\n'.join(['echo "remove COMPLEX COMPLEX.{}" >> $guest-$CB.in'.format(x) for x in range(1501, 1501+watermol_to_remove[3])[::-1]]),
	'echo "saveamberparm COMPLEX $guest-$CB.prmtop $guest-$CB.rst7" >> $guest-$CB.in',
	'echo "quit" >> $guest-$CB.in',
	'{}tleap -s -f $guest-$CB.in > $guest-$CB-tleap.log'.format(bin_path)])
	return script

def get_water_pdb():
	"""
	PRE: -
	POST: Return the PBD file for water as a string
	"""
	return 'HETATM    1  O   WAT A 1       -21.776  35.513  28.223  1.00 59.96           O'

def get_CB_mol2():
	"""
	PRE: -
	POST: Return the HF/6-31G* optimized geometry of CB7 along with its associated RESP charges in mol2 format supplied as text.
	"""
	# with open('/home/macenrola/Thesis/AMBER/explicit_water_binding/CB_RESP.mol2', 'rb') as r:
	# 	q = ''.join(r.readlines())
	# 	print q.__repr__()
	CBmol2 = '@<TRIPOS>MOLECULE\nCB7\n  126   147     1     0     0\nSMALL\nrc\n\n\n@<TRIPOS>ATOM\n      1 N1           1.1730    -4.4180     2.9510 n          2 CB7      -0.174369\n      2 C1          -0.0750    -4.6350     3.6690 c3         2 CB7      -0.092209\n      3 C2          -0.0760    -3.4910     4.7540 c3         2 CB7      -0.036585\n      4 N2           1.1720    -2.7870     4.5030 n          2 CB7      -0.197010\n      5 C3           1.9190    -3.3580     3.4710 c          2 CB7       0.676215\n      6 N3          -1.3120    -4.4060     2.9370 n          2 CB7      -0.176767\n      7 C4          -2.0520    -3.3350     3.4420 c          2 CB7       0.676759\n      8 N4          -1.3130    -2.7700     4.4840 n          2 CB7      -0.196670\n      9 O1           3.0690    -3.0470     3.1440 o          2 CB7      -0.542683\n     10 O2          -3.1940    -3.0100     3.0980 o          2 CB7      -0.542575\n     11 N5           1.1850    -0.4430     5.2690 n          2 CB7      -0.196745\n     12 C5          -0.0590    -0.0110     5.8880 c3         2 CB7      -0.010713\n     13 C6          -0.0440     1.5520     5.6860 c3         2 CB7      -0.101722\n     14 N6           1.2040     1.7890     4.9740 n          2 CB7      -0.170772\n     15 C7           1.9400     0.6200     4.7700 c          2 CB7       0.666667\n     16 N7          -1.2990    -0.4230     5.2450 n          2 CB7      -0.194436\n     17 C8          -2.0300     0.6540     4.7400 c          2 CB7       0.664340\n     18 N8          -1.2820     1.8110     4.9660 n          2 CB7      -0.166077\n     19 O3           3.0890     0.5500     4.3210 o          2 CB7      -0.540423\n     20 O4          -3.1740     0.6050     4.2740 o          2 CB7      -0.540387\n     21 N9           1.2530    -3.8140    -3.7080 n          2 CB7      -0.180228\n     22 C9           0.0120    -3.8820    -4.4670 c3         2 CB7      -0.083484\n     23 C10          0.0310    -2.5540    -5.3160 c3         2 CB7      -0.038882\n     24 N10          1.2780    -1.9150    -4.9190 n          2 CB7      -0.181216\n     25 C11          2.0090    -2.6780    -4.0050 c          2 CB7       0.674714\n     26 N11         -1.2330    -3.7870    -3.7190 n          2 CB7      -0.176292\n     27 C12         -1.9640    -2.6380    -4.0280 c          2 CB7       0.675310\n     28 N12         -1.2070    -1.8910    -4.9330 n          2 CB7      -0.184912\n     29 O5           3.1550    -2.4400    -3.6080 o          2 CB7      -0.543152\n     30 O6          -3.1100    -2.3780    -3.6470 o          2 CB7      -0.543337\n     31 N13          1.3100     4.4630    -2.8370 n          2 CB7      -0.173489\n     32 C13          0.0730     5.2300    -2.7780 c3         2 CB7      -0.092119\n     33 C14          0.0600     5.7880    -1.3030 c3         2 CB7      -0.052584\n     34 N14          1.2960     5.2700    -0.7360 n          2 CB7      -0.177124\n     35 C15          2.0420     4.5240    -1.6500 c          2 CB7       0.660386\n     36 N15         -1.1760     4.4890    -2.8690 n          2 CB7      -0.176672\n     37 C16         -1.9300     4.5500    -1.6950 c          2 CB7       0.662983\n     38 N16         -1.1900     5.2750    -0.7590 n          2 CB7      -0.178374\n     39 O7           3.1840     4.0860    -1.4750 o          2 CB7      -0.539695\n     40 O8          -3.0830     4.1300    -1.5450 o          2 CB7      -0.540146\n     41 N17          1.2380     3.8620     3.6410 n          2 CB7      -0.174910\n     42 C17         -0.0020     4.6240     3.6940 c3         2 CB7      -0.086939\n     43 C18          0.0180     5.4460     2.3490 c3         2 CB7      -0.019426\n     44 N18          1.2630     5.0330     1.7180 n          2 CB7      -0.190608\n     45 C19          1.9940     4.1350     2.4990 c          2 CB7       0.667825\n     46 N19         -1.2480     3.8750     3.6160 n          2 CB7      -0.168929\n     47 C20         -1.9790     4.1600     2.4610 c          2 CB7       0.666689\n     48 N20         -1.2220     5.0510     1.6960 n          2 CB7      -0.193126\n     49 O9           3.1390     3.7340     2.2690 o          2 CB7      -0.541045\n     50 O10         -3.1250     3.7750     2.2090 o          2 CB7      -0.540960\n     51 C21          1.8610    -0.7620    -5.5810 c3         2 CB7      -0.082124\n     52 C22         -1.8790    -1.7510     5.3500 c3         2 CB7      -0.041939\n     53 C23         -1.7600    -0.7300    -5.6060 c3         2 CB7      -0.078196\n     54 C24          1.7410    -1.7800     5.3810 c3         2 CB7      -0.039950\n     55 N21          1.1980    -5.2980     0.6480 n          2 CB7      -0.174350\n     56 C25         -0.0500    -5.9270     0.2400 c3         2 CB7      -0.082083\n     57 C26         -0.0330    -5.7850    -1.3300 c3         2 CB7      -0.048511\n     58 N22          1.2140    -5.0810    -1.5930 n          2 CB7      -0.187420\n     59 C27          1.9520    -4.8380    -0.4330 c          2 CB7       0.676758\n     60 N23         -1.2880    -5.2660     0.6260 n          2 CB7      -0.175000\n     61 C28         -2.0210    -4.8070    -0.4700 c          2 CB7       0.678034\n     62 N24         -1.2710    -5.0750    -1.6170 n          2 CB7      -0.188185\n     63 O11          3.1030    -4.3920    -0.3830 o          2 CB7      -0.543823\n     64 O12         -3.1640    -4.3380    -0.4390 o          2 CB7      -0.544116\n     65 C29          1.7460    -5.3460     1.9920 c3         2 CB7      -0.073817\n     66 C30         -1.8710    -5.3200     1.9560 c3         2 CB7      -0.071804\n     67 C31         -1.8160    -4.8710    -2.9470 c3         2 CB7      -0.068955\n     68 C32          1.7940    -4.9010    -2.9120 c3         2 CB7      -0.064415\n     69 C33          1.7860     3.1000     4.7480 c3         2 CB7      -0.067839\n     70 C34         -1.8330     3.1310     4.7180 c3         2 CB7      -0.079614\n     71 C35          1.8450     5.6670     0.5480 c3         2 CB7      -0.063406\n     72 C36         -1.7670     5.6860     0.5090 c3         2 CB7      -0.057792\n     73 N25         -1.1730     2.6990    -4.5640 n          2 CB7      -0.201096\n     74 C37          0.0770     2.5730    -5.3000 c3         2 CB7      -0.045586\n     75 C38          0.0670     1.0740    -5.7840 c3         2 CB7      -0.094638\n     76 N26         -1.1830     0.5540    -5.2510 n          2 CB7      -0.174629\n     77 C39         -1.9260     1.5230    -4.5720 c          2 CB7       0.683238\n     78 N27          1.3110     2.6730    -4.5340 n          2 CB7      -0.201172\n     79 C40          2.0440     1.4840    -4.5310 c          2 CB7       0.681940\n     80 N28          1.3000     0.5310    -5.2310 n          2 CB7      -0.172627\n     81 O13         -3.0790     1.3910    -4.1450 o          2 CB7      -0.544494\n     82 O14          3.1840     1.3310    -4.0790 o          2 CB7      -0.544389\n     83 C41         -1.7320     3.9560    -4.0990 c3         2 CB7      -0.033464\n     84 C42          1.8840     3.9190    -4.0560 c3         2 CB7      -0.033004\n     85 H1          -0.0820    -5.6420     4.1010 h2         2 CB7       0.140051\n     86 H2          -0.0870    -3.8700     5.7830 h2         2 CB7       0.128530\n     87 H3          -0.0700    -0.3130     6.9420 h2         2 CB7       0.119979\n     88 H4          -0.0390     2.1110     6.6300 h2         2 CB7       0.138668\n     89 H5           0.0050    -4.7890    -5.0820 h2         2 CB7       0.136914\n     90 H6           0.0350    -2.7310    -6.3980 h2         2 CB7       0.126956\n     91 H7           0.0900     6.0190    -3.5390 h2         2 CB7       0.138674\n     92 H8           0.0620     6.8840    -1.2540 h2         2 CB7       0.130786\n     93 H9          -0.0070     5.2580     4.5890 h2         2 CB7       0.135525\n     94 H10          0.0240     6.5310     2.5050 h2         2 CB7       0.121104\n     95 H11          1.7800    -0.9000    -6.6670 h2         2 CB7       0.108934\n     96 H12          2.9120    -0.7360    -5.2920 h2         2 CB7       0.128064\n     97 H13         -2.9300    -1.6560     5.0750 h2         2 CB7       0.119305\n     98 H14         -1.7980    -2.0880     6.3910 h2         2 CB7       0.096352\n     99 H15         -2.8160    -0.6830    -5.3370 h2         2 CB7       0.127376\n    100 H16         -1.6620    -0.8720    -6.6910 h2         2 CB7       0.108220\n    101 H17          2.7970    -1.7000     5.1220 h2         2 CB7       0.118675\n    102 H18          1.6390    -2.1200     6.4190 h2         2 CB7       0.096063\n    103 H19         -0.0660    -6.9690     0.5800 h2         2 CB7       0.136046\n    104 H20         -0.0300    -6.7490    -1.8530 h2         2 CB7       0.131588\n    105 H21          2.8040    -5.0950     1.9090 h2         2 CB7       0.125395\n    106 H22          1.6360    -6.3660     2.3810 h2         2 CB7       0.107144\n    107 H23         -2.9240    -5.0600     1.8470 h2         2 CB7       0.125043\n    108 H24         -1.7800    -6.3450     2.3370 h2         2 CB7       0.106860\n    109 H25         -2.8720    -4.6320    -2.8200 h2         2 CB7       0.125748\n    110 H26         -1.7130    -5.8030    -3.5170 h2         2 CB7       0.104799\n    111 H27          2.8520    -4.6810    -2.7610 h2         2 CB7       0.124944\n    112 H28          1.6870    -5.8370    -3.4750 h2         2 CB7       0.103666\n    113 H29          2.8410     2.9390     4.5250 h2         2 CB7       0.124366\n    114 H30          1.6890     3.6920     5.6670 h2         2 CB7       0.105736\n    115 H31         -2.8860     2.9910     4.4750 h2         2 CB7       0.127000\n    116 H32         -1.7410     3.7290     5.6340 h2         2 CB7       0.108416\n    117 H33          2.9000     5.3920     0.5370 h2         2 CB7       0.124552\n    118 H34          1.7470     6.7550     0.6520 h2         2 CB7       0.102118\n    119 H35         -1.6570     6.7730     0.6100 h2         2 CB7       0.100626\n    120 H36         -2.8250     5.4250     0.4720 h2         2 CB7       0.123211\n    121 H37          0.0950     3.2960    -6.1230 h2         2 CB7       0.132858\n    122 H38          0.0740     0.9700    -6.8760 h2         2 CB7       0.141147\n    123 H39         -2.7920     3.7780    -3.9140 h2         2 CB7       0.116000\n    124 H40         -1.6170     4.7060    -4.8930 h2         2 CB7       0.095775\n    125 H41          1.8030     4.6720    -4.8490 h2         2 CB7       0.095298\n    126 H42          2.9350     3.7200    -3.8440 h2         2 CB7       0.115858\n@<TRIPOS>BOND\n     1     1     2 1   \n     2     1     5 1   \n     3     1    65 1   \n     4     2     3 1   \n     5     2     6 1   \n     6     2    85 1   \n     7     3     4 1   \n     8     3     8 1   \n     9     3    86 1   \n    10     4     5 1   \n    11     4    54 1   \n    12     5     9 2   \n    13     6     7 1   \n    14     6    66 1   \n    15     7     8 1   \n    16     7    10 2   \n    17     8    52 1   \n    18    11    12 1   \n    19    11    15 1   \n    20    11    54 1   \n    21    12    13 1   \n    22    12    16 1   \n    23    12    87 1   \n    24    13    14 1   \n    25    13    18 1   \n    26    13    88 1   \n    27    14    15 1   \n    28    14    69 1   \n    29    15    19 2   \n    30    16    17 1   \n    31    16    52 1   \n    32    17    18 1   \n    33    17    20 2   \n    34    18    70 1   \n    35    21    22 1   \n    36    21    25 1   \n    37    21    68 1   \n    38    22    23 1   \n    39    22    26 1   \n    40    22    89 1   \n    41    23    24 1   \n    42    23    28 1   \n    43    23    90 1   \n    44    24    25 1   \n    45    24    51 1   \n    46    25    29 2   \n    47    26    27 1   \n    48    26    67 1   \n    49    27    28 1   \n    50    27    30 2   \n    51    28    53 1   \n    52    31    32 1   \n    53    31    35 1   \n    54    31    84 1   \n    55    32    33 1   \n    56    32    36 1   \n    57    32    91 1   \n    58    33    34 1   \n    59    33    38 1   \n    60    33    92 1   \n    61    34    35 1   \n    62    34    71 1   \n    63    35    39 2   \n    64    36    37 1   \n    65    36    83 1   \n    66    37    38 1   \n    67    37    40 2   \n    68    38    72 1   \n    69    41    42 1   \n    70    41    45 1   \n    71    41    69 1   \n    72    42    43 1   \n    73    42    46 1   \n    74    42    93 1   \n    75    43    44 1   \n    76    43    48 1   \n    77    43    94 1   \n    78    44    45 1   \n    79    44    71 1   \n    80    45    49 2   \n    81    46    47 1   \n    82    46    70 1   \n    83    47    48 1   \n    84    47    50 2   \n    85    48    72 1   \n    86    51    80 1   \n    87    51    95 1   \n    88    51    96 1   \n    89    52    97 1   \n    90    52    98 1   \n    91    53    76 1   \n    92    53    99 1   \n    93    53   100 1   \n    94    54   101 1   \n    95    54   102 1   \n    96    55    56 1   \n    97    55    59 1   \n    98    55    65 1   \n    99    56    57 1   \n   100    56    60 1   \n   101    56   103 1   \n   102    57    58 1   \n   103    57    62 1   \n   104    57   104 1   \n   105    58    59 1   \n   106    58    68 1   \n   107    59    63 2   \n   108    60    61 1   \n   109    60    66 1   \n   110    61    62 1   \n   111    61    64 2   \n   112    62    67 1   \n   113    65   105 1   \n   114    65   106 1   \n   115    66   107 1   \n   116    66   108 1   \n   117    67   109 1   \n   118    67   110 1   \n   119    68   111 1   \n   120    68   112 1   \n   121    69   113 1   \n   122    69   114 1   \n   123    70   115 1   \n   124    70   116 1   \n   125    71   117 1   \n   126    71   118 1   \n   127    72   119 1   \n   128    72   120 1   \n   129    73    74 1   \n   130    73    77 1   \n   131    73    83 1   \n   132    74    75 1   \n   133    74    78 1   \n   134    74   121 1   \n   135    75    76 1   \n   136    75    80 1   \n   137    75   122 1   \n   138    76    77 1   \n   139    77    81 2   \n   140    78    79 1   \n   141    78    84 1   \n   142    79    80 1   \n   143    79    82 2   \n   144    83   123 1   \n   145    83   124 1   \n   146    84   125 1   \n   147    84   126 1   \n@<TRIPOS>SUBSTRUCTURE\n     1 CB7         1 TEMP              0 ****  ****    0 ROOT\n'
	return CBmol2


def get_minimization_run_script(WaterName, GuestName, HostName, ComplexName, bin_path):
	"""
	PRE: -
	POST: Generates the script to run the minimizations
	"""
	script = []
	script.append('# Minimization with restraints of the different parts\n')
	script.append('{0}pmemd.MPI -O -i b2_min.in -p {1}.prmtop -c {1}_iso.rst7 -o {1}_min.out -r {1}_min.rst7 -ref {1}_iso.rst7'.format(bin_path, GuestName ))
	script.append('{0}pmemd.MPI -O -i cb7_min.in -p {1}.prmtop -c {1}_iso.rst7 -o {1}_min.out -r {1}_min.rst7 -ref {1}_iso.rst7'.format(bin_path, HostName))
	script.append('{0}pmemd.MPI -O -i cb7_b2_min.in -p {1}.prmtop -c {1}_iso.rst7 -o {1}_min.out -r {1}_min.rst7 -ref {1}_iso.rst7'.format(bin_path, ComplexName))

	script.append('# Minimization without restraints of the different parts\n')
	script.append('{0}pmemd.MPI -O -i all_min.in -p water.prmtop -c water_iso.rst7 -o water_min.out -r water_min.rst7'.format(bin_path))
	for name in [GuestName, HostName, ComplexName]:
		script.append('{0}pmemd.MPI -O -i all_min.in -p {1}.prmtop -c {1}_min.rst7 -o {1}_min.out -r {1}_min2.rst7'.format(bin_path, name))

	script.append('# Minimization with NVT\n')
	script.append('{0}pmemd.MPI -O -i nvt_md.in -c water_min.rst7 -p water.prmtop -o water_nvt.out -r water_nvt.rst7 -x water_nvt.nc'.format(bin_path))
	for name in [GuestName, HostName, ComplexName]:
		script.append('{0}pmemd.MPI -O -i nvt_md.in -c {1}_min2.rst7 -p {1}.prmtop -o {1}_nvt.out -r {1}_nvt.rst7 -x {1}_nvt.nc'.format(bin_path, name))

	script.append('# Minimization with MD1\n')
	for name in [WaterName, GuestName, HostName, ComplexName]:
		script.append('{0}pmemd.MPI -O -i npt_md.in -c {1}_nvt.rst7 -p {1}.prmtop -o {1}_npt.out -r {1}_npt.rst7 -x {1}_npt.nc'.format(bin_path, name))

	script.append('# Minimization with MD2\n')
	for name in [WaterName, GuestName, HostName, ComplexName]:
		script.append('{0}pmemd.MPI -O -i npt_md2.in -c {1}_npt.rst7 -p {1}.prmtop -o {1}_npt2.out -r {1}_npt2.rst7 -x {1}_npt2.nc'.format(bin_path, name))

	script.append('# Average the volumes\n')
	for name in [WaterName, GuestName, HostName, ComplexName]:
		script.append('grep VOL {}_npt2.out | head -n -2 | awk \'{{sum+=$9}}END{{printf \"%10.7f\\n\",(sum/NR)^(1/3)}}\''.format(name))

	return '\n'.join(script)




def get_minimization_route_script():
	"""
	PRE: The guest should be the residue 1 and cb residue 2
	POST: Generates the comvergence scripts named 
		  npt_md2.in
		  npt_md.in
		  nvt_md.in
		  all_min.in
		  cb7_b2_min.in
		  cb7_min.in
		  b2_min.in : min of water with restrain
	"""
	script = '\n'.join([
	'echo "Volume Equilibration MD - constant pressure" >> npt_md2.in',
	'echo " &cntrl" >> npt_md2.in',
	'echo "  imin  = 0," >> npt_md2.in',
	'echo "  irest    = 1," >> npt_md2.in',
	'echo "  ntx      = 5," >> npt_md2.in',
	'echo "  ntb      = 2," >> npt_md2.in',
	'echo "  pres0    = 1.0," >> npt_md2.in',
	'echo "  ntp      = 1," >> npt_md2.in',
	'echo "  taup     = 2.0," >> npt_md2.in',
	'echo "  cut      = 9.0," >> npt_md2.in',
	'echo "  ntr      = 0," >> npt_md2.in',
	'echo "  ntc      = 2," >> npt_md2.in',
	'echo "  ntf      = 2," >> npt_md2.in',
	'echo "  tempi    = 300.0," >> npt_md2.in',
	'echo "  temp0    = 300.0," >> npt_md2.in',
	'echo "  ntt      = 3," >> npt_md2.in',
	'echo "  gamma_ln = 5.0," >> npt_md2.in',
	'echo "  ig       = -1," >> npt_md2.in',
	'echo "  ioutfm   = 1," >> npt_md2.in',
	'echo "  iwrap    = 1," >> npt_md2.in',
	'echo "  nstlim   = 20000000," >> npt_md2.in',
	'echo "  dt       = 0.002," >> npt_md2.in',
	'echo "  ntpr = 50, ntwx = 10000, ntwr = 10000" >> npt_md2.in',
	'echo " /" >> npt_md2.in',
	'echo " &ewald" >> npt_md2.in',
	'echo "  nfft1=32," >> npt_md2.in',
	'echo "  nfft2=32," >> npt_md2.in',
	'echo "  nfft3=32," >> npt_md2.in',
	'echo "  order=4" >> npt_md2.in',
	'echo " /" >> npt_md2.in',

	'echo "Volume Equilibration MD - constant pressure" >> npt_md.in',
	'echo " &cntrl" >> npt_md.in',
	'echo "  imin  = 0," >> npt_md.in',
	'echo "  irest    = 1," >> npt_md.in',
	'echo "  ntx      = 5," >> npt_md.in',
	'echo "  ntb      = 2," >> npt_md.in',
	'echo "  pres0    = 1.0," >> npt_md.in',
	'echo "  ntp      = 1," >> npt_md.in',
	'echo "  taup     = 2.0," >> npt_md.in',
	'echo "  cut      = 10.0," >> npt_md.in',
	'echo "  ntr      = 0," >> npt_md.in',
	'echo "  ntc      = 2," >> npt_md.in',
	'echo "  ntf      = 2," >> npt_md.in',
	'echo "  tempi    = 300.0," >> npt_md.in',
	'echo "  temp0    = 300.0," >> npt_md.in',
	'echo "  ntt      = 3," >> npt_md.in',
	'echo "  gamma_ln = 5.0," >> npt_md.in',
	'echo "  ig       = -1," >> npt_md.in',
	'echo "  ioutfm   = 1," >> npt_md.in',
	'echo "  iwrap    = 1," >> npt_md.in',
	'echo "  nstlim   = 2500000," >> npt_md.in',
	'echo "  dt       = 0.002," >> npt_md.in',
	'echo "  ntpr = 1000, ntwx = 1000, ntwr = 1000" >> npt_md.in',
	'echo " /" >> npt_md.in',
	'echo " &ewald" >> npt_md.in',
	'echo "  nfft1=32," >> npt_md.in',
	'echo "  nfft2=32," >> npt_md.in',
	'echo "  nfft3=32," >> npt_md.in',
	'echo "  order=4" >> npt_md.in',
	'echo " /" >> npt_md.in',

	'echo "Heating Equilibration - constant volume" >> nvt_md.in',
	'echo " &cntrl" >> nvt_md.in',
	'echo "  imin   = 0," >> nvt_md.in',
	'echo "  irest  = 0," >> nvt_md.in',
	'echo "  ntx    = 1," >> nvt_md.in',
	'echo "  ntb    = 1," >> nvt_md.in',
	'echo "  cut    = 10.0," >> nvt_md.in',
	'echo "  ntr    = 0," >> nvt_md.in',
	'echo "  ntc    = 2," >> nvt_md.in',
	'echo "  ntf    = 2," >> nvt_md.in',
	'echo "  tempi  = 0.0," >> nvt_md.in',
	'echo "  temp0  = 300.0," >> nvt_md.in',
	'echo "  ntt    = 3," >> nvt_md.in',
	'echo "  ig     = -1," >> nvt_md.in',
	'echo "  ioutfm = 1," >> nvt_md.in',
	'echo "  iwrap  = 1," >> nvt_md.in',
	'echo "  gamma_ln = 1.0," >> nvt_md.in',
	'echo "  nstlim = 1000000, dt = 0.002" >> nvt_md.in',
	'echo "  ntpr = 1000, ntwx = 1000, ntwr = 1000" >> nvt_md.in',
	'echo " /" >> nvt_md.in',

	'echo "All atoms minimization" >> all_min.in',
	'echo " &cntrl" >> all_min.in',
	'echo "  imin   = 1," >> all_min.in',
	'echo "  maxcyc = 10000," >> all_min.in',
	'echo "  ncyc   = 500," >> all_min.in',
	'echo "  ntb    = 1," >> all_min.in',
	'echo "  ntr    = 0," >> all_min.in',
	'echo "  cut    = 10.0" >> all_min.in',
	'echo "/" >> all_min.in',

	'echo "CB7_B2: initial minimization solvent" >> cb7_b2_min.in',
	'echo " &cntrl" >> cb7_b2_min.in',
	'echo "  imin   = 1," >> cb7_b2_min.in',
	'echo "  maxcyc = 10000," >> cb7_b2_min.in',
	'echo "  ncyc   = 500," >> cb7_b2_min.in',
	'echo "  ntb    = 1," >> cb7_b2_min.in',
	'echo "  ntr    = 1," >> cb7_b2_min.in',
	'echo "  cut    = 10.0" >> cb7_b2_min.in',
	'echo "/" >> cb7_b2_min.in',
	'echo "Hold the Host (CB7) and Guest (GST) Fixed" >> cb7_b2_min.in',
	'echo "100.0" >> cb7_b2_min.in',
	'echo "RES 1 2" >> cb7_b2_min.in',
	'echo "END" >> cb7_b2_min.in',
	'echo "END" >> cb7_b2_min.in',

	'echo "CB7: initial minimization solvent" >> cb7_min.in',
	'echo " &cntrl" >> cb7_min.in',
	'echo "  imin   = 1," >> cb7_min.in',
	'echo "  maxcyc = 10000," >> cb7_min.in',
	'echo "  ncyc   = 500," >> cb7_min.in',
	'echo "  ntb    = 1," >> cb7_min.in',
	'echo "  ntr    = 1," >> cb7_min.in',
	'echo "  cut    = 10.0" >> cb7_min.in',
	'echo "/" >> cb7_min.in',
	'echo "Hold the Host (CB7) Fixed" >> cb7_min.in',
	'echo "100.0" >> cb7_min.in',
	'echo "RES 1" >> cb7_min.in',
	'echo "END" >> cb7_min.in',
	'echo "END" >> cb7_min.in',

	'echo "B2: initial minimization solvent" >> b2_min.in',
	'echo " &cntrl" >> b2_min.in',
	'echo "  imin   = 1," >> b2_min.in',
	'echo "  maxcyc = 10000," >> b2_min.in',
	'echo "  ncyc   = 500," >> b2_min.in',
	'echo "  ntb    = 1," >> b2_min.in',
	'echo "  ntr    = 1," >> b2_min.in',
	'echo "  cut    = 10.0" >> b2_min.in',
	'echo "/" >> b2_min.in',
	'echo "Hold the Guest (B2) Fixed" >> b2_min.in',
	'echo "100.0" >> b2_min.in',
	'echo "RES 1" >> b2_min.in',
	'echo "END" >> b2_min.in',
	'echo "END" >> b2_min.in'])

	# for fname in fnamelist:
	# 	with open(fname, 'rb') as r:
	# 		for line in r:
	# 			print 'echo "{}" >> {}'.format(line[:-1], fname.split('/')[-1][:-4]).__repr__() + ','
	return script

def prepare_the_system():
	"""
	PRE: The coordinates and topology files are available
	POST: Converge everything
	"""
	pass


def run_production_simu():
	"""
	PRE: The system is prepared
	POST: Runs the simulation production

	"""
	pass 

def compute_binding_energy():
	"""
	PRE: The production run is completed
	POST: Return the computed energy
	"""
	pass

def test_produce_tops(SMILES = 'C1#CC2C3C#CC4C1C1C#CC2C3C#CC41'):
	mol = Chem.MolFromSmiles(SMILES)
	mol = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol, AllChem.ETKDG())
	mol.SetProp('_Name', 'EXPANDED_CUBANE')
	print Chem.MolToMolBlock(mol)
	create_topology_file(Chem.MolToMolBlock(mol))



if __name__ == "__main__":
	print get_minimization_route_script()
	print get_minimization_run_script(WaterName, GuestName, HostName, ComplexName, bin_path)
	#test_produce_tops()
	# correct_rst_files_for_cuboid('/home/macenrola/Thesis/AMBER/explicit_water_binding/BENZENE-CB.rst7')
	# get_water_number_to_remove('/home/macenrola/Thesis/AMBER/explicit_water_binding/BENZENE-CB-tleap.log',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/water-tleap.log',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/BENZENE-tleap.log',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/CB-tleap.log',
	# 	)
	# get_minimization_script([
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/npt_md2.in.txt',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/npt_md.in.txt',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/nvt_md.in.txt',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/all_min.in.txt',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/cb7_b2_min.in.txt',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/cb7_min.in.txt',
	# 	'/home/macenrola/Thesis/AMBER/explicit_water_binding/b2_min.in.txt',
	# 	])