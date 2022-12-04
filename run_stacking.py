'''
Using this script we run the script to find the stacking interactions given a list of PDB/mmCIF files.

Requirements:
1) Linux machine
2) Python 2.7.x
3) Numpy

We need 1 argument for running the script:
1) Path to FOLDER containing the PDB/mmCIF files.

Command:
> python run_stacking.py <FOLDER>
'''
import os
import sys
import subprocess

FOLDER = sys.argv[1]

directory_list = os.listdir(FOLDER)
file_names = [x for x in directory_list if ('.pdb' in x) or ('.cif' in x)]

# p = '/usr/bin/python2.7 ' #Path of the python executable with latest version of numpy installed LOCAL.
# p = 'python ' #Python executable on the SERVER
p = 'python '
F = FOLDER
count = 1
for file_name in file_names:
	print(count),
	print(file_name.split('/')[-1])
	name = file_name.split('.')[0]
	FOLDER = F
	if not os.path.exists(FOLDER+name):
		os.system("mkdir " + FOLDER + name)
	os.system(p + " code_mmcif.py " + FOLDER + file_name)
	FOLDER = FOLDER + name + '/'
	os.system("move ring_ring " + FOLDER + name + "_ring_ring")
	os.system("move ring_ring_angles " + FOLDER + name + "_ring_ring_angles")
	os.system(p + " consolidate.py " + FOLDER + name + " > " + FOLDER + name + '_base_base')
	os.system(p + " sigma_distribution.py " + FOLDER + name + " > " + FOLDER + name + '_base_base_alphaBeta')
	os.system(p + " mcAnnotate.py " + FOLDER + name + " > " + FOLDER + name + '_base_base_mcAnnotate')
	print("---------------------------------")
	count+=1
