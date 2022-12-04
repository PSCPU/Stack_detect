'''
Using this script we consolidate all the stacking interactions into a single file with the following details:

pdb	num	nuc	chain	num	nuc	chain	consecutive/distant	cis/trans	face	topology	distance	theta	taui	tauj	sigma

Eg: 
157D	18	A	B	19	U	B	consecutive	cis	alpha-alpha	5||6, 6||6	3.8002	9.9413	15.47665	17.6212	22.5936
157D	7	U	A	8	U	A	consecutive	cis	alpha-beta	6||6	4.3531	11.9661	30.9323	36.7995	25.2653

We need 1 argument for running the above script:
1) Path to FOLDER containing the output of the stacking script for dataset.

Command:
> python consolidate_stacks.py > <OUTPUT_FILE>

'''
import os
import sys
import copy

FOLDER = sys.argv[1]
pdb_list = os.listdir(FOLDER)
pdb_list = [x for x in pdb_list if ('.pdb' not in x) and ('.cif' not in x) ]

final_store = {}


def normalize(num):
	if len(num) == 1:
		num = '000' + num
	elif len(num) == 2:
		num = '00' + num
	elif len(num) == 3:
		num = '0' + num
	return num

def generate_key(line,pdb): #Always starts with lower nucleotide number/name/chain separated by '_'
	num1,nuc1,chain1 = line[0],line[1],line[2]
	num2,nuc2,chain2 = line[3],line[4],line[5]
	if( (normalize(num2)<normalize(num1)) or ((normalize(num2)==normalize(num1)) and (nuc2<nuc1)) or ((normalize(num2)==normalize(num1)) and (nuc2==nuc1) and (chain2 < chain1))):
		num1,num2 = num2,num1
		nuc1,nuc2 = nuc2,nuc1
		chain1,chain2 = chain2,chain1
	key = num1 + '_' + nuc1 + '_' + chain1 + '_' + num2 + '_' + nuc2 + '_' + chain2 + '_' + pdb
	return key

def generate_key2(line,pdb): #Always starts with lower nucleotide number/name/chain separated by '_'
	num1,nuc1,chain1 = line[0],line[1],line[2]
	num2,nuc2,chain2 = line[4],line[5],line[6]
	if( (normalize(num2)<normalize(num1)) or ((normalize(num2)==normalize(num1)) and (nuc2<nuc1)) or ((normalize(num2)==normalize(num1)) and (nuc2==nuc1) and (chain2 < chain1))):
		num1,num2 = num2,num1
		nuc1,nuc2 = nuc2,nuc1
		chain1,chain2 = chain2,chain1
	key = num1 + '_' + nuc1 + '_' + chain1 + '_' + num2 + '_' + nuc2 + '_' + chain2 + '_' + pdb
	return key

def check_face(line,face):
	num1,nuc1,chain1 = line[0],line[1],line[2]
	num2,nuc2,chain2 = line[3],line[4],line[5]
	f1,f2 = face.split('-')
	if ((nuc1==nuc2) and f1>f2):
		f1,f2 = f2,f1
	return f1+'-'+f2

def get_parameters(b,pdb):
	store = {}
	for line in b:
		line_parts = line.split('\t')
		line_parts = list(map(str.strip,line_parts))
		key = generate_key2(line_parts[:8],pdb)
		if (key not in store):
			store[key] = {}
			store[key]['distance'] = [float(line_parts[8])]
			store[key]['theta'] = [float(line_parts[9])]
			store[key]['taui'] = [float(line_parts[10])]
			store[key]['tauj'] = [float(line_parts[11])]
			try:
				store[key]['sigma'] = [float(line_parts[12])]
			except:
				print(pdb)
				print(key)
				sys.exit(0)
		else:
			store[key]['distance'].append(float(line_parts[8]))
			store[key]['theta'].append(float(line_parts[9]))
			store[key]['taui'].append(float(line_parts[10]))
			store[key]['tauj'].append(float(line_parts[11]))
			store[key]['sigma'].append(float(line_parts[12]))
	return store

def calc_param(store,key,param):
	total = 0
	for p in store[key][param]:
		total += p
	total = total/len(store[key][param])
	return str(format(total,'.4f'))


for pdb in pdb_list:
	fo1 = open(FOLDER + pdb + '/' + pdb +'_base_base','r')
	a = fo1.readlines()
	a = [x for x in a if (('PDB Code' not in x) and ('num' not in x) and (x!='\n')) ]
	fo2 = open(FOLDER + pdb + '/' + pdb + '_ring_ring_angles','r')
	b = fo2.readlines()
	b = [x for x in b if (('PDB Code' not in x) and ('num' not in x) and (x!='\n'))]
	store_parameters = get_parameters(b,pdb)
	if (pdb not in final_store):
		final_store[pdb] = {}
	for line in a:
		line_parts = line.split('\t')
		line_parts = list(map(str.strip,line_parts))
		con_dist = line_parts[6].lower()
		ct = line_parts[7].lower()
		face = line_parts[8]
		topo = line_parts[-1]
		key = generate_key(line_parts[:6],pdb)
		face = check_face(line_parts[:6],face)
		if key not in final_store[pdb]:
			final_store[pdb][key] = {}
			final_store[pdb][key]['topology'] = topo
			final_store[pdb][key]['con_dist'] = con_dist
			final_store[pdb][key]['ct'] = ct
			final_store[pdb][key]['face'] = face
			final_store[pdb][key]['distance'] = calc_param(store_parameters,key,'distance')
			final_store[pdb][key]['theta'] = calc_param(store_parameters,key,'theta')
			final_store[pdb][key]['taui'] = calc_param(store_parameters,key,'taui')
			final_store[pdb][key]['tauj'] = calc_param(store_parameters,key,'tauj')
			final_store[pdb][key]['sigma'] = calc_param(store_parameters,key,'sigma')
fo1.close()
fo2.close()
sys.stdout = open("output.txt", "w")
print('pdb\tnum\tnuc\tchain\tnum\tnuc\tchain\tconsecutive/distant\tcis/trans\tface\ttopology\tdistance\ttheta\ttaui\ttauj\tsigma')
for pdb in pdb_list:
	for key in final_store[pdb]:
		stack = key.split('_')[:-1]
		stack = '\t'.join(stack)
		print(pdb + '\t' + stack + '\t' + final_store[pdb][key]['con_dist'] + '\t' + final_store[pdb][key]['ct'] + '\t' + final_store[pdb][key]['face'] + '\t' + final_store[pdb][key]['topology'] + '\t' + final_store[pdb][key]['distance'] + '\t' + final_store[pdb][key]['theta'] + '\t' + final_store[pdb][key]['taui'] + '\t' + final_store[pdb][key]['tauj'] + '\t' + final_store[pdb][key]['sigma'])
sys.stdout.close()
