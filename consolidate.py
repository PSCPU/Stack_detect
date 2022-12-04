import sys
import copy

file_name = sys.argv[1]
file_name = file_name.split('/')
FOLDER = '/'.join(file_name[:-1]) + '/'
pdb_list = [file_name[-1]]

store = {}

def normalize(num):
	if len(num) == 1:
		num = '000' + num
	elif len(num) == 2:
		num = '00' + num
	elif len(num) == 3:
		num = '0' + num
	return num


def generate_key(line): #Always starts with lower nucleotide number/name/chain separated by '_'
	num1,nuc1,chain1 = line[0],line[1],line[2]
	num2,nuc2,chain2 = line[4],line[5],line[6]
	if( (normalize(num2)<normalize(num1)) or ((normalize(num2)==normalize(num1)) and (nuc2<nuc1)) or ((normalize(num2)==normalize(num1)) and (nuc2==nuc1) and (chain2 < chain1))):
		num1,num2 = num2,num1
		nuc1,nuc2 = nuc2,nuc1
		chain1,chain2 = chain2,chain1
	key = num1 + '_' + nuc1 + '_' + chain1 + '_' + num2 + '_' + nuc2 + '_' + chain2
	return key

def get_face(line,face):
	num1,nuc1,chain1 = line[0],line[1],line[2]
	num2,nuc2,chain2 = line[4],line[5],line[6]
	face1,face2 = face.split('-')
	if( (normalize(num2)<normalize(num1)) or ((normalize(num2)==normalize(num1)) and (nuc2<nuc1)) or ((normalize(num2)==normalize(num1)) and (nuc2==nuc1) and (chain2 < chain1))):
		face1,face2 = face2,face1
	return face1+'-'+face2

def get_type(key): #pu_pu/pu_py/py_py
	key_parts = key.split('_')
	nuc1 = key_parts[1]
	nuc2 = key_parts[4]
	if ('A' in nuc1) or ('G' in nuc1):
		nuc1_type = 'pu'
	else:
		nuc1_type = 'py'

	if ('A' in nuc2) or ('G' in nuc2):
		nuc2_type = 'pu'
	else:
		nuc2_type = 'py'

	if nuc1_type != nuc2_type:
		return 'pu_py'
	return nuc1_type+'_'+nuc2_type

def generate_topology(key,rings):
	temp = {}
	temp['5'] = ''
	temp['6'] = ''
	key_parts = key.split('_')
	nuc1 = '_'.join(key_parts[:3])
	nuc2 = '_'.join(key_parts[3:])
	for ring in rings:
		ring_parts = ring.split('_')
		r1 = '_'.join(ring_parts[:4])
		r2 = '_'.join(ring_parts[4:])
		if(nuc1 in r1):
			temp[r1[-1]] += r2[-1]
		elif(nuc1  in r2):
			temp[r2[-1]] += r1[-1]
	temp['5'] = ''.join(sorted(temp['5']))
	temp['6'] = ''.join(sorted(temp['6']))
	if(temp['5'] != '' and temp['6'] != ''):
		return '5||'+temp['5'] + ', 6||' + temp['6']
	elif(temp['5'] != '' and temp['6'] == ''):
		return '5||'+temp['5']
	elif(temp['5']=='' and temp['6']!=''):
		return '6||'+temp['6']


for pdb in pdb_list:
	fo = open(FOLDER + pdb + '_ring_ring','r')
	a = fo.readlines()
	a = [x for x in a if (x!='\n' and 'num' not in x and 'PDB' not in x)]
	if pdb not in store:
		store[pdb] = {}
	for line in a:
		line_parts = line.split('\t')
		line_parts = list(map(str.strip, line_parts))
		key = generate_key(line_parts[:8])
		if key not in store[pdb]:
			store[pdb][key] = {}
			store[pdb][key]['con_dist'] = line_parts[8]
			store[pdb][key]['cis_trans'] = line_parts[9]
			store[pdb][key]['face'] = get_face(line_parts[:8],line_parts[10])
			store[pdb][key]['type'] = get_type(key)
			store[pdb][key]['rings'] = ['_'.join(line_parts[:8])]
		else:
			store[pdb][key]['rings'].append('_'.join(line_parts[:8]))


for pdb in store:
	print ('------------------------------PDB Code: '+ pdb +'------------------------------')
	print ('num\tnuc\tchain\tnum\tnuc\tchain\tconsecutive/non-consecutive\tC1\' orientation\tfaces\tinteracting bases\ttopology')
	for key in store[pdb]:
		store[pdb][key]['topology'] = generate_topology(key,store[pdb][key]['rings'])
		store[pdb][key]['consolidated'] = '\t'.join(key.split('_')) + '\t' + store[pdb][key]['con_dist'] + '\t' + store[pdb][key]['cis_trans'] + '\t' + store[pdb][key]['face'] + '\t' + store[pdb][key]['type'] + '\t' + store[pdb][key]['topology']
		print (store[pdb][key]['consolidated'])


