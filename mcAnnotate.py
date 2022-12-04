import os
import sys
import copy
from utils import mcAnnotate_FR3D_notation,convert_modified_bases

file_name = sys.argv[1]
file_name = file_name.split('/')
FOLDER = '/'.join(file_name[:-1]) + '/'
pdb_list = [file_name[-1]]

final = []

for pdb in pdb_list:
	fo = open(FOLDER + '/' + pdb+'_base_base_alphaBeta','r')
	a = fo.readlines()
	a = [x for x in a if (('PDB' not in x) and ('num' not in x) and (x!='\n')) ]
	for line in a:
		line_parts = line.split('\t')
		line_parts = list(map(str.strip,line_parts))
		nucs = [line_parts[1].strip(),line_parts[4].strip()]
		basic_nucs = convert_modified_bases.convert_modified_bases(copy.deepcopy(nucs))
		faces = line_parts[8].strip()
		mcAnnotate_notation,FR3D_notation = mcAnnotate_FR3D_notation.mcAnnotate_FR3D_notation(basic_nucs,faces)
		line_parts[8] = mcAnnotate_notation
		line_parts.append(FR3D_notation)
		s = ('\t').join(line_parts)
		final.append(s)

print('------------------------------PDB Code: '+ pdb +'------------------------------')
print('num\tnuc\tchain\tnum\tnuc\tchain\tconsecutive/distant\tcis/trans\tMcAnnotate\t\ttopology\tdistance\ttheta\ttaui\ttauj\tsigma\tFR3D')
for line in final:
	print(line)