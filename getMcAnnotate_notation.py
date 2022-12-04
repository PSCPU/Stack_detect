'''
Using this script we generate the mcAnnotate and FR3D version of annotation.

We need to pass 1 argument for running the script:
1) A file containing the list of stacks interactions.

Command:
> python mcAnnotate_stacks.py <PATH_TO_STACKS_FILE>
'''

import os
import sys
import copy
from utils import mcAnnotate_FR3D_notation,convert_modified_bases

file_name = sys.argv[1]
fo = open(file_name,'r')
a = fo.readlines()
a = [x for x in a if (('PDB Code' not in x) and ('num' not in x) and (x!='\n')) ]

final = []
for line in a:
	line_parts = line.split('\t')
	line_parts = list(map(str.strip,line_parts))
	nucs = [line_parts[2],line_parts[5]]
	basic_nucs = convert_modified_bases.convert_modified_bases(copy.deepcopy(nucs))
	faces = line_parts[9]
	mcAnnotate_notation,FR3D_notation = mcAnnotate_FR3D_notation.mcAnnotate_FR3D_notation(basic_nucs,faces)
	line_parts[9] = mcAnnotate_notation
	line_parts.append(FR3D_notation)
	s = ('\t').join(line_parts)
	final.append(s)
sys.stdout = open("final_output.txt", "w")
print ('pdb\tnum\tnuc\tchain\tnum\tnuc\tchain\tconsecutive/distant\tcis/trans\tMcAnnotate\t\ttopology\tdistance\ttheta\ttaui\ttauj\tsigma\tFR3D')
for line in final:
	print (line)
sys.stdout.close()

