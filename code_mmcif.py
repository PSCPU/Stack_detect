"""
R1 = SIGMA(ri*cos(i))/6 where 'ri' is the postion vector of each atom wrt center of the ring and N is the total number of atoms in the ring 
R2 = SIGMA(ri*sin(i))/6 where 'ri' is the postion vector of each atom wrt center of the ring and N is the total number of atoms in the ring 
N = R1 x R2 i.e the cross product of R1 and R2
"""
import sys
import math
from utils import dihedral,base_step
import re
import copy
from operator import sub,add
pi = math.pi

atom = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', "C1'"] #List of all atoms whose coordinates should be extracted
nuc = ['G', 'C', 'A', 'U', '1MA','A2M','M2G','H2U','7MG','2MU','OMG','5MC','OMC','PSU','4SU','5MU','OMU','2MG', '5AA','YYG'] #List of the 4 possible bases in an RNA molecule
atom5 = ['C4', 'C5', 'N7', 'C8', 'N9'] #List of the atoms in the 5 membered ring of a purine
atom6 = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'] #List of the atoms in the 6 membered ring of a purine or a pyrimidine
atom9 = ['N1', 'C2', 'N3', 'C4', 'N9', 'C8', 'N7', 'C5', 'C6']
nucAtoms = {} #Stores the nucleotides with all the atom and its (x,y,z) coordinates. key --> <nucleotide num>_<nucleotide name>_<chain name>
nucRings = {} #Stores the nucleotide rings(5,6) deatils - position vectors of each atoms, ring center, R1, R2 and normal(N). key --> <nucleotide num>_<nucleotide name>_<chain name>_<ring>
ringRingIntr = {} #Stores the ring-ring stacking interactions along with the parameters - consecutive/non-consecutive, cis/trans, face, distance(d), taui and tauj, theta
pu9Ring = {} #Stores the purine 9-membered ring deatils - position vectors of each atoms, ring center, R1, R2 and normal(N). key --> <nucleotide num>_<nucleotide name>_<chain name>_9
baseBaseStepParams = {} #Stores the base-base stacking interaction step parameters - shift, slide, rise, twist, roll, tilt
processed = [] #Stores the visited pairs of nucleotides
pur = ['A','G','1MA','A2M','M2G','G7M','OMG','2MG','5AA','YYG','7MG']
pyri = ['C','U','H2U','5MC','OMC','PSU','4SU','5MU','OMU','2MU'] 
#Digitize the numbers to 4 digits
def digitize(num):
	if len(num) == 1:
		num = '000' + num
	elif len(num) == 2:
		num = '00' + num
	elif len(num) == 3:
		num = '0' + num
	return num


#Cross product of two vectors; v1 and v2 are of type list
def cross(v1,v2):
	v1 = [float(x) for x in v1]
	v2 = [float(x) for x in v2]
	x = v1[1]*v2[2] - v1[2]*v2[1]
	y = v1[2]*v2[0] - v1[0]*v2[2]
	z = v1[0]*v2[1] - v1[1]*v2[0]
	return [format(x,'.4f'),format(y,'.4f'),format(z,'.4f')]


def cal_mod(v): #Calculate magnitude of a vector
	return math.sqrt((float(v[0])*float(v[0]))+(float(v[1])*float(v[1]))+(float(v[2])*float(v[2])))

def angle(v1,v2):	# Calculate angle between two vectors, (ax+cy+ez) and (bx+dy+fz)
	v1 = [float(x) for x in v1]
	v2 = [float(x) for x in v2]
	x=v1[0]*v2[0]
	y=v1[1]*v2[1]
	z=v1[2]*v2[2]
	modi=cal_mod(v1)
	modj=cal_mod(v2)
	ang = math.degrees(math.acos(((x+y+z)/modi)/modj))
	return format(ang,'.4f')

def normalize(v): #Normalize a vector v = v/|v|
	v = [float(x) for x in v]
	mod_v = cal_mod(v)
	x=v[0]/mod_v
	y=v[1]/mod_v
	z=v[2]/mod_v
	return [format(x,'.4f'),format(y,'.4f'),format(z,'.4f')]


def project(v,n): #Project a vector (v) on a plane with normal (n). Projection vector = v - (v.n)n
	v = [float(x) for x in v]
	n = [float(x) for x in n]
	d = v[0]*n[0]+v[1]*n[1]+v[2]*n[2]
	x = v[0] - d*n[0]
	y = v[1] - d*n[1]
	z = v[2] - d*n[2]
	return [format(x,'.4f'),format(y,'.4f'),format(z,'.4f')]

"""
-------------------------------------------------------------------------------------------------------------------------------
Begining of Part1 - Extracting Nucleotides from PDB
-------------------------------------------------------------------------------------------------------------------------------
"""
print ("Extracting PDB....")
file_name = sys.argv[1]
file_extn = file_name.split('.')[-1].lower()
fo = open(file_name,'r+')
lines = fo.readlines()

if (file_extn=="pdb"):
	#Loop to extract the coordinates of the atoms in pdb file. Coordinates of only the atoms present in the "atom" list are extracted.
	for l in lines:
		if l.strip()!="END" and l!="" and l.strip()!="ENDMOL": #Termination condition for reading the pdb file since there are files with same molecule overlapped with different chain name.
			l = l.strip("\n")
			if (l.startswith('ATOM') or l.startswith('HETATM')) and l[12:16].strip() in atom and l[17:20].strip() in nuc:
				key = l[22:27].strip()+'_'+l[17:20].strip()+'_'+l[20:22].strip()
				if key not in nucAtoms:
					nucAtoms[key] = {}
				nucAtoms[key][l[12:16].strip()] = [l[30:38].strip(),l[38:46].strip(),l[46:54].strip()]
	fo.close()
elif (file_extn=="cif"):
	#Loop to extract the coordinates of the atoms in pdb file. Coordinates of only the atoms present in the "atom" list are extracted.
	for l in lines:
		if l.strip()!="END" and l!="" and l.strip()!="ENDMOL": #Termination condition for reading the pdb file since there are files with same molecule overlapped with different chain name.
			l = l.strip("\n")
			if l.startswith('ATOM') or l.startswith('HETATM'):
				t = l.split()
				if (t[-2].strip().strip('"') in atom and t[-4].strip() in  nuc):
					key = t[-5].strip() + '_' + t[-4].strip() + '_' + t[-3].strip()
					if key not in nucAtoms:
						nucAtoms[key] = {}
					nucAtoms[key][t[-2].strip().strip('"')] = [t[10].strip(), t[11].strip(), t[12].strip()]
	fo.close()	

#If incomplete nucleotides are present in the pdb/cif file then delete them.
for key in nucAtoms.keys():
	nuc = key.split('_')[1]
	if nuc in pur:
		if len(nucAtoms[key])!=10:
			del nucAtoms[key]
	elif nuc in pyri:
		if len(nucAtoms[key])!=7:
			del nucAtoms[key]

#Loop to calucate the center of the ring, R1, R2 and normal vector to define the mean plane.
for key in nucAtoms.keys():
	x1=y1=z1=x2=y2=z2=x5=y5=z5=x6=y6=z6=x9=y9=z9=0.0
	nuc = key.split('_')[1]
	if nuc in pur:
		key6 = key + '_6' #Key for storing 6 membered ring:
		key5 = key + '_5' #Key for storing the 5 membered ring
		key9 = key + '_9' #Key for storing the 9 membered ring 
		nucRings[key6]={}
		nucRings[key5]={}
		pu9Ring[key9]={} 

		for atm in atom6: #Adding the x, y, z coordinates of the atoms in 6 membered ring of purines
			x6+=float(nucAtoms[key][atm][0])
			y6+=float(nucAtoms[key][atm][1])
			z6+=float(nucAtoms[key][atm][2])
		for atm in atom5: #Adding the x, y, z coordinates of the atoms in 5 membered ring of purines
			x5+=float(nucAtoms[key][atm][0])
			y5+=float(nucAtoms[key][atm][1])
			z5+=float(nucAtoms[key][atm][2])
		for atm in atom9: #Adding the x, y, z coordinates of the atoms in 9 membered ring of purines 
			x9+=float(nucAtoms[key][atm][0]) 
			y9+=float(nucAtoms[key][atm][1]) 
			z9+=float(nucAtoms[key][atm][2]) 
			
		nucRings[key6]['center'] = [format(x6/6,'.4f'),format(y6/6,'.4f'),format(z6/6,'.4f')] #Storing the center of 6 membered ring
		nucRings[key5]['center'] = [format(x5/5,'.4f'),format(y5/5,'.4f'),format(z5/5,'.4f')] #Storing the center of 5 membered ring
		pu9Ring[key9]['center'] = [format(x9/9,'.4f'),format(y9/9,'.4f'),format(z9/9,'.4f')] #Storing the center of 9 membered ring 
		nucRings[key6]['gly_center'] = nucRings[key5]['center'] 
		nucRings[key5]['gly_center'] = nucRings[key5]['center'] 
		pu9Ring[key9]['gly_center'] = nucRings[key5]['center'] 

		#Storing the position vectors of the atoms, center, r1, r2 and normal(n) of the 6-membered ring of purines
		for atm in atom6: 
			nucRings[key6][atm] = list(map(sub, [float(x) for x in nucAtoms[key][atm]], [float(x) for x in nucRings[key6]['center']])) #Storing the position vector of atoms with center as origin in 6 membered ring of purines
			x1=x1+(float(nucRings[key6][atm][0])*math.sin((2*pi*atom6.index(atm))/6))
			y1=y1+(float(nucRings[key6][atm][1])*math.sin((2*pi*atom6.index(atm))/6))
			z1=z1+(float(nucRings[key6][atm][2])*math.sin((2*pi*atom6.index(atm))/6))
			x2=x2+(float(nucRings[key6][atm][0])*math.cos((2*pi*atom6.index(atm))/6))
			y2=y2+(float(nucRings[key6][atm][1])*math.cos((2*pi*atom6.index(atm))/6))
			z2=z2+(float(nucRings[key6][atm][2])*math.cos((2*pi*atom6.index(atm))/6))
		nucRings[key6]['r1'] = normalize([format(x1,'.4f'),format(y1,'.4f'),format(z1,'.4f')])
		nucRings[key6]['r2'] = normalize([format(x2,'.4f'),format(y2,'.4f'),format(z2,'.4f')])
		nucRings[key6]['n'] = normalize(cross(nucRings[key6]['r1'],nucRings[key6]['r2']))
		x1=y1=z1=x2=y2=z2=0.0
		#Storing the position vectors of the atoms, center, r1, r2 and normal(n) of the 5-membered ring of purines
		for atm in atom5:
			nucRings[key5][atm] = list(map(sub, [float(x) for x in nucAtoms[key][atm]], [float(x) for x in nucRings[key5]['center']]))
			x1+=(float(nucRings[key5][atm][0])*math.sin((2*pi*atom5.index(atm))/5))
			y1+=(float(nucRings[key5][atm][1])*math.sin((2*pi*atom5.index(atm))/5))
			z1+=(float(nucRings[key5][atm][2])*math.sin((2*pi*atom5.index(atm))/5))
			x2+=(float(nucRings[key5][atm][0])*math.cos((2*pi*atom5.index(atm))/5))
			y2+=(float(nucRings[key5][atm][1])*math.cos((2*pi*atom5.index(atm))/5))
			z2+=(float(nucRings[key5][atm][2])*math.cos((2*pi*atom5.index(atm))/5))

		nucRings[key5]['r1'] = normalize([format(x1,'.4f'),format(y1,'.4f'),format(z1,'.4f')])
		nucRings[key5]['r2'] = normalize([format(x2,'.4f'),format(y2,'.4f'),format(z2,'.4f')])
		nucRings[key5]['n'] = normalize(cross(nucRings[key5]['r2'],nucRings[key5]['r1']))

		#Storing the position vectors of the atoms, center, r1, r2 and normal(n) of the 9-membered ring of purines
		for atm in atom9:
			pu9Ring[key9][atm] = map(sub,map(float,nucAtoms[key][atm]),map(float,pu9Ring[key9]['center'])) #Storing the position vector of atoms with center as origin in 9 membered ring of purines
			pu9Ring[key9][atm] = list(map(sub, [float(x) for x in nucAtoms[key][atm]], [float(x) for x in pu9Ring[key9]['center']]))

			x1+=(float(pu9Ring[key9][atm][0])*math.sin((2*pi*atom9.index(atm))/9))
			y1+=(float(pu9Ring[key9][atm][1])*math.sin((2*pi*atom9.index(atm))/9))
			z1+=(float(pu9Ring[key9][atm][2])*math.sin((2*pi*atom9.index(atm))/9))
			x2+=(float(pu9Ring[key9][atm][0])*math.cos((2*pi*atom9.index(atm))/9))
			y2+=(float(pu9Ring[key9][atm][1])*math.cos((2*pi*atom9.index(atm))/9))
			z2+=(float(pu9Ring[key9][atm][2])*math.cos((2*pi*atom9.index(atm))/9))

		pu9Ring[key9]['r1'] = normalize([format(x1,'.4f'),format(y1,'.4f'),format(z1,'.4f')])
		pu9Ring[key9]['r2'] = normalize([format(x2,'.4f'),format(y2,'.4f'),format(z2,'.4f')])
		pu9Ring[key9]['n'] = normalize(cross(pu9Ring[key9]['r2'],pu9Ring[key9]['r1']))
		

	# elif nuc=='U' or nuc=='C':
	elif nuc in pyri:
		key6 = key+'_6'
		nucRings[key6]={}
		for atm in atom6:
			x6+=float(nucAtoms[key][atm][0])
			y6+=float(nucAtoms[key][atm][1])
			z6+=float(nucAtoms[key][atm][2])
		nucRings[key6]['center'] = [format(x6/6,'.4f'),format(y6/6,'.4f'),format(z6/6,'.4f')] #Storing the center of pyrimidines
		nucRings[key6]['gly_center'] = nucRings[key6]['center']
		#Storing the position vectors of the atoms, center, r1, r2 and normal(n) of the 6-membered ring of pyrimidines
		for atm in atom6:
			nucRings[key6][atm] = map(sub,map(float,nucAtoms[key][atm]),map(float,nucRings[key6]['center'])) #Storing the position vector of atoms with center as origin of pyrimidines
			nucRings[key6][atm] = list(map(sub, [float(x) for x in nucAtoms[key][atm]], [float(x) for x in nucRings[key6]['center']]))

			x1+=(float(nucRings[key6][atm][0])*math.sin((2*pi*atom6.index(atm))/6))
			y1+=(float(nucRings[key6][atm][1])*math.sin((2*pi*atom6.index(atm))/6))
			z1+=(float(nucRings[key6][atm][2])*math.sin((2*pi*atom6.index(atm))/6))
			x2+=(float(nucRings[key6][atm][0])*math.cos((2*pi*atom6.index(atm))/6))
			y2+=(float(nucRings[key6][atm][1])*math.cos((2*pi*atom6.index(atm))/6))
			z2+=(float(nucRings[key6][atm][2])*math.cos((2*pi*atom6.index(atm))/6))
		nucRings[key6]['r1'] = normalize([format(x1,'.4f'),format(y1,'.4f'),format(z1,'.4f')])
		nucRings[key6]['r2'] = normalize([format(x2,'.4f'),format(y2,'.4f'),format(z2,'.4f')])
		nucRings[key6]['n'] = normalize(cross(nucRings[key6]['r1'],nucRings[key6]['r2']))

"""
-------------------------------------------------------------------------------------------------------------------------------
Begining of Part2 - Determining Stacks And Calculating Parameters
-------------------------------------------------------------------------------------------------------------------------------
"""
print ("Determining Stacks....")
def con_dist(key1,key2): #To check whether the two nucleotides are consecutive or distant
	key1 = key1.split('_')[0]
	key1 = float(re.split('[a-zA-z]+',key1)[0]) # Some nucleotide numbers are alphanumeric. Considering only the numeric part.
	key2 = key2.split('_')[0]
	key2 = float(re.split('[a-zA-z]+',key2)[0]) # Some nucleotide numbers are alphanumeric. Considering only the numeric part.
	if key1==key2:
		return "inter-chain"
	elif(abs(key1-key2)>1):
		return "non-consecutive"
	else:
		return "consecutive"

"""
If the C1'-center1-center2-C1' torsion angle is < 90 then its CIS else its TRANS.
"""
def cis_trans(key1,key2): #To check whether the two nucleotides are in CIS/TRANS orientation
	# center1 = map(float,nucRings[key1]['gly_center'])
	center1 = [float(x) for x in nucRings[key1]['gly_center']]
	# center2 = map(float,nucRings[key2]['gly_center'])
	center2 = [float(x) for x in nucRings[key2]['gly_center']]
	# c1_coord1 = map(float,nucAtoms[key1[:-2]]["C1'"])
	c1_coord1 = [float(x) for x in nucAtoms[key1[:-2]]["C1'"]]
	# c1_coord2 = map(float,nucAtoms[key2[:-2]]["C1'"])
	c1_coord2 = [float(x) for x in nucAtoms[key2[:-2]]["C1'"]]
	ang = dihedral.dihedral(c1_coord1,center1,center2,c1_coord2)
	if(ang>0.0 and ang<90.0):
		return ['CIS',format(ang,'.4f')]
	else:
		return ['TRANS',format(ang,'.4f')]

"""
Sit on 1st plane and watch 2nd plane for clockwise/counter-clockwise direction of numbering. 
If clockwise its the alpha face else its the beta face.
Silimarly sit on the 2nd plane and watch the 1st plane for alpha/beta face.
R1 x R2 runs in clockwise direction. Therefore R1 x R2 gives normal from the alpha face. 
If the center1->center2 vector makes an angle > 90 with the normal of ring1 then its beta else alpha.  

Alpha face is the one from which the normal comes out.
Beta face is the one from which the other side of alpha face.
"""
def face_orientation(key1,key2):
	dvec_c2_c1 = map(sub,map(float,nucRings[key1]['center']),map(float,nucRings[key2]['center'])) #vector direction from c1 to c2
	a = angle(nucRings[key1]['n'],dvec_c2_c1)
	if float(a)>90.0:
		f1='alpha'
	else:
		f1='beta'
	dvec_c1_c2 = map(sub,map(float,nucRings[key2]['center']),map(float,nucRings[key1]['center'])) #vector direction from c2 to c1
	a = angle(nucRings[key2]['n'],dvec_c1_c2)
	if float(a)>90.0:
		f2='alpha'
	else:
		f2='beta'
	return [f1,f2]	

def check_stack(key1,key2):
	if key1[:-2]!=key2[:-2]: #If both dont represent the same nucleotide.
		# dvec = map(sub,map(float,nucRings[key2]['center']),map(float,nucRings[key1]['center']))
		dvec = list(map(sub,[float(x) for x in nucRings[key2]['center']],[float(x) for x in nucRings[key1]['center']]))
		d = cal_mod(dvec)

		if((d<4.5) and (d>0)): # The distance between the two centers cannot be less than radius of H atom(0.5A). [488d has overlapped residues]
			theta = angle(nucRings[key1]['n'],nucRings[key2]['n']) #Theta is the angle between the two normals.
			theta = format(min(float(theta),180.0-float(theta)),'.4f') #theta must lie between 90 to -90
			if(float(theta)<23.0):
			#if (float(theta) <= 90 and float(theta) >= 60):
				taui = angle(nucRings[key1]['n'],dvec) #taui is the angle between distance vector and the normal of 1st ring.
				taui = format(min(float(taui),180.0-float(taui)),'.4f') #taui must lie between 90 to -90
				tauj = angle(nucRings[key2]['n'],dvec) #tauj is the angle between disance vector and the normal of 2nd ring.
				tauj = format(min(float(tauj),180.0-float(tauj)),'.4f') #tauj must lie between 90 to -90
				if(float(taui)<40.0 and float(tauj)<40.0):
					c_d = con_dist(key1,key2)
					c_t = cis_trans(key1,key2)
					f = face_orientation(key1,key2)
					return [True,format(d,'.4f'),theta,taui,tauj,c_d,c_t,f]
	return [False]

"""
y-axis for pyrimidine is in the direction of projection of vector joining center and N3 atom on the plane.
y-axis for 5/6 membered ring of purine is in the direction of projection of vector joining center and N1 atom on the corresponding plane.
"""
# def get_yaxis(key):
# 	temp_key = key[:-1] + '6'
# 	if(temp_key.split('_')[1]=='A' or temp_key.split('_')[1]=='G'):
# 		yaxis = normalize(project(nucRings[temp_key]['N1'],nucRings[temp_key]['n']))
# 	else:
# 		yaxis = normalize(project(nucRings[temp_key]['N3'],nucRings[temp_key]['n']))

# 	if(temp_key != key):
# 		yaxis = normalize(project(map(sub,map(float,nucAtoms[key[:-2]]['N1']),map(float,nucRings[key]['center'])),nucRings[key]['n']))
# 	return yaxis


def get_yaxis(key):
	# temp_key = key[:-1] + '6'
	# if(key.split('_')[1]=='A' or key.split('_')[1]=='G'):
	if (key.split('_')[1] in pur):
		yaxis = normalize(project(pu9Ring[key]['N1'],pu9Ring[key]['n']))
	else:
		yaxis = normalize(project(nucRings[key]['N3'],nucRings[key]['n']))

	# if(temp_key != key):
	# 	yaxis = normalize(project(map(sub,map(float,nucAtoms[key[:-2]]['N1']),map(float,nucRings[key]['center'])),nucRings[key]['n']))
	return yaxis
"""
x-axis is the cross product of normal(z-axis) and y-axis.
"""
def get_xaxis(key):
	# if(key.split('_')[1]=='A' or key.split('_')[1]=='G'):
	if (key.split('_')[1] in pur):
		yaxis = map(float,pu9Ring[key]['yaxis'])
		n = map(float,pu9Ring[key]['n'])
	else:
		yaxis = map(float,nucRings[key]['yaxis'])
		n = map(float,nucRings[key]['n'])
	xaxis = normalize(cross(yaxis,n))
	return xaxis


def get_step_parameters(key91, key92):
	if (sorted([key91[:-2], key92[:-2]]) in processed):
		return

	if(key91.split('_')[1] in pur):
		key91 = key91[:-1] + '9' 
	if(key92.split('_')[1] in pur):
		key92 = key92[:-1] + '9' 

	if key91 not in baseBaseStepParams:
		baseBaseStepParams[key91] = {}
	if key92 not in baseBaseStepParams:
		baseBaseStepParams[key92] = {}
	if key92 not in baseBaseStepParams[key91]:
		baseBaseStepParams[key91][key92] = {}
	if key91 not in baseBaseStepParams[key92]:
		baseBaseStepParams[key92][key91] = {}

	if ((key91.split('_')[1] in pur) and (key92.split('_')[1] in pur)):
		pu9Ring[key91]['yaxis'] = get_yaxis(key91)
		pu9Ring[key91]['xaxis'] = get_xaxis(key91)
		pu9Ring[key92]['yaxis'] = get_yaxis(key92)
		pu9Ring[key92]['xaxis'] = get_xaxis(key92)	

		t = base_step.base_step(copy.deepcopy(pu9Ring[key91]), copy.deepcopy(pu9Ring[key92]))
	elif ((key91.split('_')[1] in pur) and (key92.split('_')[1] not in pur)):
		pu9Ring[key91]['yaxis'] = get_yaxis(key91)
		pu9Ring[key91]['xaxis'] = get_xaxis(key91)
		nucRings[key92]['yaxis'] = get_yaxis(key92)
		nucRings[key92]['xaxis'] = get_xaxis(key92)	
		
		t = base_step.base_step(copy.deepcopy(pu9Ring[key91]), copy.deepcopy(nucRings[key92]))
	elif ((key91.split('_')[1] not in pur) and (key92.split('_')[1] in pur)):
		nucRings[key91]['yaxis'] = get_yaxis(key91)
		nucRings[key91]['xaxis'] = get_xaxis(key91)
		pu9Ring[key92]['yaxis'] = get_yaxis(key92)
		pu9Ring[key92]['xaxis'] = get_xaxis(key92)	
		
		t = base_step.base_step(copy.deepcopy(nucRings[key91]), copy.deepcopy(pu9Ring[key92]))
	elif ((key91.split('_')[1] not in pur) and (key92.split('_')[1] not in pur)):
		nucRings[key91]['yaxis'] = get_yaxis(key91)
		nucRings[key91]['xaxis'] = get_xaxis(key91)
		nucRings[key92]['yaxis'] = get_yaxis(key92)
		nucRings[key92]['xaxis'] = get_xaxis(key92)	
		
		t = base_step.base_step(copy.deepcopy(nucRings[key91]), copy.deepcopy(nucRings[key92]))

	for params in t.keys():
		baseBaseStepParams[key91][key92][params] = t[params]
		baseBaseStepParams[key92][key91][params] = t[params]

	processed.append(sorted([key91[:-2], key92[:-2]]))




# Loop to check if two rings are stacked. If they are stacked then calucate the center-center distance, theta, taui and tauj and base step parameters
key_list = list(nucRings.keys())
for i in range(len(key_list)):
	key1 = key_list[i]
	for j in range(i+1,len(key_list)):
		key2 = key_list[j]
		checker = check_stack(key1,key2) #Function to check whether two rings are stacking.
		if checker[0]==True:
			p={}
			p['distance'] = checker[1]
			p['theta'] = checker[2]
			p['taui'] = checker[3]
			p['tauj'] = checker[4]
			p['c_d'] = checker[5]
			p['c_t'] = checker[6]
			p['f'] = checker[7]
			# print key1 + '\t' + key2 + '\t' + format(checker[1],'.4f') + '\t' + checker[2] + '\t' + checker[3] + '\t' + checker[4] + '\t' + c_d + '\t' + c_t[0] + '\t' + c_t[1] +'\t' + f[0] + ' - ' + f[1] 
			if key1 not in ringRingIntr:
				ringRingIntr[key1] = {}
			if key2 not in ringRingIntr:
				ringRingIntr[key2] = {}
			if key2 not in ringRingIntr[key1]:
				ringRingIntr[key1][key2] = {}
			if key1 not in ringRingIntr[key2]:
				ringRingIntr[key2][key1] = {}


			# t = base_step.base_step(copy.deepcopy(nucRings[key1]),copy.deepcopy(nucRings[key2])) #Using deepcopy to send a copy of the dictionary, not the original dict.
			# p.update(t) # .update() adds the dictionary t to p. 
			for params in p.keys():
				ringRingIntr[key1][key2][params] = p[params] 
				ringRingIntr[key2][key1][params] = p[params] 
			ringRingIntr[key2][key1]['f'] = ringRingIntr[key2][key1]['f'][::-1]

			get_step_parameters(key1,key2)

"""
-------------------------------------------------------------------------------------------------------------------------------
PRINTING:
Print the rings along with the center(s), R1, R2 and normal.
Print the C1', N1 and N3 coordinates of each nucleotide
-------------------------------------------------------------------------------------------------------------------------------
"""
print ("Printing to File....")

def get_nomenclature(r1,r2):
	num1,nuc1,chain1,ring1 = r1.split('_')
	num2,nuc2,chain2,ring2 = r2.split('_')
	# if (nuc1>nuc2) or (nuc1==nuc2 and ring1>ring2) or(nuc1==nuc2 and ring1==ring2 and digitize(num1)>digitize(num2)) or (nuc1==nuc2 and ring1==ring2 and num1==num2 and chain1>chain2):
	if (digitize(num1)>digitize(num2)):
		r1,r2 = r2,r1
		num1,nuc1,chain1,ring1 = r1.split('_')
		num2,nuc2,chain2,ring2 = r2.split('_')
		# nuc1,nuc2 = nuc2,nuc1
		# ring1,ring2 = ring2,ring1
	try:
		nomenclature = 	nuc1+'('+ring1+')-'+ringRingIntr[r1][r2]['f'][0]+'||'+nuc2+'('+ring2+')-'+ringRingIntr[r1][r2]['f'][1]
	except:
		print (r1,r2)
		print (ringRingIntr[r1])
		print (ringRingIntr[r2])
		sys.exit(0)
	return nomenclature

repeat_list = [] # To check if a stacked pair has already been visited or not.

fo = open('ring_ring','w')
fo.write('------------------------------PDB Code: '+ file_name.split('/')[-1].split('.')[0] +'------------------------------\n')
fo.write('num\tnuc\tchain\tring\tnum\tnuc\tchain\tring\tconsecutive/non-consecutive\tC1\' orientation\tfaces\n')
for r1 in ringRingIntr.keys():
	for r2 in ringRingIntr[r1].keys():
		if sorted([r1,r2]) not in repeat_list:
			k1 = r1
			k2 = r2
			num1,nuc1,chain1,ring1 = k1.split('_')
			num2,nuc2,chain2,ring2 = k2.split('_')
			if (digitize(num1)>digitize(num2)):
				k1,k2 = k2,k1
				num1,nuc1,chain1,ring1 = k1.split('_')
				num2,nuc2,chain2,ring2 = k2.split('_')
			nomenclature = get_nomenclature(k1,k2)
			fo.write(num1+'\t'+nuc1+'\t'+chain1+'\t'+ring1+'\t'+num2+'\t'+nuc2+'\t'+chain2+'\t'+ring2+'\t'+ringRingIntr[k1][k2]['c_d']+'\t'+ringRingIntr[k1][k2]['c_t'][0]+'\t'+ringRingIntr[k1][k2]['f'][0]+'-'+ringRingIntr[k1][k2]['f'][1]+'\t'+nomenclature+'\n')
			repeat_list.append(sorted([k1,k2]))
fo.close()

repeat_list = []
fo = open('ring_ring_angles','w')
fo.write('------------------------------PDB Code: '+ file_name.split('/')[-1].split('.')[0] +'------------------------------\n')
fo.write('num\tnuc\tchain\tring\tnum\tnuc\tchain\tring\tdistance\ttheta\ttaui\ttauj\tsigma\n')
for r1 in list(ringRingIntr.keys()):
	for r2 in list(ringRingIntr[r1].keys()):
		if sorted([r1,r2]) not in repeat_list:
			k1 = r1
			k2 = r2
			num1,nuc1,chain1,ring1 = k1.split('_')
			num2,nuc2,chain2,ring2 = k2.split('_')
			if (digitize(num1)>digitize(num2)):
			# if (digitize(num1)>digitize(num2) or (digitize(num1)==digitize(num2) and nuc1 > nuc2) or (digitize(num1)==digitize(num2) and nuc1==nuc2 and chain1 > chain2)):
				k1,k2 = k2,k1
				num1,nuc1,chain1,ring1 = k1.split('_')
				num2,nuc2,chain2,ring2 = k2.split('_')
			fo.write(num1+'\t'+nuc1+'\t'+chain1+'\t'+ring1+'\t'+num2+'\t'+nuc2+'\t'+chain2+'\t'+ring2+'\t'+ringRingIntr[k1][k2]['distance']+'\t'+ringRingIntr[k1][k2]['theta']+'\t'+ringRingIntr[k1][k2]['taui']+'\t'+ringRingIntr[k1][k2]['tauj']+'\t'+ringRingIntr[k1][k2]['c_t'][1]+'\n')
			# fo.write(num1+'\t'+nuc1+'\t'+chain1+'\t'+ring1+'\t'+num2+'\t'+nuc2+'\t'+chain2+'\t'+ring2+'\t'+ringRingIntr[k1][k2]['distance']+'\t'+ringRingIntr[k1][k2]['theta']+'\t'+ringRingIntr[k1][k2]['taui']+'\t'+ringRingIntr[k1][k2]['tauj']+'\n')
			repeat_list.append(sorted([k1,k2]))
fo.close()

repeat_list = []

#fo = file('base_base_step_parameters','w')
#fo.write('------------------------------PDB Code: '+ file_name.split('/')[-1].split('.')[0] +'------------------------------\n')
#fo.write('num\tnuc\tchain\tnum\tnuc\tchain\tshift\tslide\trise\ttwist\troll\ttilt\n')
#for r1 in baseBaseStepParams.keys():
#	for r2 in baseBaseStepParams[r1].keys():
#		if sorted([r1,r2]) not in repeat_list:
#			k1 = r1
#			k2 = r2
#			# print k1, k2
#			num1,nuc1,chain1,ring1 = k1.split('_')
#			num2,nuc2,chain2,ring2 = k2.split('_')
#			if (digitize(num1)>digitize(num2)):
#			# if (digitize(num1)>digitize(num2) or (digitize(num1)==digitize(num2) and nuc1>nuc2) or (digitize(num1)==digitize(num2) and nuc1==nuc2 and chain1>chain2)):
#				k1,k2 = k2,k1
#				num1,nuc1,chain1,ring1 = k1.split('_')
#				num2,nuc2,chain2,ring2 = k2.split('_')
#			fo.write(num1+'\t'+nuc1+'\t'+chain1+'\t'+num2+'\t'+nuc2+'\t'+chain2+'\t'+baseBaseStepParams[k1][k2]['shift']+'\t'+baseBaseStepParams[k1][k2]['slide']+'\t'+baseBaseStepParams[k1][k2]['rise']+'\t'+baseBaseStepParams[k1][k2]['twist']+'\t'+baseBaseStepParams[k1][k2]['roll']+'\t'+baseBaseStepParams[k1][k2]['tilt']+'\n')
#			repeat_list.append(sorted([k1,k2]))
#o.close()#
print ("Successfully Comepleted....")
