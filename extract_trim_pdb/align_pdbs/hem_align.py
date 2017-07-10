import numpy as np
import sys
import os

def rotation_matrix(axis, theta):
    
    # Use the Euler-Rodrigues formula to establish a rotation of theta radians about axis
    
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


	
infile = open(sys.argv[1],'r')
lines = infile.readlines()
infile.close()

# get coordinates of heme's nitrogens
for line in lines:
	if line[0:6] == 'HETATM' and line[17:20] == 'HEM':	
		if line[12:16].strip() == 'NA':
			NA = np.asarray([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
			#print('NA:', NA)
		if line[12:16].strip() == 'NB':
			NB = np.asarray([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
			#print('NB:', NB)
		if line[12:16].strip() == 'NC':
			NC = np.asarray([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
			#print('NA:', NA)
		if line[12:16].strip() == 'ND':
			ND = np.asarray([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
			#print('NA:', NA)

centroid = (NA + NB + NC + ND) / 4

# find 1st rotation axis and angle: CNA x CNB - to be parallel to z-axis
CNA = np.asarray(NA) - np.asarray(centroid)
CNB = np.asarray(NB) - np.asarray(centroid)
z = np.cross(CNA, CNB)
axis1 = np.cross(z, [0, 0, 1])
angle1 = np.arccos(np.dot(z, [0, 0, 1]) / np.sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]))

# find 2nd rotation: CNA to be parallel to x - axis
x = np.dot(rotation_matrix(axis1, angle1), CNA)
axis2 = [0, 0, 1]
angle2 = np.arccos(np.dot(x, [1, 0, 0]) / np.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]))
new_x = np.dot(rotation_matrix(axis2, angle2), x)

# check that new_x is parallel to x-axis
theta = np.arccos(np.dot(new_x, [1, 0, 0]) / np.sqrt(new_x[0] * new_x[0] + new_x[1] * new_x[1] + new_x[2] * new_x[2]))

if not np.abs(theta) < 0.001:
	angle2 = - angle2
	new_x = np.dot(rotation_matrix(axis2, angle2), x)
theta = np.arccos(np.dot(new_x, [1, 0, 0]) / np.sqrt(new_x[0] * new_x[0] + new_x[1] * new_x[1] + new_x[2] * new_x[2]))
print(theta)

# get single rotation matrix
R = np.dot(rotation_matrix(axis2, angle2), rotation_matrix(axis1, angle1))
#print(np.dot(R, CNA))


# write new .pdb file
pdb_code = os.path.basename(sys.argv[1][:-4])
outfile = pdb_code + '_trans.pdb'
f = open(outfile, 'w+')
for line in lines:
	if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM':
		coords = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
		coords = np.asarray(coords) - np.asarray(centroid)
		coords = np.dot(R, coords)
		entry = line[0:30]
		for coordinate in coords:
			entry += ('        ' + str(np.round(coordinate, 3)))[-8:]
		entry += line[54:]
		f.write(entry)
	else:
		f.write(line)
f.close()

