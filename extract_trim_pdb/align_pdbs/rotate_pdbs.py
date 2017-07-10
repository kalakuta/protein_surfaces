'''
Read two pdb files with bound heme HETATM records;
translate each file to be centred at the centorid of
the heme's four nitrogen atoms and find the rmsd of the
optimal rotation aligning the heme molecules (based
on the four nitrogen atoms
Method and notation (matrix names) as per
Theobald - Acta Cryst (2005) A61 pp478-480

'''

import numpy as np
import sys
import os

def alignment_matrix(a, b):

	# takes two ordered sets of vectors a and b, and returns the rotation matrix
	# which best aligns them by rotating all the vectors in set b,
	# in the sense of minimising the sum of squared distances between their endpoints
	

	m = []
	p = []
	for i in range(a.shape[0]):
		xm = a[i][0] - b[i][0]
		xp = a[i][0] + b[i][0]
		ym = a[i][1] - b[i][1]
		yp = a[i][1] + b[i][1]
		zm = a[i][2] - b[i][2]
		zp = a[i][2] + b[i][2]
		m.append([xm,ym,zm])
		p.append([xp,yp,zp])

	e11 = 0
	e12 = 0
	e13 = 0
	e14 = 0
	e22 = 0
	e23 = 0
	e24 = 0
	e33 = 0
	e34 = 0
	e44 = 0

	for i in range(a.shape[0]):
		e11 += m[i][0] ** 2 + m[i][1] ** 2 + m[i][2] ** 2
		e12 += p[i][1] * m[i][2] - m[i][1] * p[i][2]
		e13 += m[i][0] * p[i][2] - p[i][0] * m[i][2]
		e14 += p[i][0] * m[i][1] - m[i][0] * p[i][1]
		e22 += p[i][1] ** 2 + p[i][2] ** 2 + m[i][0] ** 2
		e23 += m[i][0] * m[i][1] - p[i][0] * p[i][1]
		e24 += m[i][0] * m[i][2] - p[i][0] * p[i][2]
		e33 += p[i][0] ** 2 + p[i][2] ** 2 + m[i][1] ** 2
		e34 += m[i][1] * m[i][2] - p[i][1] * p[i][2]
		e44 += 	p[i][0] ** 2 + p[i][1] ** 2 + m[i][2] ** 2

	ematrix = np.matrix([[e11,e12,e13,e14],[e12,e22,e23,e24],[e13,e23,e33,e34],[e14,e24,e34,e44]])
	w , v = np.linalg.eig(ematrix)

	minw = w[0]
	min_index = 0
	for i in range(4):
		if w[i] < minw:
			minw = w[i]
			min_index = i

	evector = (v[:,min_index])
	theta = 2 * np.arccos(float(evector[0]))
	axis = [float(evector[1]),float(evector[2]),float(evector[3])]
	if theta != 0.0:
		R = rotation_matrix(axis,theta)
	else:
		R = np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]])
	return(R)
	#return(axis, theta)

def mult_q(q1, q2):

	# quaternion multiplication

	w1, x1, y1, z1 = q1
	w2, x2, y2, z2 = q2
	w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
	x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
	y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
	z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
	return(w, x, y, z)
	
	
def quaternion_rotate(v, q):

	v = np.array([0, v[0], v[1], v[2]])
	q_inv = q[0], -q[1], -q[2], -q[3]
	v_new = mult_q(mult_q(q, v), q_inv)[1:]
	
	return(v_new)

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
	

def rmsd_rotation(a, b):

	# find the minimum rmsd rotation for two sets of coordinates

	M = np.dot(b.transpose(),a)
	Sxx, Sxy, Sxz = M[0]
	Syx, Syy, Syz = M[1]
	Szx, Szy, Szz = M[2]
	key_matrix = [[Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
				  [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
				  [Szx - Sxz, Sxy + Syx, Syy - Szz - Sxx, Syz + Szy],
				  [Sxy - Syx, Szx + Sxz, Syz + Szy, Szz - Sxx - Syy]]
	rotation = np.linalg.eigh(key_matrix)[1][-1]
	return(rotation)

pdbA = open('1C2RA_heme_centred.pdb', 'r')
pdb_a = pdbA.readlines()
pdbA.close()

A = np.zeros(shape=(4,3))
for line in pdb_a:
	if line[0:6] == 'HETATM' and line[17:20] == 'HEM':
		if line[13:15] == 'NA':
			A[0] = [float(i) for i in line[32:56].split()]
		if line[13:15] == 'NB':
			A[1] = [float(i) for i in line[32:56].split()]
		if line[13:15] == 'NC':
			A[2] = [float(i) for i in line[32:56].split()]
		if line[13:15] == 'ND':
			A[3] = [float(i) for i in line[32:56].split()]	


pdb_list = os.listdir('translated_pdbs')
for pdbB in pdb_list:
	file = open('translated_pdbs/' + pdbB, 'r')
	pdb_b = file.readlines()
	file.close()
	
	# get Nx3 matrix of coordinates of heme nitrogens in each file
	
	B = np.zeros(shape=(4,3))

	for line in pdb_b:
		if line[0:6] == 'HETATM' and line[17:20] == 'HEM':
			if line[13:15] == 'NA':
				B[0] = [float(i) for i in line[32:56].split()]
			if line[13:15] == 'NB':
				B[1] = [float(i) for i in line[32:56].split()]
			if line[13:15] == 'NC':
				B[2] = [float(i) for i in line[32:56].split()]
			if line[13:15] == 'ND':
				B[3] = [float(i) for i in line[32:56].split()]	


	rotation_quaternion = rmsd_rotation(B, A)
	r1, r2, r3, r4 = rotation_quaternion
	
	theta = 2 * np.arccos(rotation_quaternion[0])
	axis = rotation_quaternion[1:]
	R = alignment_matrix(B, A)
	
	
	# write new .pdb file
	pdb_code = os.path.basename(pdbB[:-4])
	outfile = 'rotated_pdbs/' + pdb_code + '_rot.pdb'
	f = open(outfile, 'w+')
	for line in pdb_b:
		if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM':
			coords = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
			coords = np.dot(R, coords)
			entry = line[0:30]
			for coordinate in coords:
				entry += ('        ' + str(np.round(coordinate, 3)))[-8:]
			entry += line[54:]
			f.write(entry)
		else:
			f.write(line)
	f.close()

	

