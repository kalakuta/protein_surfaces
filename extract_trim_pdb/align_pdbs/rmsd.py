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


# find the minimum rmsd for two sets of coordinates

def rmsd(a, b):

	# takes two ordered sets of vectors a and b (as Nx3 matrices),
	# and returns the rotation matrix
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
	rmsd = np.linalg.eigh(ematrix)[0][0]
	
	return(rmsd)


pdb_rmsd = {}
pdb_list = os.listdir('cropped_pdbs')
for pdbA in pdb_list:
	file = open('cropped_pdbs/' + pdbA, 'r')
	pdb_a = file.readlines()
	file.close()
	pdb_rmsd[pdbA] = 0.0
	
	# get Nx3 matrix of coordinates of heme nitrogens in each file
	# translate them so their centroid is (0,  0, 0)

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

	centroid = (A[0] + A[1] + A[2] + A[3]) / 4		

	for i in range(4):
		A[i] -= centroid
	
	for pdbB in pdb_list:
		file = open('cropped_pdbs/' + pdbB, 'r')
		pdb_b = file.readlines()
		file.close() 
		
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

		centroid = (B[0] + B[1] + B[2] + B[3]) / 4		

		for i in range(4):
			B[i] -= centroid
		
		pdb_rmsd[pdbA] += rmsd(B, A)

results = open('results.txt', 'w+')	
results.write('pairwise total min rsmd values for 96 heme binding proteins\n')
results.write('with alignments by the heme nitrogen atoms\n\n')
		
min_rmsd = 100.0
best_pdb = ''
for pdb in pdb_rmsd:
	if pdb_rmsd[pdb] < min_rmsd:
		min_rmsd = pdb_rmsd[pdb]
		best_pdb = pdb
results.write('structure with best overall rmsd: ' + best_pdb + ' ' + str(min_rmsd / len(pdb_rmsd)) + '\n\n')
for pdb in pdb_rmsd:
	results.write(pdb + ': ' + str(pdb_rmsd[pdb]) + '\n')
results.close()