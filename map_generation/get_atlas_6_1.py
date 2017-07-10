#! /usr/bin/env python3
# generates maps covering (with overlap) a protein
# by selecting a random point; generating a map using a distance a centred here
# more blah
import re
import sys
import math
import hydroscores as h
import hbondscores as hb
import assignhex as ah
import numpy
from collections import defaultdict
import scipy
import scipy.sparse.csgraph
numpy.set_printoptions(threshold=numpy.nan)


def rotation_matrix(axis, theta):
    
    # Use the Euler-Rodrigues formula to establish a rotation of theta radians about axis
    
    axis = numpy.asarray(axis)
    theta = numpy.asarray(theta)
    axis = axis/math.sqrt(numpy.dot(axis, axis))
    a = math.cos(theta/2)
    b, c, d = -axis*math.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return numpy.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

#read .dms file to generate list of points and normals

filename = sys.argv[1]
domain = filename[:-4]		# strip '.dms' from filename to give pdb domain
txt = open(filename,"r")
surface = {}
direction = {}
type = {}
base = {}
seq = {}
density = {}
i = 0
while True:
	line = txt.readline()
	if not line: break
	if 'A\n' not in line:
		a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11 = line.split()
		surface[i] = [float(a4),float(a5),float(a6)]
		direction[i] = [float(a9),float(a10),float(a11)]
		type[i] = a3
		base[i] = a1
		seq[i] = a2
		density[i] = float(a8)
		i += 1
txt.close()

points_list = list(range(i))

while len(points_list) > 0:
	print(len(points_list))
	centre_point = numpy.random.choice(points_list)
	
	# reduce points_list set to points more than 8 angstroms from centre point
	for point in surface:
		a = numpy.array(surface[point])
		b = numpy.array(surface[centre_point])
		if numpy.linalg.norm(a-b) < 7:
			if point in points_list:
				points_list.remove(point)

	# reduce data set to points within 17 angstroms of centre point for centre cell
	# and 5 outer rings
	
	close_points = {}
	for point in surface:
		a = numpy.array(surface[point])
		b = numpy.array(surface[centre_point])
		if numpy.linalg.norm(a-b) < 17:
			close_points[point] = surface[point]

			
	# translate all points in close_points to have centre_point as origin
	origin = surface[centre_point]
	for point in close_points:
		close_points[point] = numpy.asarray(close_points[point]) - numpy.asarray(origin)

	#print(surface[centre_point])
		

	#create index for close_points
	index = {}
	n = 0
	for point in close_points:
		index[n] = point
		n += 1

		
	# populate a 2d dictionary (i.e graph) of distances within 17 angstroms of target point
	distances = numpy.empty((0,len(index)))
	for i in index:
		key1 = index[i]
		row=[]
		a = numpy.array(close_points[key1])
		for j in index:
			key2 = index[j]
			b = numpy.array(close_points[index[j]])
			dist = numpy.linalg.norm(a-b)		#calculates the euclidean distance between two points
			if dist < 1.5:						#so points within 1.5 angstroms are connected in the graph
				row=numpy.append(row,dist)
			else:
				row=numpy.append(row,0)			#on conversion to sparse matrix will represent no connection in graph
		distances=numpy.vstack([distances,row])

	#create sparse matrix with zeroes replaced by blanks

	sparse_distances = scipy.sparse.csr_matrix(distances)

	#create matrix of geodesic distances

	geodesic_distances = scipy.sparse.csgraph.dijkstra(sparse_distances)

	#find the index of the centre point
	centre_index = 0
	for i in index:
		if index[i] == centre_point:
			centre_index = i

	#extract only the points which are connected in the graph to the centre point and their geodesic distances (as geo_dist)
	#since a sphere centred at centre_point may well intersect surface at several disjoint sections

	connected_points = {}
	geo_dist = {}
	for i in index:
		if geodesic_distances.item((centre_index,i)) < 17:
			geo_dist[i] = geodesic_distances.item((centre_index,i))
			connected_points[i] = close_points[index[i]]

			
	# get average of unit vectors connecting atom_centre to point to use as the average normal for connected_points
	normal = numpy.asarray([0,0,0])
	for i in connected_points:
		normal = normal + numpy.divide(direction[index[i]] , len(connected_points))
	#print(normal)

	# rotate all points in connected_points so that normal is parallel to (0,0,1) - later, points will be projected on x-y plane
	# also rotate direction vector of each point as orientation

	# create rotation axis perpendicular to normal and to (0,0,1)
	axis = [normal[1],-normal[0],0]			# vector (a,-b,0) is perpendicular to (a,b,c) and to (0,0,1)
	# calculate required angle of rotation
	theta = math.acos(normal[2] / math.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]))
	# get rotation matrix
	r = rotation_matrix(axis, theta)
	# rotate each point with matrix r
	orientation = {}
	if normal[0] != 0 or normal[1] != 0:			#only use r if normal is not already parallel to (0,0,1)
		for point in connected_points:
			connected_points[point] = numpy.dot(r, numpy.asarray(connected_points[point]))
			orientation[point] = numpy.dot(r, numpy.asarray(direction[point]))
			
	print(orientation[centre_index])
				


	# project each point in reduced set onto x-y plane
	c_coords = {}
	for point in connected_points:
		c_coords[point] = [connected_points[point][0],connected_points[point][1]]
		
	# scale coordinates by the geodesic distance from centre_point
	geo_coords = {}

	for point in connected_points:
		a = c_coords[point][0]
		b = c_coords[point][1]
		geo_d = geo_dist[point]
		euc_dist = math.sqrt(a*a + b*b)
		if euc_dist != 0:
			sf = geo_d / euc_dist
		else:
			sf = 0
		geo_coords[point] = [a*sf,b*sf]



	# assign to a hexagonal grid - each point is in a particular hexagonal cell

	# size (= diameter of circumscribed circle) of each hexagonal cell

	cell_size = 3

	# assign each point in connected_points a hexagonal cell by establishing
	# which cell centre it's closest to

	hex_coords = {}
	for point in connected_points:
		hex_coords[point] = ah.assign(geo_coords[point][0],geo_coords[point][1],cell_size)

	#print(hex_coords)	



	# get hydrophobicity score for each point and weight according to sampling 'density' 
	# (dms file gives area associated with each point as 'density' so scale by this area)

	hydro_scores = {}
	for point in connected_points:
		lookup = base[index[point]] + type[index[point]]
		hydro_scores[point] = h.hscore(lookup) * density[index[point]]
		
	# assign each hexagonal cell a total hydrophobicity score

	cell_hydro_scores = {}
	for point in connected_points:
		if hex_coords[point] not in cell_hydro_scores:
			cell_hydro_scores[hex_coords[point]] = hydro_scores[point]
		else:
			cell_hydro_scores[hex_coords[point]] += hydro_scores[point]

	# get h-bonding scores for each point and weight according to sampling 'density' 
	# (dms file gives area associated with each point as 'density' so scale by this area)

	hbond_scores = {}
	for point in connected_points:
		lookup = base[index[point]] + type[index[point]]
		hbond_scores[point] = (hb.bondscore(lookup)[0] * density[index[point]] , hb.bondscore(lookup)[1] * density[index[point]])
		
	# assign each hexagonal cell a total h-bond donation score

	cell_donation_scores = {}
	for point in connected_points:
		if hex_coords[point] not in cell_donation_scores:
			cell_donation_scores[hex_coords[point]] = hbond_scores[point][0]
		else:
			cell_donation_scores[hex_coords[point]] += hbond_scores[point][0]

	# assign each hexagonal cell a total h-bond acceptance score

	cell_acceptance_scores = {}
	for point in connected_points:
		if hex_coords[point] not in cell_acceptance_scores:
			cell_acceptance_scores[hex_coords[point]] = hbond_scores[point][1]
		else:
			cell_acceptance_scores[hex_coords[point]] += hbond_scores[point][1]

			
			
	# standardise orientation for each point by converting to unit vector and scaling for density

	for point in connected_points:
		a = orientation[point]
		orientation[point] = density[point] * a / numpy.sqrt(numpy.dot(a,a))

	# assign each hexagonal cell an orientation unit vector (total of orientation for each point,
	# result converted to a unit vector)

	cell_orientations = {}
	for point in connected_points:
		if hex_coords[point] not in cell_orientations:
			cell_orientations[hex_coords[point]] = orientation[point]
		else:
			cell_orientations[hex_coords[point]] += orientation[point]

	# convert each cell's orientation vector to a unit vector.

	for point in cell_orientations:
		n = cell_orientations[point]
		cell_orientations[point] = n / numpy.sqrt(numpy.dot(n, n))




	# generate a list of cells with scores; write to .map file
	# creates a map with 8 concentric hexagonal rings

	s = int(cell_size)
	cell_list = [(x,y) for x in range(-5*s,6*s,s) for y in range(-5*s,6*s,s) if x+y < 6*s and x+y > -6*s]

	aa_lookup ={
	'A':'ALA',
	'R':'ARG',
	'N':'ASN',
	'D':'ASP',
	'C':'CYS',
	'E':'GLU',
	'Q':'GLN',
	'G':'GLY',
	'H':'HIS',
	'I':'ILE',
	'L':'LEU',
	'K':'LYS',
	'M':'MET',
	'F':'PHE',
	'P':'PRO',
	'S':'SER',
	'T':'THR',
	'W':'TRP',
	'Y':'TYR',
	'V':'VAL'
	}

	aa_rev_lookup = {
	'ALA':'A',
	'ARG':'R',
	'ASN':'N',
	'ASP':'D',
	'CYS':'C',
	'GLU':'E',
	'GLN':'Q',
	'GLY':'G',
	'HIS':'H',
	'ILE':'I',
	'LEU':'L',
	'LYS':'K',
	'MET':'M',
	'PHE':'F',
	'PRO':'P',
	'SER':'S',
	'THR':'T',
	'TRP':'W',
	'TYR':'Y',
	'VAL':'V'
	}


	# get amino acid mix for each cell
	cell_aa_makeup = {}
	for cell in cell_list:
		cell_aa_proportions = {}
		for aa in aa_lookup:
			cell_aa_proportions[aa] = 0.0
		count = 0
		for point in connected_points:
			if hex_coords[point] == cell:
				count += 1
		#print(count)
		if count != 0:
			for point in connected_points:
				if hex_coords[point] == cell:
					aa_type = base[point]
					aa_letter = aa_rev_lookup[aa_type]
					cell_aa_proportions[aa_letter] += 1.0 / count
		aa_entry = ' '
		for aa in cell_aa_proportions:
			if cell_aa_proportions[aa] != 0:
				aa_entry += aa + ' ' + str(cell_aa_proportions[aa]) + ', '
		cell_aa_makeup[cell] = aa_entry
				
			
			

	file = open('./runs/%s_%s.map' % (domain,centre_point), 'w+')
	file.write(domain + '\t')
	file.write(str(centre_point) + '\t')
	file.write(str(cell_size) + '\t')
	file.write(seq[centre_point] + '\t')
	file.write(base[centre_point] + '\t')
	file.write(type[centre_point] + '\n')

	for cell in cell_list:
		file.write(str(cell))
		file.write(' : ')
		count = 0
		for point in connected_points:
			if hex_coords[point] == cell:
				count += 1
		file.write(str(count) + ' : ')
		if count != 0:
			file.write(str(cell_hydro_scores[cell]) + ' : ')
			file.write(str(cell_acceptance_scores[cell]) + ' : ')
			file.write(str(cell_donation_scores[cell]) + ' : ')
			file.write(str(cell_orientations[cell][0]) + ', ' + str(cell_orientations[cell][1]) + ', ' + str(cell_orientations[cell][2]) + ' : ') 
			file.write(cell_aa_makeup[cell])
		else:
			file.write('0.0 : 0.0 : 0.0 : 0.0, 0.0, 1.0 : NULL')
		file.write('\n')

	file.close()
