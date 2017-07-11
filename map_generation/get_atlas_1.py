#! /usr/bin/env python3

# program to generate an atlas covering an entire protein surface
# test building geodesic distances
import sys
import numpy as np
import time
import scipy
import scipy.sparse.csgraph
import scipy.spatial.distance

start = (time.time())


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

def step(x):
	if x == 0: return(1)
	else:
		return(x / abs(x))


filename = sys.argv[1]
txt = open(filename,"r")
surface = np.empty((0,3))
direction = np.empty((0,3))
atom = {}
base = {}
seq = {}
density = {}
i = 0
while True:
	line = txt.readline()
	if not line: break
	if 'A\n' not in line:
		a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11 = line.split()
		surface = np.vstack([surface, [float(a4),float(a5),float(a6)]])
		direction  = np.vstack([direction, [float(a9),float(a10),float(a11)]])
		atom[i] = a3
		base[i] = a1
		seq[i] = a2
		density[i] = float(a8)
		i += 1
txt.close()

points_list = list(range(i))
size = i


#generate matrix of distances between all pairs of surface points
distances = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(surface))

#replace distances over 1.5 with 0 [no direct connection in graph]
distances[distances > 1.5] = 0.

#create scipy sparse matrix with zeroes replaced by blanks
sparse_distances = scipy.sparse.csr_matrix(distances)

#create matrix of geodesic distances
geodesic_distances = scipy.sparse.csgraph.dijkstra(sparse_distances)

##	print('dijsktra done. Total run time = ' + str(time.time() - start) + ' seconds')

maps = 0
while len(points_list) > 0:
	maps += 1
	#print(len(points_list))
	#print(time.time()-start)
	centre_point = np.random.choice(points_list)
	origin = surface[centre_point]
	
	# reduce points_list set to points more than 7 angstroms from centre point
	for j in range(size):
		a = surface[j]
		b = surface[centre_point]
		if np.linalg.norm(a-b) < 8:
			if j in points_list:
				points_list.remove(j)
	
	# create map centred at centre_point
	
	# create array of points within 23 angstroms (geodesic) of centre_point
	
	close_points = np.empty((0,3))
	orientations = np.empty((0,3))
	close_points_index = {}
	k=0
	for j in range(size):
		if geodesic_distances[centre_point][j] < 23:
			close_points = np.vstack([close_points, surface[j]])
			orientations = np.vstack([orientations, direction[j]])
			close_points_index[k] = j
			k += 1
	
	
	centre_index = -1
	for k in range(len(close_points_index)):
		if close_points_index[k] == centre_point:
			centre_index = k
			
	#translate all points in close_points to have centre_point as origin
	for j in range(len(close_points_index)):
		close_points[j] = close_points[j] - origin

	#find average normal of points in close_points
	
	normal = np.asarray([0,0,0])
	for k in range(len(close_points_index)):
		normal = normal + np.divide(direction[close_points_index[k]] , len(close_points_index))
	
	# rotate all points in connected_points so that normal is parallel to (0,0,1)
	# later, points will be projected on x-y plane
	# also rotate direction vector of each point as orientation

	# create rotation axis perpendicular to normal and to (0,0,1)
	axis = [normal[1],-normal[0],0]			# vector (a,-b,0) is perpendicular to (a,b,c) and to (0,0,1)
	# calculate required angle of rotation
	theta = np.arccos(normal[2] / np.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]))
	# get rotation matrix
	r = rotation_matrix(axis, theta)
	
	#rotate close_points and orientations
	
	close_points = np.transpose(np.dot(r, np.transpose(close_points)))
	orientations = np.transpose(np.dot(r, np.transpose(orientations)))
	
	centre_index = -1
	for k in range(len(close_points_index)):
		if close_points_index[k] == centre_point:
			centre_index = k
			
	print(orientations[centre_index])
	
	
	

	# project each point in reduced set onto x-y plane
	# and scale by geodesic distance from centre_point
	geo_coords = np.empty((0, 2))
	for k in range(len(close_points_index)):	
		a, b = close_points[k][0], close_points[k][1]
		euc_dist = np.sqrt(a * a + b * b)
		geo_dist = geodesic_distances[centre_point][close_points_index[k]]
		if euc_dist != 0:
			sf = geo_dist / euc_dist
		else:
			sf = 0
		geo_coords = np.vstack([geo_coords, [a*sf,b*sf]])		
		


	# assign to a hexagonal grid - each point is in a particular hexagonal cell

	# size (= diameter of circumscribed circle) of each hexagonal cell

	cell_size = 3

	# assign each point in connected_points a hexagonal cell by establishing
	# which cell centre it's closest to
'''
	hex_coords = {}
	for point in connected_points:
		hex_coords[point] = ah.assign(geo_coords[point][0],geo_coords[point][1],cell_size)

	#print(hex_coords)	



	# get hydrophobicity score for each point and weight according to sampling 'density' 
	# (dms file gives area associated with each point as 'density' so scale by this area)

	hydro_scores = {}
	for point in connected_points:
		lookup = base[index[point]] + atom[index[point]]
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
		lookup = base[index[point]] + atom[index[point]]
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
	
	weighted_orientation = {}
	for point in connected_points:
		a = orientation[point]
		weighted_orientation[point] = density[point] * a / np.sqrt(np.dot(a,a))

	# assign each hexagonal cell an orientation unit vector (total of orientation for each point,
	# result converted to a unit vector)

	cell_orientations = {}
	for point in connected_points:
		if hex_coords[point] not in cell_orientations:
			cell_orientations[hex_coords[point]] = weighted_orientation[point]
		else:
			cell_orientations[hex_coords[point]] += weighted_orientation[point]

	# convert each cell's orientation vector to a unit vector.

	for point in cell_orientations:
		n = cell_orientations[point]
		cell_orientations[point] = n / np.sqrt(np.dot(n, n))

	# find a curvature score for each point in connected_points, using cell_orientations[(0,0)]
	# as the central normal and (0,0,0) as the centre point
	
	curvatures = {}
	n0 = cell_orientations[(0,0)]
	for point in connected_points:
		d = geo_dist[point]
		if d == 0:
			curvatures[point] = 0
		else:
			n1 = orientation[point]
			curvatures[point] = np.linalg.norm(n1-n0) * step(np.linalg.norm(n1 + connected_points[point] - n0) - d) / d
		
	# assign each hexagonal cell a total curvature score

	cell_curvature_scores = {}
	for point in connected_points:
		if hex_coords[point] not in cell_curvature_scores:
			cell_curvature_scores[hex_coords[point]] = curvatures[point]
		else:
			cell_curvature_scores[hex_coords[point]] += curvatures[point]

			


	# generate a list of cells with scores; write to .map file
	# creates a map with 6 concentric hexagonal rings including centre cell
	
	s = int(cell_size)
	cell_list = [(x,y) for x in range(-7*s,8*s,s) for y in range(-7*s,8*s,s) if x+y < 8*s and x+y > -8*s]

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
		cell_aa_count = {}
		for aa in aa_lookup:
			cell_aa_count[aa] = 0
		for point in connected_points:
			if hex_coords[point] == cell:
				aa_atom = base[point]
				aa_letter = aa_rev_lookup[aa_atom]
				cell_aa_count[aa_letter] += 1
		aa_entry = ' '
		for aa in cell_aa_count:
			if cell_aa_count[aa] != 0:
				aa_entry += aa + ' ' + str(cell_aa_count[aa]) + ', '
		cell_aa_makeup[cell] = aa_entry
				
			
			

	file = open('./runs/%s_%s.map' % (domain,centre_point), 'w+')
	file.write(domain + '\t')
	file.write(str(centre_point) + '\t')
	file.write(str(cell_size) + '\t')
	file.write(seq[centre_point] + '\t')
	file.write(base[centre_point] + '\t')
	file.write(atom[centre_point] + '\t')
	file.write(str(surface[centre_point][0]) + '\t' + str(surface[centre_point][1]) + '\t' + str(surface[centre_point][2]) + '\t')
	for entry in atom_list:
		file.write(entry + '\t')
	file.write('\n')
	

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
			file.write(str(cell_curvature_scores[cell]) + ' : ')
			file.write(cell_aa_makeup[cell])
		else:
			file.write('0.0 : 0.0 : 0.0 : 0.0, 0.0, 1.0 : 0.0 : NULL')
		file.write('\n')

	file.close()
'''
		
print(maps, time.time()-start)
	
	
	


