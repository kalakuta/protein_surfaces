#! /usr/bin/env python3

# program to generate an atlas covering an entire protein surface
# test building geodesic distances
import sys
import os
import numpy as np
np.seterr(all='raise')
import hydroscores as h
import hbondscores as hb
import assignhex as ah
import time
import scipy.sparse.csgraph
import scipy.spatial.distance

start = (time.time())


def rotation_matrix(axis, theta):
    
    # Use the Euler-Rodrigues formula to establish
	# a rotation of theta radians about axis
    
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

def rotate_onto_z(vector):
	if vector[0] == 0 and vector[1] == 0 and vector[2] >= 0:
		return np.identity(3)
	elif vector[0] == 0 and vector[1] == 0 and vector[2] < 0:
		return -np.identity(3)
	else:
		#gives the rotation matrix to rotate a given vector onto [0, 0, 1]
		axis = [vector[1],-vector[0],0]			# vector (a,-b,0) is perpendicular to (a,b,c) and to (0,0,1)
		theta = np.arccos(vector[2] / np.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]))
		return rotation_matrix(axis, theta)

def step(x):
	if x == 0: return(1)
	else:
		return(x / abs(x))


def unit(vector):
	length = np.sqrt(np.dot(vector, vector))
	if length == 0:
		return vector
	else:
		return(vector / length)


def deflection(p1, n1, p2, n2, p3):
	# calculates the angle of deflection when travelling the path
	# from point p1 -> p2 in the plane parallel to n1 and then the path
	# from p2 -> p3 in the plane parallel to n2
	n1 = unit(n1)
	n2 = unit(n2)
	r1 = unit(p2 - p1)
	r2 = unit(p3 - p2)
	try:
		angle = np.arccos(np.dot(unit(np.cross(n1, r1)),unit(np.cross(n2, r2))))
	except FloatingPointError:
		if np.dot(unit(np.cross(n1, r1)),unit(np.cross(n2, r2))) >= 1:
			angle = 0.
		elif np.dot(unit(np.cross(n1, r1)),unit(np.cross(n2, r2))) <= -1:
			angle = np.pi
	sign = 0
	if np.dot(np.cross(n1, r1), r2) >= 0:
		sign = 1
	else:
		sign = -1
	return(angle * sign)


datestamp = time.strftime("%d_%m_%y_%H:%M:%S")
filename = sys.argv[1]
domain = filename[:-4]
os.mkdir('runs/%s_%s' % (domain, datestamp))
txt = open(filename,"r")
surface = np.empty((0,3))
direction = np.empty((0,3))
atom = {}
base = {}
seq = {}
density = {}
i = 0

# read in surface records from .dms file

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
geodesic_distances, predecessors = scipy.sparse.csgraph.dijkstra(sparse_distances, return_predecessors=True)
print(time.time() - start)


maps = 0
#while len(points_list) > 0:
for num in range(10):
	maps += 1
	centre_point = np.random.choice(points_list)
	print(centre_point)
	origin = surface[centre_point]
	print(origin)
	
	# reduce points_list set to points more than 8 angstroms (geodesic) from centre point
	# increase this if possible?  In order to reduce size of atlas...
	for j in range(size):
		if geodesic_distances[centre_point, j] <= 8:
			if j in points_list:
				points_list.remove(j)
	
	# create map centred at centre_point
	
	# create array of the k points within 23 angstroms (geodesic) of centre_point
	
	close_points = np.empty((0,3))
	orientations = np.empty((0,3))
	close_points_index = {}
	reverse_index = {}
	k=0
	for j in range(size):
		if geodesic_distances[centre_point][j] < 23:
			close_points = np.vstack([close_points, surface[j]])
			orientations = np.vstack([orientations, direction[j]])
			close_points_index[k] = j
			reverse_index[j] = k
			k += 1

	# translate all points so centre_point has coordinates (0, 0, 0)

	close_points -= origin

	# rotate all points and normals so that centre_point's normal is
	# parallel to (0, 0, 1)
	
	r = rotate_onto_z(direction[centre_point])
	close_points = np.transpose(np.dot(r, np.transpose(close_points)))
	orientations = np.transpose(np.dot(r, np.transpose(orientations)))

	geo_coords = np.empty((0, 2))
	ratios = []
	for i in range(k):
		# get geodesic and euclidean distances from centre_point to point i
		geo = geodesic_distances[centre_point, close_points_index[i]]
		euc = scipy.spatial.distance.pdist(np.vstack([origin, surface[close_points_index[i]]]))
		# get geodesic path to point i
		path = []
		predecessor = close_points_index[i]
		while predecessor != centre_point:
			path = [predecessor] + path
			predecessor = predecessors[centre_point][predecessor]
		path = [centre_point] + path


		# index path by close_points
		path = [reverse_index[point] for point in path]


		path_points = np.empty((0, 3))
		path_directions = np.empty((0, 3))
		for j in range(len(path)):
			path_points = np.vstack([path_points, close_points[path[j]]])
			path_directions = np.vstack([path_directions, orientations[path[j]]])
		#print(path_points)
		coord = np.asarray([0., 0.])
		if len(path) > 1:
			theta = 0.
			previous = np.asarray([-1., 0., 0.])
			prev_norm = np.asarray([0., 0., 1.])
			for l in range(len(path) - 1):
				g = geodesic_distances[close_points_index[path[l]], close_points_index[path[l+1]]]
				current = path_points[l]
				current_norm = path_directions[l]
				next = path_points[l + 1]
				theta += deflection(previous, prev_norm, current, current_norm, next)
				coord += np.asarray([g * np.cos(theta), g* np.sin(theta)])
				previous = current
				previous_norm = current_norm
		geo_coords = np.vstack([geo_coords, [coord]])
		# get planar euclidean distance of new coordinates
		x = coord[0]
		y = coord[1]
		new_euc = np.sqrt(x * x + y * y)
		#if np.isnan(new_euc):
		#	print(previous, prev_norm, current, current_norm, next)
		#print(geo, euc[0], new_euc)
		if geo != 0:
			ratios.append(new_euc / geo)
		else:
			ratios.append(1.)
	#print(geo_coords)
	print(np.mean(ratios))







	#create a list of the residues and atoms present in map
	atom_list = []
	for k in range(len(close_points_index)):
		res_seq_atom = base[close_points_index[k]] + ':' + seq[close_points_index[k]] + ':' + atom[close_points_index[k]]
		if res_seq_atom not in atom_list:
			atom_list.append(res_seq_atom)


	# assign to a hexagonal grid - each point is in a particular hexagonal cell

	# size (= diameter of circumscribed circle) of each hexagonal cell

	cell_size = 3

	# assign each point in connected_points a hexagonal cell by establishing
	# which cell centre it's closest to

	hex_coords = {}
	for k in range(len(close_points_index)):
		hex_coords[k] = ah.assign(geo_coords[k][0], geo_coords[k][1], cell_size)


	# get hydrophobicity score for each point and weight according to sampling 'density'
	# get h-bonding scores for each point and weight according to sampling 'density'
	# (dms file gives area associated with each point as 'density' so scale by this area)

	hydro_scores = {}
	hbond_scores = {}
	for k in range(len(close_points_index)):
		j = close_points_index[k]
		dens = density[j]
		lookup = base[j] + atom[j]
		hydro_scores[k] = h.hscore(lookup) * dens
		hbond_scores[k] = (hb.bondscore(lookup)[0] * dens , hb.bondscore(lookup)[1] * dens)
		


	# assign each hexagonal cell a total hydrophobicity score
	# assign each hexagonal cell a total h-bond donation score
	# assign each hexagonal cell a total h-bond acceptance score
		
	cell_hydro_scores = {}
	cell_donation_scores = {}
	cell_acceptance_scores = {}
		
	for k in range(len(close_points_index)):
		cell_k = hex_coords[k]
		if cell_k not in cell_hydro_scores:
			cell_hydro_scores[cell_k] = hydro_scores[k]
		else:
			cell_hydro_scores[cell_k] += hydro_scores[k]

		if cell_k not in cell_donation_scores:
			cell_donation_scores[cell_k] = hbond_scores[k][0]
		else:
			cell_donation_scores[cell_k] += hbond_scores[k][0]

		if cell_k not in cell_acceptance_scores:
			cell_acceptance_scores[cell_k] = hbond_scores[k][1]
		else:
			cell_acceptance_scores[cell_k] += hbond_scores[k][1]
			
					
	# standardise orientation for each point by converting to unit vector and scaling for density

	weighted_orientations = {}
	for k in range(len(close_points_index)):
		a = orientations[k]
		weighted_orientations[k] = density[close_points_index[k]] * a / np.sqrt(np.dot(a,a))
	


	# assign each hexagonal cell an orientation unit vector (total of orientation for each point,
	# result converted to a unit vector)


	cell_orientations = {}
	for k in range(len(close_points_index)):
		cell_k = hex_coords[k]
		if cell_k not in cell_orientations:
			cell_orientations[cell_k] = weighted_orientations[k]
		else:
			cell_orientations[cell_k] += weighted_orientations[k]

	# convert each cell's orientation vector to a unit vector.

	for point in cell_orientations:
		n = cell_orientations[point]
		norm = np.sqrt(np.dot(n, n))
		if norm == 0:
			cell_orientations[point] = [0, 0, 1]
		else:
			cell_orientations[point] = n / norm

	# find a curvature score for each point in connected_points, using n0
	# as the central normal and (0,0,0) as the centre point
	# scale for density
	n0 = cell_orientations[(0, 0)]
	curvatures = {}
	for k in range(len(close_points_index)):
		point = close_points[k]
		d = np.linalg.norm(point)
		dens = density[close_points_index[k]]
		if d == 0:
			curvatures[k] = 0
		else:
			n1 = orientations[k]
			curvatures[k] = dens * np.linalg.norm(n1-n0) * step(np.linalg.norm(n1 + point - n0) - d) / d
	


	# assign each hexagonal cell a total curvature score

	cell_curvature_scores = {}
	for k in range(len(close_points_index)):
		cell_k = hex_coords[k]
		if cell_k not in cell_curvature_scores:
			cell_curvature_scores[cell_k] = curvatures[k]
		else:
			cell_curvature_scores[cell_k] += curvatures[k]

	


	# generate a list of cells with scores; write to .map file
	# creates a map with 8 concentric hexagonal rings including centre cell
	
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
		for k in range(len(close_points_index)):
			if hex_coords[k] == cell:
				aa_atom = base[close_points_index[k]]
				aa_letter = aa_rev_lookup[aa_atom]
				cell_aa_count[aa_letter] += 1
		aa_entry = ' '
		for aa in cell_aa_count:
			if cell_aa_count[aa] != 0:
				aa_entry += aa + ' ' + str(cell_aa_count[aa]) + ', '
		cell_aa_makeup[cell] = aa_entry
				
			
			

	file = open('./runs/%s_%s/%s_%s.map' % (domain, datestamp, domain, centre_point), 'w+')
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
		for k in range(len(close_points_index)):
			if hex_coords[k] == cell:
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

		
print(maps, time.time()-start)

