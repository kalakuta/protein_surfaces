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
while len(points_list) > 0:
#for num in range(10):
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
	paths = {}
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
		paths[i] = path

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




	file = open('./runs/%s_%s/%s_%s.mp' % (domain, datestamp, domain, centre_point), 'w+')
	file.write(domain + '\t')
	file.write(str(centre_point) + '\t')
	file.write(seq[centre_point] + '\t')
	file.write(base[centre_point] + '\t')
	file.write(atom[centre_point] + '\t')
	file.write('\n')
	for i in range(k):
		point = close_points_index[i]
		file.write(str(point) + '\t')
		file.write(str(surface[point][0]) + '\t')
		file.write(str(surface[point][1]) + '\t')
		file.write(str(surface[point][2]) + '\t')
		file.write(str(geo_coords[i][0]) + '\t')
		file.write(str(geo_coords[i][1]) + '\t')
		file.write(base[point] + '\t' + seq[point] + '\t' + atom[point] + '\t')
		for l in range(len(paths[i])):
			file.write(str(paths[i][l]) + ' ')
		file.write('\n')

	file.close()

		
print(maps, time.time()-start)

