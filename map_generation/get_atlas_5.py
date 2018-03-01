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


def angle(vec1, vec2):
	vec1 = np.asarray(vec1)
	vec2 = np.asarray(vec2)
	mag1 = np.sqrt(np.dot(vec1, vec1))
	mag2 = np.sqrt(np.dot(vec2, vec2))
	try:
		return(np.arccos(max(min(1., (np.dot(vec1, vec2))/(mag1 * mag2)), -1.)))
	except FloatingPointError:
		print(vec1, vec2)
		return(0.0)


def least_squares_plane(points_list):
	sx, sy, sz = 0, 0, 0
	sxx, syy, szz = 0, 0, 0
	sxy, syz, szx = 0, 0, 0
	for point in points_list:
		x, y, z = [point[i] for i in range(3)]
		sx += x
		sy += y
		sz += z
		sxx += x * x
		syy += y * y
		szz += z * z
		sxy += x * y
		syz += y * z
		szx += z * x
	A = np.asarray([[sxx, sxy, szx], [sxy, syy, syz], [szx, syz, szz]])
	normal = np.dot(np.linalg.inv(A), [sx, sy, sz])
	d = 1 / np.sqrt(np.dot(normal, normal))
	normal = normal * d
	return[normal, d]

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
	# r = rotate_onto_z(direction[centre_point])

	# find points close to centre_point (within two angstroms on geodesic path)
	# find the average normal for these points using 'plane of best fit'
	# check direction is outward (negate if cosine is negative)
	
	neighbouring_points = [centre_point]
	for i in range(size):
		if geodesic_distances[centre_point][i] < 2:
			neighbouring_points.append(i)
	neighbour_list = []
	for point in neighbouring_points:
		neighbour_list.append(surface[point])
	
	#get normal of plane of best fit

	normal = least_squares_plane(neighbour_list)[0]
	
	#get rotation matrix onto this normal

	r = rotate_onto_z(normal)	
	
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

		for j in range(len(path)):
			path_points = np.vstack([path_points, close_points[path[j]]])

		#print(path_points)
		coord = np.asarray([0., 0.])
		if len(path) > 1:
			theta = 0.
			for l in range(len(path) - 1):
				g = geodesic_distances[close_points_index[path[l]], close_points_index[path[l+1]]]
				next = path_points[l + 1]
				yaw = np.angle(next[0] + next[1] * 1j)
				theta += yaw
				coord += np.asarray([g * np.cos(theta), g * np.sin(theta)])
				if next[0] == 0. and next[1] == 0.:
					pitch = np.pi / 2
				else:
					pitch = angle(next, [next[0], next[1], 0])
				if next[2] < 0:
					pitch = -pitch
				#translate origin onto 'next'
				path_points -= next
				#rotate x-axis onto projection of 'current-next' using yaw
				r1 = rotation_matrix([0, 0, 1], -yaw)
				#then rotate x-axis onto 'current-next' using pitch
				r2 = rotation_matrix([0, 1, 0], -pitch)
				R = np.dot(r2, r1)
				path_points = np.transpose(np.dot(R, np.transpose(path_points)))
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

