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
	
print(maps, time.time()-start)
	
	
	


