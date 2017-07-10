#! /usr/bin/env python3
# test building geodesic distances
import sys
import numpy
import time
import scipy
import scipy.sparse.csgraph
import scipy.spatial.distance

start = (time.time())

filename = sys.argv[1]
txt = open(filename,"r")
surface = numpy.empty((0,3))
direction = {}
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
		surface = numpy.vstack([surface, [float(a4),float(a5),float(a6)]])
		direction[i] = [float(a9),float(a10),float(a11)]
		atom[i] = a3
		base[i] = a1
		seq[i] = a2
		density[i] = float(a8)
		i += 1
txt.close()

#generate matrix of distances between all pairs of surface points
distances = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(surface))

#replace distances over 1.5 with 0 [no direct connection in graph]
distances[distances > 1.5] = 0.

#create scipy sparse matrix with zeroes replaced by blanks
sparse_distances = scipy.sparse.csr_matrix(distances)

#create matrix of geodesic distances
geodesic_distances = scipy.sparse.csgraph.dijkstra(sparse_distances)

print('dijsktra done. Total run time = ' + str(time.time() - start) + ' seconds')


