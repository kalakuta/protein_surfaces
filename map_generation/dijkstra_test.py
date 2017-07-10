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
import time
from collections import defaultdict
import scipy
import scipy.sparse.csgraph

start = (time.time())

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
		surface = [i,float(a4),float(a5),float(a6)]
		direction[i] = [float(a9),float(a10),float(a11)]
		type[i] = a3
		base[i] = a1
		seq[i] = a2
		density[i] = float(a8)
		i += 1
txt.close()
print(surface)
'''
points_list = list(range(i))

centre_point = numpy.random.choice(points_list)
origin = numpy.asarray(surface[centre_point])

close_points = {}

for point in surface:
	close_points[point] = numpy.asarray(surface[point]) - origin

#print(surface[centre_point])
	

#create index for close_points
index = {}
n = 0
for point in close_points:
	index[n] = point
	n += 1

close_points_array = {}
for i in index:
	close_points_array[i] = numpy.array(close_points[index[i]])	
'''

distances = numpy.empty((0,len(surface)))

for point1 in surface.values():
	row=[]
	a, b, c = (point1[k] for k in range(3))
	for point2 in surface.values():
		p, q, r = (point2[k] for k in range(3))
		dist2 = (p - a) ** 2 + (q - b) ** 2	+ (r - c) ** 2
		row = numpy.append(row,dist2)
		
print(time.time()-start)

#create sparse matrix with zeroes replaced by blanks
'''
sparse_distances = scipy.sparse.csr_matrix(distances)

#create matrix of geodesic distances

geodesic_distances = scipy.sparse.csgraph.dijkstra(sparse_distances)

print('dijsktra done. Run time = ' + str(time.time() - start) + ' seconds')
'''

