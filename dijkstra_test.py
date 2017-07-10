import re
import sys
import math
import hydroscores as h
import hbondscores as hb
import assignhex as ah
import numpy as np
import scipy as cp
import scipy.sparse.csgraph

distances = np.empty((0,5))
row1 = [0, 1, 2, 1, 3]
distances=np.vstack([distances,row1])
row2 = [1, 0, 5, 0, 9]
distances=np.vstack([distances,row2])
row3 = [2, 5, 0, 7, 11]
distances=np.vstack([distances,row3])
row4 = [1, 0, 7, 0, 0]
distances=np.vstack([distances,row4])
row5 = [3, 9, 11, 0, 0]
distances=np.vstack([distances,row5])
sparse_distances = scipy.sparse.csr_matrix(distances)

geodesic_distances = scipy.sparse.csgraph.dijkstra(sparse_distances, indices=1)


print(sparse_distances)
print(geodesic_distances)

for i in range(5):
	print(geodesic_distances.item(i))