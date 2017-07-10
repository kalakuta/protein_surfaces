# assigns a point in x-y coordinates to a hexagonal grid cell
import math
import numpy as np

def assign(x,y,size):
	x = x / size
	y = y / size
	point = np.asarray([x,y])
	x1 = x - y / math.sqrt(3)
	y1 = 2 * y / math.sqrt(3)
	p = math.floor(x1)
	q = math.floor(y1)
	a = np.asarray([p + q/2 , q * math.sqrt(3)/2])
	b = np.asarray([p + q/2 + 1 , q * math.sqrt(3) / 2])
	c = np.asarray([p + q/2 + 1.5 , (q + 1) * math.sqrt(3) / 2])
	d = np.asarray([p + q/2 + 0.5 , (q + 1) * math.sqrt(3) / 2])
	min_dist = math.sqrt(3)
	#nearest_point = a
	for k in [a,b,c,d]:
		if np.linalg.norm(point - k) < min_dist:
			min_dist = np.linalg.norm(point - k)
			nearest_point = k
	p = p * size
	q = q * size
	if np.array_equal(nearest_point,a):
		return((p,q))
	if np.array_equal(nearest_point,b):
		return((p+size,q))
	if np.array_equal(nearest_point,c):
		return((p+size,q+size))
	if np.array_equal(nearest_point,d):
		return((p,q+size))	