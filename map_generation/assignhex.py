# assigns a point in x-y coordinates to a hexagonal grid cell
import numpy as np

def assign(x,y,size):
	x = x / size
	y = y / size
	point = np.asarray([x,y])
	x1 = x - y / np.sqrt(3)
	y1 = 2 * y / np.sqrt(3)
	p = np.floor(x1)
	q = np.floor(y1)
	a = np.asarray([p + q/2 , q * np.sqrt(3)/2])
	b = np.asarray([p + q/2 + 1 , q * np.sqrt(3) / 2])
	c = np.asarray([p + q/2 + 1.5 , (q + 1) * np.sqrt(3) / 2])
	d = np.asarray([p + q/2 + 0.5 , (q + 1) * np.sqrt(3) / 2])
	min_dist = 10
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