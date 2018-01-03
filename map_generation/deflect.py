import numpy as np

def unit(vector):
	length = np.sqrt(np.dot(vector, vector))
	if length == 0:
		return vector
	else:
		return(vector / length)


def deflection(p1, n1, p2, n2, p3):
	n1 = unit(n1)
	n2 = unit(n2)
	r1 = unit(p2 - p1)
	r2 = unit(p3 - p2)
	angle = np.arccos(np.dot(unit(np.cross(n1, r1)),unit(np.cross(n2, r2))))
	sign = 0
	if np.dot(np.cross(n1, r1), r2) >= 0:
		sign = 1
	else:
		sign = -1
	return(angle * sign)


p1 = np.asarray([0, 0, 0])
n1 = np.asarray([0, 0, 1])
p2 = np.asarray([12, 4, -31])
print(np.cross(n1, unit(p2 - p1)))
n2 = np.asarray([0, 0, 1])
p3 = np.asarray([0, -5, -10])
print(unit(np.cross(n2, unit(p3 - p2))))
print(deflection(p1, n1, p2, n2, p3))
