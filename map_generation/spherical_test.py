import random
import numpy as np

def random_spherical():
    phi = random.random()*np.pi
    theta = random.random() * 2 * np.pi
    x, y, z = np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)
    return([x, y, z])
    
def spherical_deflection(p1, p2, p3):
    angle = np.arccos(np.dot(np.cross(p1, p2),np.cross(p2, p3)))
    sense = 0
    if np.dot(np.cross(p1, p2), p3) < 0:
        sense = 1
    elif np.dot(np.cross(p1, p2), p3) > 0:
        sense = -1
    return(angle * sense)


points = [[0., 0., -1.]]
for j in range(10):
	points.append(random_spherical())
path_points = np.empty((0, 3))
current_point = [0., 0., -1.]
while len(points) > 0:
    path_points = np.vstack([path_points, current_point])
    points.remove(current_point)
    min_length = 10.
    min_point = [10., 10., 10.]
    for point in points:
        vector = np.asarray(point) - np.asarray(current_point)
        if np.dot(vector,vector) < min_length:
            min_length = np.dot(vector, vector)
            min_point = point
    current_point = min_point

print(path_points)
        
    

'''

if len(path) == 1:
	geo_coords = np.vstack([geo_coords, [0., 0.]])
else:
	coord = np.asarray([0., 0.])
	for j in range(len(path) - 1):
		path_points = path_points - surface[path[j]]
		r = rotate_onto_z(direction[path[j]])
		path_directions = np.transpose(np.dot(r, np.transpose(path_directions)))
		g = geodesic_distances[path[j]][path[j+1]]
		euc = np.sqrt(np.dot(path_points[j+1], path_points[j+1]))
		coord += (g / euc) * np.asarray([path_points[j+1][0],path_points[j+1][1]])
	geo_coords = np.vstack([geo_coords, [coord]])
print(geo_coords)

'''

    
'''
p = [0., 0., -1.]
a = random_spherical()
b = random_spherical()


a = [0.5133411855906812, 0.33398446282342947, -0.79052843451004284]
b = [0.26131419870465467, 0.58525057562501437, -0.76759146248895549]

print(p)
print(a)
print(b)
print(spherical_deflection(p, a, b))


a = [1., 0., 0.]
b = [0.1, -0.1, .9899494936]


pa = np.asarray(a) - np.asarray(p)
pa_length = np.sqrt(np.dot(pa, pa))

ab = np.asarray(b) - np.asarray(a)
ab_length = np.sqrt(np.dot(ab, ab))


p_new = [0., 0.]
a_new = np.asarray([a[0], a[1]])
a_new = pa_length * a_new / np.sqrt(np.dot(a_new, a_new))

a_direction = np.arctan(a[1] / a[0])
if a[0] < 0:
    a_direction = a_direction + np.pi
    
angle = a_direction + spherical_deflection(p, a, b)
    
ab_new = np.asarray([ab_length * np.cos(angle), ab_length * np.sin(angle)])

b_new = a_new + ab_new

print(p_new)
print(a_new)
print(b_new)

'''

