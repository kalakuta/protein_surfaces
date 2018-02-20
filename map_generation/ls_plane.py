
import numpy as np

def unit(normal):
	normal = np.asarray(normal)
	if np.dot(normal, normal) == 0:
		return(normal)
	else:
		return(normal / np.sqrt(np.dot(normal, normal)))

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

		
points = [[3, 5, 2], [8, 17, -1], [5, 5, 5], [4, -2, 9], [-1, -3, 17], [-2, -13, -23], [0, -1, -11]]
print(least_squares_plane(points))

# write a .kin file with points and least squares plane plotted

file = open('ls_plane.kin' , 'w+')
file.write('@text\n')
file.write('This file shows tests method for finding plane of best fit.' + '\n\n')
file.write('@kinemage\n\n')


file.write('@balllist {points} color=red radius=0.5\n')
for point in points:
	file.write(str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2]) + '\n')



#write all surface points to file
file.write('@dotlist {best fit plane} color=cyan\n')
a, b, c, d = least_squares_plane(points)[0][0], least_squares_plane(points)[0][1], least_squares_plane(points)[0][2], least_squares_plane(points)[1]
plane_points = [[x + .5 * y, y * np.sqrt(3) / 2] for x in np.arange(-49.5, 50.5, 1) for y in np.arange(-49.5, 50.5, 1)]
for point in plane_points:
	z = (d - a * point[0] - b * point[1]) / c
	file.write(str(point[0]) + ' ' + str(point[1]) + ' ' + str(z) + '\n')

file.close()


		
