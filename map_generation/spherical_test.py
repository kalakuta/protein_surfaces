import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

def rotate_onto_z(vector):
	if vector[0] == 0 and vector[1] == 0:
		return np.identity(3)
	else:
		#gives the rotation matrix to rotate a given vector onto [0, 0, 1]
		axis = [vector[1],-vector[0],0]			# vector (a,-b,0) is perpendicular to (a,b,c) and to (0,0,1)
		theta = np.arccos(vector[2] / np.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]))
		return rotation_matrix(axis, theta)



points = [[0., 0., -1.]]
for j in range(10):
	points.append(random_spherical())
path_points = np.empty((0, 3))
path_points = np.vstack([path_points, [-1., 0., 0.]])
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

# get 2d coordinates
coords = np.empty((0, 2))
current = np.asarray([0., 0.])
angle = 0.
for i in range(10):
    coords = np.vstack([coords, current])
    vector = path_points[i + 2, :] - path_points[i + 1, :]
    length = np.sqrt(np.dot(vector, vector))
    angle += spherical_deflection(path_points[i], path_points[i + 1], path_points[i + 2])
    current += np.asarray([length * np.cos(angle), length * np.sin(angle)])
coords = np.vstack([coords, current])
print(coords)
    

x = coords[:, 0]
y = coords[:, 1]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_xlim(-5,5)
ax1.set_ylim(-5,5)
ax1.scatter(x, y)
ax1.plot(x, y)
#fig1.show()
fig1.savefig("2Dpath.png")

xs = path_points[1:, 0]
ys = path_points[1:, 1]
zs = path_points[1:, 2]

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(xs, ys, zs)
ax2.plot(xs, ys, zs)
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax2.plot_surface(x, y, z, color='#FFE6FF')

fig2.savefig("3Dpath.png")


path_points = path_points[1:,:]
normals = - path_points

coords2 = np.empty((0,2))
#current_point = [0., 0.]
coord = [0., 0.]
coords2 = np.vstack([coords2, coord])
for j in range(10):
	path_points = path_points - path_points[j]
	r = rotate_onto_z(normals[j])
	path_points = np.transpose(np.dot(r, np.transpose(path_points)))
	normals = np.transpose(np.dot(r, np.transpose(normals)))
	euc = np.sqrt(np.dot(path_points[j + 1], path_points[j + 1]))
	proj = np.sqrt(np.dot(path_points[j + 1][:2], path_points[j + 1][:2]))
	coord += (euc / proj) * np.asarray([path_points[j + 1][0],path_points[j + 1][1]])
	coords2 = np.vstack([coords2, coord])

print(coords2)

x = coords2[:, 0]
y = coords2[:, 1]


fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.set_xlim(-5,5)
ax3.set_ylim(-5,5)
ax3.scatter(x, y)
ax3.plot(x, y)
fig3.savefig("2Dpath_2.png")
    
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

