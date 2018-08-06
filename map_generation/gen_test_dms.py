import numpy as np

def get_random_point():
    theta = np.random.random()* 2 * np.pi
    phi = np.random.random() * np.pi
    return([30 * np.cos(theta) * np.sin(phi), 30 * np.sin(theta) * np.sin(phi), 30 * np.cos(phi)])

def get_rand_normal():
	phi = np.random.random() * 0.5
	theta = np.random.random() * 2 * np.pi
	return([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
	
def min_dist(p, p_list):
    min_dist = 1000.
    for point in p_list:
        vect = [p[i] - point[i] for i in range(3)]
        dist = np.dot(vect,vect)
        if dist < min_dist:
            min_dist = dist
    return(np.sqrt(min_dist))

'''
#points on a hexagonal lattice
points_list = [[x + .5 * y, y * np.sqrt(3) / 2] for x in np.arange(-49.5, 50.5, 1) for y in np.arange(-49.5, 50.5, 1)]

writefile = open('hex_lattice.dms', 'w+')
for i in range(len(points_list)):
	m = points_list[i][0]
	n = points_list[i][1]
	normal = get_rand_normal()
	writefile.write('GLY   1A    N  ' + str(m)[:5] + '\t' + str(n)[:5] + '\t0.0 SR0  ' + '0.9\t0.0\t0.0\t1.0\n')
writefile.close()


#points on hex lattice, normals randomly peturbed from [0, 0, 1]
points_list = [[x + .5 * y, y * np.sqrt(3) / 2] for x in np.arange(-49.5, 50.5, 1) for y in np.arange(-49.5, 50.5, 1)]
writefile = open('hex_lattice_variable_normal.dms', 'w+')
for i in range(len(points_list)):
	m = points_list[i][0]
	n = points_list[i][1]
	norm = get_rand_normal()
	writefile.write('GLY   1A    N  ' + str(m)[:5] + '\t' + str(n)[:5] + '\t0.0 SR0  ' + '0.9\t' + str(norm[0])[:5] + '\t' + str(norm[1])[:5] + '\t' +str(norm[2])[:5] + '\n')
writefile.close()
'''

#points on a hexagonal lattice with random peturbation
points_list = [[x + .5 * y + np.random.normal(0,0.18), np.random.normal(1,0.18) + y * np.sqrt(3) / 2] for x in np.arange(-49.5, 50.5, 0.7) for y in np.arange(-49.5, 50.5, 0.7)]

writefile = open('rand_hex_lattice.dms', 'w+')
for i in range(len(points_list)):
	m = points_list[i][0]
	n = points_list[i][1]
	writefile.write('GLY   1A    N  ' + str(m)[:5] + '\t' + str(n)[:5] + '\t0.0 SR0  ' + '0.9\t0.0\t0.0\t1.0\n')
writefile.close()
	
'''
# random points on a sphere
points_list = []
for i in range(10000):
    print(i)
    while True:
        point = get_random_point()
        if min_dist(point, points_list) > 0.7:
            points_list.append(point)
            break

writefile = open('rand_sphere.dms', 'w+')
for i in range(len(points_list)):
	x = points_list[i][0]
	y = points_list[i][1]
	z = points_list[i][2]
	writefile.write('GLY   1A    N  ' + str(x)[:5] + '\t' + str(y)[:5] + '\t' + str(z)[:5] + '\t0.9'+ str(x)[:5] + '\t' + str(y)[:5] + '\t' + str(z)[:5] + '\n')
writefile.close()
'''

'''	
# more uniform points on a sphere
points_list = []
points_list.append([0., 0., 1.])
for lat in range(1, 195, 1):
    print(len(points_list))
    phi = np.pi * lat / 194
    if lat % 2 == 1:
        for theta in np.arange(0, 2 * np.pi, 1 / (30 * np.sin(phi))):
            point = ([30 * np.cos(theta) * np.sin(phi), 30 * np.sin(theta) * np.sin(phi), 30 * np.cos(phi)])
            point = [np.round(value, 3) for value in point]
            if min_dist(point, points_list) > 0.9:
                points_list.append(point)
    elif lat % 2 == 0:
        start = 1 / (60 *np.sin(phi))
        for theta in np.arange(start, 2 * np.pi, 1 / (30 * np.sin(phi))):
            point = ([30 * np.cos(theta) * np.sin(phi), 30 * np.sin(theta) * np.sin(phi), 30 * np.cos(phi)])
            point = [np.round(value, 3) for value in point]
            if min_dist(point, points_list) > 0.9:
                points_list.append(point)

       
writefile = open('sphere_lattice.dms', 'w+')
for i in range(len(points_list)):
	x = (str(points_list[i][0]) + '000')[:5]
	y = (str(points_list[i][1]) + '000')[:5]
	z = (str(points_list[i][2]) + '000')[:5]
	nx = (str(points_list[i][0] / 30) + '000')[:5]
	ny = (str(points_list[i][1] / 30) + '000')[:5]
	nz = (str(points_list[i][2] / 30) + '000')[:5]
	writefile.write('GLY   1A    N  ' + x + '\t' + y + '\t' + z + '\t0.9\t' + nx + '\t' + ny + '\t' + nz + '\n')
writefile.close()
'''
