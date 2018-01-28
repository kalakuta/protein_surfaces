import numpy as np


#generate 10000 random 'surface' points

x = np.random.uniform(-50, 50, 10000)
y = np.random.uniform(-50, 50, 10000)

writefile = open('plane_uniform.dms', 'w+')

for i in range(10000):
	writefile.write('GLY   1A    N  ' + str(x[i])[:5] + '\t' + str(y[i])[:5] + '\t0.0 SR0  ' + '0.9\t0.0\t0.0\t1.0\n')

writefile.close()

points_list = [[x + .5 * y, y * np.sqrt(3) / 2] for x in np.arange(-49.5, 50.5, 1) for y in np.arange(-49.5, 50.5, 1)]

writefile = open('hex_lattice.dms', 'w+')
for i in range(len(points_list)):
	m = points_list[i][0]
	n = points_list[i][1]
	writefile.write('GLY   1A    N  ' + str(m)[:5] + '\t' + str(n)[:5] + '\t0.0 SR0  ' + '0.9\t0.0\t0.0\t1.0\n')


