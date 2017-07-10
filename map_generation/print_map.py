import sys
import numpy as np
import matplotlib.pyplot as plt


#generate 2d plot using polar coordinates: angle is original direction, length is geodesic distance


map_file = sys.argv[1]
domain = map_file.split('\\')[-1][:-4]		# strip '.map' from filename to give pdb domain
txt = open(map_file,"r")
cell_size = float(txt.readline().split()[2])

map = {}
while True:
	line = txt.readline()
	if not line: break
	c = line.split(':')[0]
	c = c.strip('(')
	c = c.strip()
	c = c.strip(')')
	cell = (int(c.split(', ')[0]),int(c.split(', ')[1]))
	value = float(line.split(':')[2])
	map[cell]=value


'''
# plot individual points
plt.plot()
plt.xlim(-20,20)
plt.ylim(-20,20)
plt.gca().set_aspect('equal',adjustable='box')

for point in connected_points:
	col = 'green'
	if hydro_scores[point] < 0:
		col = 'red'
	p = geo_coords[point][0]
	q = geo_coords[point][1]
	plt.plot(p,q,'ro',c = col)



#plt.show()
plt.savefig('./runs/2D_map%s_%s.png' % (centre_point,rotation_angle))
plt.close()
'''

#plot hexagonal grid, coloured by hydrophobicity score for cell

plt.plot()
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.gca().set_aspect('equal',adjustable='box')
plt.suptitle(domain)


for cell in map:
	g = min(1,(0.5 + map[cell] / 25))
	if g < 0.0:
		g = 0.0
	r = 1 - g
	p = cell[0] + 0.5 * cell[1]
	q = cell[1] * np.sqrt(3) / 2
	plt.plot(p,q,'rh',markersize = 6.5 * cell_size,c = (r,g,0))



plt.show()
#plt.savefig('./runs/hex_map%s.png' % domain)