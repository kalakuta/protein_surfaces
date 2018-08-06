import sys
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



mapfile = sys.argv[1]
txt = open(mapfile,'r')
line = txt.readline()
pdb, centre, chain, aa, atom = line.split()[:5]
normal = [float(x) for x in line.split()[-3:]]
print(normal)

coords3D = {}
coords2D = {}
paths = {}

while True:
	line = txt.readline().split()
	if not line: break
	coords3D[line[0]] = [line[1], line[2], line[3]]
	coords2D[line[0]] = [float(line[4]), float(line[5])]
	paths[line[0]] = line[9:]

txt.close()

point = random.choice(list(coords3D.keys()))
#point = '6669'
path = paths[point]
file = open('./path_files/%s_%s_%s.kin' % (pdb, centre, point), 'w+')
file.write('@text\n')
file.write('This file shows part of the surface of ' + pdb + ' centred at the point ' + centre + '(' + chain + ', ' + aa + ', ' + atom + ') and the path to point ' + point + '.\n\n')
file.write('@kinemage\n\n')

#write all surface points to file
file.write('@balllist {whole_surface} color=black radius = 0.1\n')
for key in coords3D:
	file.write(coords3D[key][0] + ' ' + coords3D[key][1] + ' ' + coords3D[key][2] + '\n')
#highlight centre_point
file.write('@balllist {centre_point} color=red radius=0.2\n')
file.write(coords3D[centre][0] + ' ' + coords3D[centre][1] + ' ' + coords3D[centre][2] + '\n')
#highlight points on path
file.write('@balllist {path points} color=blue radius = 0.2\n')
length = len(path)
for i in range(length - 1):
	file.write(coords3D[path[i + 1]][0] + ' ' + coords3D[path[i + 1]][1] + ' ' + coords3D[path[i + 1]][2] + '\n')
#write path from centre to point
file.write('@vectorlist {path_to_point} color=green\n')
file.write(coords3D[path[0]][0] + ' ' + coords3D[path[0]][1] + ' ' + coords3D[path[0]][2] + ' ')
for i in range(length - 1):
	file.write(coords3D[path[i + 1]][0] + ' ' + coords3D[path[i + 1]][1] + ' ' + coords3D[path[i + 1]][2] + '\n')
	file.write(coords3D[path[i + 1]][0] + ' ' + coords3D[path[i + 1]][1] + ' ' + coords3D[path[i + 1]][2] + ' ')
file.write(coords3D[path[length - 1]][0] + ' ' + coords3D[path[length - 1]][1] + ' ' + coords3D[path[length - 1]][2] + '\n')
#write normal vector
file.write('@vectorlist {centre normal} color=red\n')
file.write(coords3D[centre][0] + ' ' + coords3D[centre][1] + ' ' + coords3D[centre][2])
for i in range(3):
	file.write(' ' + str(float(coords3D[centre][i]) + normal[i]))
file.write('\n')

file.close()

path_coords = [coords2D[item] for item in path]
x = [item[0] for item in path_coords]
y = [item[1] for item in path_coords]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_xlim(-20,20)
ax1.set_ylim(-20,20)
ax1.scatter(x, y)
ax1.plot(x, y)
#fig1.show()
fig1.savefig('./path_files/%s_%s_%s.png' % (pdb, centre, point))





