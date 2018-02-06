import sys
import random
mapfile = sys.argv[1]
txt = open(mapfile,'r')
line = txt.readline()
pdb, centre, chain, aa, atom = line.split()

coords3D = {}
coords2D = {}
paths = {}

while True:
	line = txt.readline().split()
	if not line: break
	coords3D[line[0]] = [line[1], line[2], line[3]]
	coords2D[line[0]] = [line[4], line[5]]
	paths[line[0]] = line[9:]

txt.close()

point = random.choice(list(coords3D.keys()))
print(paths[point])
file = open('./path_files/%s_%s_%s.kin' % (pdb, centre, point), 'w+')
file.write('@text\n')
file.write('This file shows part of the surface of ' + pdb + ' centred at the point ' + centre + '(' + chain + ', ' + aa + ', ' + atom + ') and the path to point ' + point + '.\n\n')
file.write('@kinemage\n\n')

#write all surface points to file
file.write('@dotlist {whole_surface} color=gray\n')
for key in coords3D:
	file.write(coords3D[key][0] + ' ' + coords3D[key][1] + ' ' + coords3D[key][2] + '\n')
	
#write path from centre to point 
file.write('@vectorlist {path_to_point} color=green\n')
file.write(coords3D[paths[point][0]][0] + ' ' + coords3D[paths[point][0]][1] + ' ' + coords3D[paths[point][0]][2] + ' ') 
length = len(paths[point])
for i in range(length - 1):
	file.write(coords3D[paths[point][i + 1]][0] + ' ' + coords3D[paths[point][i + 1]][1] + ' ' + coords3D[paths[point][i + 1]][2] + '\n')
	file.write(coords3D[paths[point][i + 1]][0] + ' ' + coords3D[paths[point][i + 1]][1] + ' ' + coords3D[paths[point][i + 1]][2] + ' ')
file.write(coords3D[paths[point][length - 1]][0] + ' ' + coords3D[paths[point][length - 1]][1] + ' ' + coords3D[paths[point][length - 1]][2] + '\n')







