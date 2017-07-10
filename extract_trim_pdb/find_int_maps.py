import sys
import os


atlas = sys.argv[1]
map_list = os.listdir(atlas)

for filename in map_list:
	file = open(atlas + filename, "r")
	file.readline()
	empty_cells = 0
	while True:
		line = file.readline()
		if not line: break
		aa = line.split(':')[7].strip()
		if aa == 'NULL':
			empty_cells += 1
	if empty_cells > 50:
		print(filename)

