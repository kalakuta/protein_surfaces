# reads two cell_scores files and finds their hydro_distance score minimised over
# rotations of 0, 15, 30, 45, ...,345 and translations to two rings out from the centre

import math
import numpy as np
import sys
import map_convert as mc
import rotate as rt



def distance(a,b):
	# dist(x,y) calculates the mean absolute difference between values in lists a and b
	def dist(x,y):
		d = 0.0
		n = len(x)
		for i in range(n):
			d += (np.abs(x[i] - y[i])) / n
		return(d)


	# convert input files to list of values
		
	map1 = {}
	filename = a
	txt = open(filename,"r")
	line = txt.readline()
	cell_size = float(line.split()[2])
	while True:
		line = txt.readline()
		if not line: break
		c = line.split(':')[0]
		c = c.strip('(')
		c = c.strip()
		c = c.strip(')')
		cell = (int(c.split(', ')[0]),int(c.split(', ')[1]))
		hydro_value = float(line.split(':')[2])
		map1[cell] = hydro_value

	txt.close()

	map2 = {}
	filename = b
	txt = open(filename,"r")
	line = txt.readline()
	cell_size2 = float(line.split()[2])
	if cell_size2 != cell_size:
		sys.exit('cell sizes must be equal')
	while True:
		line = txt.readline()
		if not line: break
		c = line.split(':')[0]
		c = c.strip('(')
		c = c.strip()
		c = c.strip(')')
		cell = (int(c.split(', ')[0]),int(c.split(', ')[1]))
		hydro_value = float(line.split(':')[2])
		map2[cell] = hydro_value

	txt.close()

	map_list2 = mc.convert(map2,int(cell_size))

	min_score = 10000.0
	min_angle = 0
	min_offset = (0, 0)

	s = int(cell_size)
	# offset_list contains the translations of the base map which are run through in the alignment process
	#offset_list = [(0, 0)]
	offset_list = [(x,y) for x in range(-2*s,3*s,s) for y in range(-2*s,3*s,s) if x+y < 3*s and x+y > -3*s]

	for offset in offset_list:
		offset_map1 = {}
		for cell in map1:
			a = cell[0] + offset[0]
			b = cell[1] + offset[1]
			offset_cell = (a,b)
			offset_map1[offset_cell] = map1[cell]
			
		map_list1 = mc.convert(offset_map1,s)

		for i in range(0, 360, 15):
			(q, r) = divmod(i,60)
			base_map = map_list1
			# rotate base_map by 60 q times
			for j in range(q):
				base_map = rt.rotate60(base_map)
			# rotate base_map by r degress
			if r == 15:
				base_map = rt.rotate15(base_map)
			if r == 30:
				base_map = rt.rotate30(base_map)
			if r == 45:
				base_map = rt.rotate45(base_map)
			m1 = base_map[:61]
			m2 = map_list2[:61]
			d = dist(m1,m2)
			if d < min_score:
				min_score = d
				min_angle = i
				min_offset = offset

	return([min_score, min_angle, min_offset])


	



