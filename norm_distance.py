# reads two cell_scores files, best aligns the second with the first
# and finds their norm_distance score minimised over
# rotations of 0, 15, 30, 45, ...,345 and translations to two rings out from the centre

import math
import numpy as np
import sys
import map_convert as mc
import rotate_norm as rt
import align_normals as al

# dot(x,y) calculates the mean of the dot products between vectors in lists x and y
def dot(x,y):
	d = 0.0
	n = len(x)
	for i in range(n):
		d += (np.dot(x[i],y[i])) / n
	return(d)

def distance(a,b):
	
	# convert input files to lists of normal vectors
		
	map1 = {}
	filename = a
	txt = open(filename,"r")
	line = txt.readline()
	cell_size = float(line.split()[2])
	while True:
		line = txt.readline()
		if not line: break
		cell_count = int(line.split(':')[1].strip())
		c = line.split(':')[0]
		c = c.strip()
		c = c.strip('(')
		c = c.strip(')')
		cell = (int(c.split(',')[0]),int(c.split(',')[1]))
		normal_string = line.split(':')[5]
		if cell_count == 0:
			normal = [0.0,0.0,1.0]
		else:	
			normal = [float(normal_string.split(',')[0]), float(normal_string.split(',')[1]), float(normal_string.split(',')[2])]
		map1[cell] = normal

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
		cell_count = int(line.split(':')[1].strip())
		c = line.split(':')[0]
		c = c.strip()
		c = c.strip('(')
		c = c.strip(')')
		cell = (int(c.split(', ')[0]),int(c.split(', ')[1]))
		normal_string = line.split(':')[5]
		if cell_count == 0:
			normal = [0.0,0.0,1.0]
		else:	
			normal = [float(normal_string.split(',')[0]), float(normal_string.split(',')[1]), float(normal_string.split(',')[2])]
		map2[cell] = normal

	txt.close()
	# find the best alignment for sets of normals and rotate map2 by this alignment
	
	

	max_score = 0.0
	max_angle = 0
	max_offset = (0, 0)

	s = int(cell_size)
	
	
	# offset_list contains the translations of the base map which are run through in the alignment process
	offset_list = [(x,y) for x in range(-2*s,3*s,s) for y in range(-2*s,3*s,s) if x+y < 3*s and x+y > -3*s]
	
	for offset in offset_list:			# run through possible translations
		offset_map1 = {}
		for cell in map1:
			x = cell[0] + offset[0]
			y = cell[1] + offset[1]
			offset_cell = (x, y)
			offset_map1[offset_cell] = map1[cell]
		
		map_list1 = mc.convert(offset_map1,s)
		map_list2 = mc.convert(map2,s)	

		for i in range(0, 360, 15):		# run through possible rotations
			(q, r) = divmod(i,60)
			base_map = map_list1
			# rotate base_map by 60 degrees q times
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
			# get alignment matrix to best rotate m2 onto m1
			R = al.align(m1,m2)
			for k in range(len(m2)):
				m2[k] = np.dot(R,m2[k])
			d = dot(m1,m2)
			if d > max_score:
				max_score = d
				max_angle = i
				max_offset = offset

	return([max_score, max_angle, max_offset])


	



