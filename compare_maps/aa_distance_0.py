# functions for comparing amino acid make-up of two cells
# order A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
import rotate_norm as rt
import map_convert as mc
import sys

aa_order = {
'A': 0,
'R': 1,
'N': 2,
'D': 3,
'C': 4,
'Q': 5,
'E': 6,
'G': 7,
'H': 8,
'I': 9,
'L': 10,
'K': 11,
'M': 12,
'F': 13,
'P': 14,
'S': 15,
'T': 16,
'W': 17,
'Y': 18,
'V': 19,
'X': 20
}

def aa_list(p):
	a = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	p = (p.rstrip()).rstrip(',')

	aa_lista = p.split(',')
	for i in range(len(aa_lista)):
		aa_lista[i] = (aa_lista[i].lstrip()).split()
		aa = aa_lista[i][0]
		if aa == 'NULL':
			a[20] = 1.0
		else:
			aa_index = aa_order[aa]
			aa_score = float(aa_lista[i][1])
			a[aa_index] += aa_score
	return(a)

	
def identity_score(a_list,b_list):
	score = 0.0
	# check for identical amino acids, score and then reduce each set of proportions by the overlapping amount
	for i in range(len(a_list)):
		for j in range(21):
			score += min(a_list[i][j], b_list[i][j])

	return(score)
	

def distance(a,b):
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
		aa_makeup = line.split(':')[7]
		map1[cell] = aa_list(aa_makeup)

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
		aa_makeup = line.split(':')[7]
		map2[cell] = aa_list(aa_makeup)

	txt.close()

	map_list2 = mc.convert(map2,int(cell_size))

	max_score = 0.0
	max_angle = 0
	max_offset = (0, 0)

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
			
		map_list1 = mc.convert(offset_map1,cell_size)

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
			d = identity_score(m1,m2)
			if d > max_score:
				max_score = d
				max_angle = i
				max_offset = offset

	return([max_score, max_angle, max_offset])


	





