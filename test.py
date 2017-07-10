import sys
import map_convert as mc
from time import sleep
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

	
map1 = {}
filename = sys.argv[1]
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

s = int(cell_size)
	# offset_list contains the translations of the base map which are run through in the alignment process
	#offset_list = [(-3,0)]
offset_list = [(x,y) for x in range(-2*s,3*s,s) for y in range(-2*s,3*s,s) if x+y < 3*s and x+y > -3*s]

for offset in offset_list:
	offset_map1 = {}
	for cell in map1:
		a = cell[0] + offset[0]
		b = cell[1] + offset[1]
		offset_cell = (a,b)
		offset_map1[offset_cell] = map1[cell]
	map_list1 = mc.convert(offset_map1,cell_size)
	print(map_list1)

	
