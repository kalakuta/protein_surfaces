# converts map file output to list in a standard form
# for comparisons
a_ring = [(0,0)]
b_ring = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
c_ring = [(2,0),(1,1),(0,2),(-1,2),(-2,2),(-2,1),(-2,0),(-1,-1),(0,-2),(1,-2),(2,-2),(2,-1)]
d_ring = [(3,0),(2,1),(1,2),(0,3),(-1,3),(-2,3),(-3,3),(-3,2),(-3,1),(-3,0),(-2,-1),(-1,-2),(0,-3),(1,-3),(2,-3),(3,-3),(3,-2),(3,-1)]
e_ring = [(4,0),(3,1),(2,2),(1,3),(0,4),(-1,4),(-2,4),(-3,4),(-4,4),(-4,3),(-4,2),(-4,1),(-4,0),(-3,-1),(-2,-2),(-1,-3),(0,-4),(1,-4),(2,-4),(3,-4),(4,-4),(4,-3),(4,-2),(4,-1)]
f_ring = [(5,0),(4,1),(3,2),(2,3),(1,4),(0,5),(-1,5),(-2,5),(-3,5),(-4,5),(-5,5),(-5,4),(-5,3),(-5,2),(-5,1),(-5,0),(-4,-1),(-3,-2),(-2,-3),(-1,-4),(0,-5),(1,-5),(2,-5),(3,-5),(4,-5),(5,-5),(5,-4),(5,-3),(5,-2),(5,-1)]

def convert(map,cell_size):

	standardised_cells = {}
	for cell in map:
		a = cell[0] / cell_size
		b = cell[1] / cell_size
		new_cell = (int(a),int(b))
		standardised_cells[new_cell] = map[cell]
	
	score_list = []
	for cell in a_ring:
		score_list.append(standardised_cells[cell])
	for cell in b_ring:
		score_list.append(standardised_cells[cell])
	for cell in c_ring:
		score_list.append(standardised_cells[cell])
	for cell in d_ring:
		score_list.append(standardised_cells[cell])
	for cell in e_ring:
		score_list.append(standardised_cells[cell])
	for cell in f_ring:
		score_list.append(standardised_cells[cell])
	
	return(score_list)
		


