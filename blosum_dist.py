# functions for comparing amino acid make-up of two cells
# order A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V

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
'V': 19
}

blosum45 = [[5.0, -2.0, -1.0, -2.0, -1.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0, 0.0, -2.0, -2.0, 0.0],
[-2.0, 7.0, 0.0, -1.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 3.0, -1.0, -2.0, -2.0, -1.0, -1.0, -2.0, -1.0, -2.0],
[-1.0, 0.0, 6.0, 2.0, -2.0, 0.0, 0.0, 0.0, 1.0, -2.0, -3.0, 0.0, -2.0, -2.0, -2.0, 1.0, 0.0, -4.0, -2.0, -3.0],
[-2.0, -1.0, 2.0, 7.0, -3.0, 0.0, 2.0, -1.0, 0.0, -4.0, -3.0, 0.0, -3.0, -4.0, -1.0, 0.0, -1.0, -4.0, -2.0, -3.0],
[-1.0, -3.0, -2.0, -3.0, 12.0, -3.0, -3.0, -3.0, -3.0, -3.0, -2.0, -3.0, -2.0, -2.0, -4.0, -1.0, -1.0, -5.0, -3.0, -1.0],
[-1.0, 1.0, 0.0, 0.0, -3.0, 6.0, 2.0, -2.0, 1.0, -2.0, -2.0, 1.0, 0.0, -4.0, -1.0, 0.0, -1.0, -2.0, -1.0, -3.0],
[-1.0, 0.0, 0.0, 2.0, -3.0, 2.0, 6.0, -2.0, 0.0, -3.0, -2.0, 1.0, -2.0, -3.0, 0.0, 0.0, -1.0, -3.0, -2.0, -3.0],
[0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 7.0, -2.0, -4.0, -3.0, -2.0, -2.0, -3.0, -2.0, 0.0, -2.0, -2.0, -3.0, -3.0],
[-2.0, 0.0, 1.0, 0.0, -3.0, 1.0, 0.0, -2.0, 10.0, -3.0, -2.0, -1.0, 0.0, -2.0, -2.0, -1.0, -2.0, -3.0, 2.0, -3.0],
[-1.0, -3.0, -2.0, -4.0, -3.0, -2.0, -3.0, -4.0, -3.0, 5.0, 2.0, -3.0, 2.0, 0.0, -2.0, -2.0, -1.0, -2.0, 0.0, 3.0],
[-1.0, -2.0, -3.0, -3.0, -2.0, -2.0, -2.0, -3.0, -2.0, 2.0, 5.0, -3.0, 2.0, 1.0, -3.0, -3.0, -1.0, -2.0, 0.0, 1.0],
[-1.0, 3.0, 0.0, 0.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -3.0, 5.0, -1.0, -3.0, -1.0, -1.0, -1.0, -2.0, -1.0, -2.0],
[-1.0, -1.0, -2.0, -3.0, -2.0, 0.0, -2.0, -2.0, 0.0, 2.0, 2.0, -1.0, 6.0, 0.0, -2.0, -2.0, -1.0, -2.0, 0.0, 1.0],
[-2.0, -2.0, -2.0, -4.0, -2.0, -4.0, -3.0, -3.0, -2.0, 0.0, 1.0, -3.0, 0.0, 8.0, -3.0, -2.0, -1.0, 1.0, 3.0, 0.0],
[-1.0, -2.0, -2.0, -1.0, -4.0, -1.0, 0.0, -2.0, -2.0, -2.0, -3.0, -1.0, -2.0, -3.0, 9.0, -1.0, -1.0, -3.0, -3.0, -3.0],
[1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -3.0, -1.0, -2.0, -2.0, -1.0, 4.0, 2.0, -4.0, -2.0, -1.0],
[0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 2.0, 5.0, -3.0, -1.0, 0.0],
[-2.0, -2.0, -4.0, -4.0, -5.0, -2.0, -3.0, -2.0, -3.0, -2.0, -2.0, -2.0, -2.0, 1.0, -3.0, -4.0, -3.0, 15.0, 3.0, -3.0],
[-2.0, -1.0, -2.0, -2.0, -3.0, -1.0, -2.0, -3.0, 2.0, 0.0, 0.0, -1.0, 0.0, 3.0, -3.0, -2.0, -1.0, 3.0, 8.0, -1.0],
[0.0, -2.0, -3.0, -3.0, -1.0, -3.0, -3.0, -3.0, -3.0, 3.0, 1.0, -2.0, 1.0, 0.0, -3.0, -1.0, 0.0, -3.0, -1.0, 5.0]]

# takes two strings of amino acid cell make-up p and q and returns a similarity score
# based on the BLOSUM45 matrix
def blosum(p,q):
	a = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	b = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	p = (p.rstrip()).rstrip(',')
	q = (q.rstrip()).rstrip(',')

	aa_lista = p.split(',')
	for i in range(len(aa_lista)):
		aa_lista[i] = (aa_lista[i].lstrip()).split()
		aa = aa_lista[i][0]
		aa_index = aa_order[aa]
		aa_score = float(aa_lista[i][1])
		a[aa_index] += aa_score

	aa_listb = q.split(',')
	for i in range(len(aa_listb)):
		aa_listb[i] = (aa_listb[i].lstrip()).split()
		aa = aa_listb[i][0]
		aa_index = aa_order[aa]
		aa_score = float(aa_listb[i][1])
		b[aa_index] += aa_score

	score = 0.0
	# check for identical amino acids, score and then reduce each set of proportions by the overlapping amount

	non_id = 1.0
	for i in range(20):
		id = min(a[i],b[i])
		score += id * blosum45[i][i]
		a[i] = a[i] - id
		b[i] = b[i] - id
		non_id = non_id - id

	# score the non-identical parts proportionally

	for i in range(20):
		for j in range(20):
			score += (a[i] * b[j] * blosum45[i][j])/non_id
			
	return(score)
	





