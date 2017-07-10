'''
Reads a pdb file and removes all but ATOM records

'''

import os

pdb_list = os.listdir('rotated_pdbs')
for pdbB in pdb_list:
	file = open('rotated_pdbs/' + pdbB, 'r')
	pdb_b = file.readlines()
	file.close()

	# write new .pdb file
	pdb_code = os.path.basename(pdbB[:5])
	outfile = 'final_pdbs/' + pdb_code + '_final.pdb'
	f = open(outfile, 'w+')
	for line in pdb_b:
		if line[0:6].strip() == 'ATOM':
			f.write(line)
	f.close()

	

