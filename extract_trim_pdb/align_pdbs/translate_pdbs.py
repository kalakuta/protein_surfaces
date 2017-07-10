import numpy as np
import sys
import os

pdb_list = os.listdir('cropped_pdbs')
for pdb in pdb_list:	
	file = open('cropped_pdbs/' + pdb,'r')
	lines = file.readlines()
	file.close()

	# get coordinates of heme's nitrogens
	for line in lines:
		if line[0:6] == 'HETATM' and line[17:20] == 'HEM':	
			if line[12:16].strip() == 'NA':
				NA = np.asarray([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
				#print('NA:', NA)
			if line[12:16].strip() == 'NB':
				NB = np.asarray([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
				#print('NB:', NB)
			if line[12:16].strip() == 'NC':
				NC = np.asarray([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
				#print('NA:', NA)
			if line[12:16].strip() == 'ND':
				ND = np.asarray([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
				#print('NA:', NA)

	centroid = (NA + NB + NC + ND) / 4

	# write new .pdb file
	pdb_code = os.path.basename(pdb[:-4])
	outfile = pdb_code + '_trans.pdb'
	f = open(outfile, 'w+')
	for line in lines:
		if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM':
			coords = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
			coords = np.asarray(coords) - np.asarray(centroid)
			entry = line[0:30]
			for coordinate in coords:
				entry += ('        ' + str(np.round(coordinate, 3)))[-8:]
			entry += line[54:]
			f.write(entry)
		else:
			f.write(line)
	f.close()

