#! /usr/bin/env python3
# selects pdb domains from list and crops the relevant pdb files,
# recording only ATOM records in the principle conformation where alternatives given,
# excluding hydrogen atoms

import sys

filename = sys.argv[1]
pdb_file = open(filename,"r")
lines = pdb_file.readlines()
pdb_file.close()
domain = filename[0:4]
print(domain) 
outfile = '/d/mw4/u/mp001/' + domain + '.pdb'
domainpdb = open(outfile,"w+")
for text in lines:
	if text.split()[0] == 'ATOM':	# select only ATOM records
		if text[77] == 'H':			# exlude hydrogen atoms
			continue
		if text[16] != ' ':
			if text[16] != 'A':
				continue
		domainpdb.write(text)
domainpdb.close()
