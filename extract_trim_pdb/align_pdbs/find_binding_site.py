'''
Find the heme binding residues from a pdb file

'''

import sys
import os


pdb_file = open(sys.argv[1],"r")
lines = pdb_file.readlines()
pdb_file.close()


site_identifiers = []

for i in range(len(lines)):
	if lines[i][0:10] != "REMARK 800":
		continue
	if len(lines[i].split()) < 3:
		continue
	else:
		if lines[i].split()[2] == "SITE_IDENTIFIER:":
			site_identifier = lines[i].split()[3]
			for j in range(3):
				line=lines[i+j]
				if " HEM " in line and "SITE_DESCRIPTION" in line:
					site_identifiers.append(site_identifier)

site_lists = {}
			
for site in site_identifiers:
	for line in lines:
		if line.split()[0] == "SITE":
			if line.split()[2] == site:
				if site not in site_lists:
					site_lists[site] = line.split()[4:]
				else:
					site_lists[site] = site_lists[site] + line.split()[4:]

'''
make format match .map files
'''			

		
	

