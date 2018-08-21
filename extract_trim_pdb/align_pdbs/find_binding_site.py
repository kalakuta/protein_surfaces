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
			if line[11:14] == site:
				records = line[18:]
				for i in range(len(records) // 11):
					record = records[i * 11: i * 11 + 11]
					if record != "           ":
						aa = record[0:3]
						seq = record[5:9].strip() + record[4] 
					else: continue
					if aa == 'HOH': continue
					if site not in site_lists:
						site_lists[site] = [[aa, seq]]
					else:
						site_lists[site].append([aa, seq])

for key in site_lists:
	print(key)
	for i in range(len(site_lists[key])):
		print(site_lists[key][i][0], site_lists[key][i][1])


		
	

