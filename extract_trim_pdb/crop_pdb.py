#! /usr/bin/env python3
# selects pdb domains from list and crops the relevant pdb files,
# recording only ATOM records in the principle conformation where alternatives given,,
# excluding hydrogen atoms

import sys

filename = sys.argv[1]
domain_file = open(filename,"r")
while True:
        line = domain_file.readline()
        if not line: break
        domain = line.split()[0]
        chain = domain[4]
        residues = line.split(None,1)[1]
        residues = residues.split()
        residues_list = []
        for i in range(len(residues)):
                a,b = int(residues[i].split(':')[0]),int(residues[i].split(':')[1])
                for j in range(a,b+1):
                        residues_list.append(j)
        pdb = domain[:4]
        filename = '/d/mw4/u/mp001/pdb/pdb' + pdb + '.ent'
        pdb_file = open(filename,"r")
        lines = pdb_file.readlines()
        pdb_file.close()
        filename = '/d/mw4/u/mp001/non_homologous_set/pdb/' + domain + '.pdb'
        domainpdb = open(filename,"w+")
        for text in lines:
                if text.split()[0] == 'ATOM':
                        if text[14] == 'H':
                                continue
                        if text[16] != ' ':
                                if text[16] != 'A':
                                        continue
                        if text.split()[4] == chain:
                                try:
                                        y = int(text.split()[5])
                                except:
                                        y = 'blah'
                                if y in residues_list:
                                        domainpdb.write(text)
                        elif text.split()[4][0] == chain:
                                try:
                                        y = int(text.split()[4][1:])
                                except:
                                        y = 'blah'
                                if y in residues_list:
                                        domainpdb.write(text)
        domainpdb.close()

domain_file.close()