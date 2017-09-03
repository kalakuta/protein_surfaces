import os
import sys

for filename in os.listdir(sys.argv[1]):
	pdb_code = filename[:-10]
	cmd = 'dms ' + sys.argv[1] + filename + ' -n -a -d 0.2 -o ' + 'dms/' + pdb_code + '.dms'
	os.system(cmd)
	#print(cmd)
