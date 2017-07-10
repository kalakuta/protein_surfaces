import os

for filename in os.listdir('/d/mw4/u/mp001/test'):
	pdb_code = filename[:-4]
	base_dir = '/d/mw4/u/mp001/test/'
	cmd = '/d/user5/mp001/dms/bin/dms ' + base_dir + filename + ' -n -a -d 0.2 -o ' + base_dir + 'dms/' + pdb_code + '.dms'
	os.system(cmd)
	#print(cmd)
