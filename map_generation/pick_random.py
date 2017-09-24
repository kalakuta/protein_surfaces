import os
import sys
import random

pdb_list = os.listdir(sys.argv[1])

for i in range(5):
	print(random.choice(pdb_list))
