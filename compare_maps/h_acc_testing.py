from h_acc_distance import distance as d
import sys
from p_scores import h_acc_p_score as p
filename = sys.argv[1]
file = open(filename,"r")
file.readline()
distributions = {}
while True:
	line = file.readline()
	if not line: break
	map1 = line.split(',')[0]
	map2 = line.split(',')[1]
	score = p(float(line.split(',')[3]))
	if map1 not in distributions:
		distributions[map1] = score / 658
	else:
		distributions[map1] += score /658
	if map2 not in distributions:
		distributions[map2] = score / 658
	else:
		distributions[map2] += score /658

for i in range(100):
	j = i / 100
	count = 0
	for map in distributions:
		if distributions[map] < j + 0.01:
			count += 1
	print(j, count)

for map in distributions:
		if distributions[map] < 0.26:
			print(map)
