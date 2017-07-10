import sys
import os
from numpy import random as rand
from curvature_distance import distance as cu
from hydro_distance import distance as hy
from h_acc_distance import distance as ha
from h_don_distance import distance as hd
from aa_distance import distance as aa

map_list = os.listdir('/d/mw4/u/mp001/non_homologous_set/maps')
file = open('/d/mw4/u/mp001/non_homologous_set/distances.txt', 'w+')
file.write('map1, map2, rotation, hydro_distance, h_acc_distance, h_don_distance, curvature_distance, aa_distance\n')
#run through all pairs of distinct maps in non_homologous_set
for i in range(len(map_list)):
	for j in range(i+1, len(map_list)):
		rotation = 15 * rand.randint(0, 23)
		map1 = '/d/mw4/u/mp001/non_homologous_set/maps/' +  map_list[i]
		map2 = '/d/mw4/u/mp001/non_homologous_set/maps/' +  map_list[j]
		file.write(map_list[i] + ',' + map_list[j] + ',')
		file.write(str(rotation) + ',')
		file.write(str(hy(map1, map2, rotation)) + ',')
		file.write(str(ha(map1, map2, rotation)) + ',')
		file.write(str(hd(map1, map2, rotation)) + ',')
		file.write(str(cu(map1, map2, rotation)) + ',')
		file.write(str(aa(map1, map2, rotation)))
		file.write('\n')
file.close()

