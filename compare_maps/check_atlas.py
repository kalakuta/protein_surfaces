# uses distance.py to run through each file in atlas_path
# to find best match

import sys
import os
import numpy as np
import p_scores as p
from curvature_distance import distance as cu
from hydro_distance import distance as hy
from h_acc_distance import distance as ha
from h_don_distance import distance as hd
from aa_distance import distance as aa




map1_path = sys.argv[1]
map1_list = os.listdir(map1_path)
atlas_path = sys.argv[2]
offset_list = [(x,y) for x in range(-2,3) for y in range(-2,3) if x+y < 3 and x+y > -3]
for map1 in map1_list:	# check all maps in map1_path
	map1 = map1_path + map1
	txt = open(map1, "r")
	line1 = txt.readline().split()
	txt.close()
	filename = 'atlas_check' + map1[-9:-4] + '.csv'
	file = open(filename, "w+")
	for item in line1:
		file.write(item + ",")
	file.write("\n")
	file.write("\n")
	for map in os.listdir(atlas_path):
		txt = open(atlas_path + map, "r")
		atom = txt.readline().split()[3:6] # residue number, aa type, atom
		txt.close()
		for offset in offset_list:
			for rotation in range(0, 360, 15):
				hy_score = hy(atlas_path + map , map1, offset, rotation)
				ha_score = ha(atlas_path + map , map1, offset, rotation)
				hd_score = hd(atlas_path + map , map1, offset, rotation)
				cu_score = cu(atlas_path + map , map1, offset, rotation)
				aa_score = aa(atlas_path + map , map1, offset, rotation)
				if p.hydro_p_score != 0.0: #0.0 indicates a perfect score, but error thrown on log odds score
					p_hy = -np.log(p.hydro_p_score(hy_score))
				else:
					p_hy = 1000.0 #very large value - exceeding any obtained for non-exact matches
				if p.h_acc_p_score != 0.0: #0.0 indicates a perfect score
					p_ha = -np.log(p.h_acc_p_score(ha_score))
				else:
					p_ha = 1000.0 #very large value - exceeding any obtained for non-exact matches
				if p.h_don_p_score != 0.0: #0.0 indicates a perfect score
					p_hd = -np.log(p.h_don_p_score(hd_score))
				else:
					p_hd = 1000.0 #very large value - exceeding any obtained for non-exact matches
				if p.curvature_p_score != 0.0: #0.0 indicates a perfect score
					p_cu = -np.log(p.curvature_p_score(cu_score))
				else:
					p_cu = 1000.0 #very large value - exceeding any obtained for non-exact matches
				if p.aa_p_score != 0.0: #0.0 indicates a perfect score
					p_aa = -np.log(p.aa_p_score(aa_score))
				else:
					p_aa = 1000.0 #very large value - exceeding any obtained for non-exact matches
				overall_score = p_hy + p_ha + p_hd + p_cu
				if overall_score > 5.0:
						file.write(str(p_hy) + ',' + str(p_ha) + ',' + str(p_hd) + ',' + str(p_cu) + ',' + str(p_aa) + ',' + str(overall_score) + ',' + str(map) + ',' + str(offset) + ',' + str(rotation) + ',' + atom[0] + ',' + atom [1] + ',' + atom[2] + '\n')
	file.close()		