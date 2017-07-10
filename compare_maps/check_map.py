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


map1 = sys.argv[1]
map2 = sys.argv[2]
atlas_path = sys.argv[2]
#offset_list = [(x,y) for x in range(-2,3) for y in range(-2,3) if x+y < 3 and x+y > -3]
offset_list = [(0, 0)]

filename = 'map_check_aa_' + map1[-9:-4] + map2[-14:-4] + '.csv'
file = open(filename, "w+")

for offset in offset_list:
	for rotation in range(0, 360, 15):
		#hy_score = hy(map2 , map1, offset, rotation)
		#ha_score = ha(map2 , map1, offset, rotation)
		#hd_score = hd(map2 , map1, offset, rotation)
		#cu_score = cu(map2 , map1, offset, rotation)
		aa_score = aa(map2 , map1, offset, rotation)
		#if p.hydro_p_score(hy_score) != 0.0: #0.0 indicates a perfect score, but error thrown on log odds score
		#	p_hy = -np.log(p.hydro_p_score(hy_score))
		#else:
		#	p_hy = 1000.0 #very large value - exceeding any obtained for non-exact matches
		
		#if p.h_acc_p_score(ha_score) != 0.0: #0.0 indicates a perfect score
		#	p_ha = -np.log(p.h_acc_p_score(ha_score))
		#else:
		#	p_ha = 1000.0 #very large value - exceeding any obtained for non-exact matches
		
		#if p.h_don_p_score(hd_score) != 0.0: #0.0 indicates a perfect score
		#	p_hd = -np.log(p.h_don_p_score(hd_score))
		#else:
		#	p_hd = 1000.0 #very large value - exceeding any obtained for non-exact matches
		
		#if p.curvature_p_score(cu_score) != 0.0: #0.0 indicates a perfect score
		#	p_cu = -np.log(p.curvature_p_score(cu_score))
		#else:
		#	p_cu = 1000.0 #very large value - exceeding any obtained for non-exact matches
		
		if p.aa_p_score(aa_score) != 0.0: #0.0 indicates a perfect score
			p_aa = -np.log(p.aa_p_score(aa_score))
		else:
			p_aa = 1000.0 #very large value - exceeding any obtained for non-exact matches
		
		file.write('aa,' + str(p_aa) + ',' + str(rotation) + '\n')
file.close()		