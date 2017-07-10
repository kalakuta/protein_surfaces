import sys
from curvature_distance import distance as cd
from hydro_distance import distance as hd
from h_acc_distance import distance as ad
from h_don_distance import distance as dd
from aa_distance import distance as ad

a = sys.argv[1]
b = sys.argv[2]
print(hd(a,b))
print(ad(a,b))
print(dd(a,b))
print(cd(a,b))
print(dd(a,b))