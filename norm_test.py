import norm_distance as nd
import sys
from os import listdir


b = sys.argv[1]
max_d = 0.0
for a in listdir('./runs/tosend/run8'):
	x = './runs/tosend/run8/' + str(a)
	try:
		dist = nd.distance(x, b)[0]
		if  dist > max_d:
			max_d = nd.distance(x, b)[0]
		print(a, nd.distance(x, b))
	except:
		continue
print(max_d)