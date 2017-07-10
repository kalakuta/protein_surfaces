import matplotlib.pyplot as plt
import math


plt.plot()
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.gca().set_aspect('equal',adjustable='box')
#plt.suptitle(domain)
s = 3
cell_list = [(x,y) for x in range(-5*s,6*s,s) for y in range(-5*s,6*s,s) if x+y < 6*s and x+y > -6*s]

for cell in cell_list:
	
	p = cell[0] + 0.5 * cell[1]
	q = cell[1] * math.sqrt(3) / 2
	plt.plot(p,q,'h',markersize=22, color=(0.9,0.9,1))



plt.show()
#plt.savefig('./runs/hex_map%s.png' % domain)