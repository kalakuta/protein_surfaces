import sys
txt = open(sys.argv[1],'r')
i=0
j=0
points = {}
while True:
	line = txt.readline()
	i += 1
	if not line: break
	if 'A\n' not in line:
		a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11 = line.split()
		points[j] = [i, [float(a4),float(a5),float(a6)]]
		j += 1
for j in points:
	if points[j][0] == 1688:
		print(str(j), points[j])
	if points[j][0] == 1689:
		print(str(j), points[j])
	if points[j][0] == 4053:
		print(str(j), points[j])
	
