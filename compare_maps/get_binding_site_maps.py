import os

check_list = [25, 28, 36, 40, 55, 57, 59, 68, 70, 86, 88, 105, 107, 123, 125, 128, 130, 133]

for filename in os.listdir('1x8qA'):
	txt = open(('1x8qA/' + filename), 'r')
	line = txt.readline().split('\t')
	base = int(line[3][:-1])
	type = line[4]
	if base in check_list:
		print(type, base, filename)
	txt.close()