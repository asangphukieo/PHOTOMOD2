import sys

def check_value_correct(file1,file2,file_combind):
	sum_file1=0
	for i in open(file1):
		if '#' not in i:
			i=i.replace('\n','')
			line1=i.split(',')
			del line1[0]
			del line1[len(line1)-1]
			
			sum_line1=0
			for j in line1:
				sum_line1+=int(j)
			
			sum_file1+=sum_line1

	sum_file2=0
	for i in open(file2):
		if '#' not in i:
			i=i.replace('\n','')
			line2=i.split(',')
			del line2[0]
			del line2[len(line2)-1]
			
			sum_line2=0
			for j in line2:
				sum_line2+=int(j)

			sum_file2+=sum_line2

	print 'File 1 => ',sum_file1
	print 'File 2 => ',sum_file2

	sum_file3=0
	for i in open(file_combind):
		if '#' not in i:
			i=i.replace('\n','')
			line3=i.split(',')
			del line3[0]
			del line3[len(line3)-1]
			
			sum_line3=0
			for j in line3:
				sum_line3+=int(j)

			sum_file3+=sum_line3

	print 'Sum of file1 and file2 => ',(sum_file1+sum_file2)
	print 'File combind => ',sum_file3
	if int(sum_file1+sum_file2) != int(sum_file3):
		print "FIle input and ouput are not equal!!. Stop process."
		sys.exit()

check_value_correct(sys.argv[1],sys.argv[2],sys.argv[3])
