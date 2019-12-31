file_input="dump.seq_e10.mci.I12"
num=1
output=open('assign_cluster_num.out','w')
for i in open(file_input):
	i=i.replace('\n','')
	o=i.split()
	
	for j in o:
		output.write(j+'\t'+str(num)+'\n')

	num+=1
output.close()
