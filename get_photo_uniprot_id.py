import os

input_file1='get_photo_gene.out2'
input_file2="blast_photo_cutoff.out"

os.system('rm get_photo_gene.out3')
for j in open(input_file1):
	j=j.replace('\n','')
	os.system('grep "'+j+'" '+input_file2+' >> get_photo_gene.out3')

for i in open('get_photo_gene.out3'):
	i=i.replace('\n','')
	col=i.split('\t')
	prot_id=col[1].split('|')
	print col[0]+'\t'+prot_id[1]
