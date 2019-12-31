import sys
import commands
import pandas as pd

def check_value_correct(file1,file2,file_combind):
	sum_file1=0
	for i in open(file1):
		if '.C.' not in i:
			i=i.replace('\n','')
			line1=i.split(',')
			#del line1[0]
			del line1[len(line1)-1]
			
			sum_line1=0
			for j in line1:
				sum_line1+=int(j)
			
			sum_file1+=sum_line1

	sum_file2=0
	for i in open(file2):
		if '.C.' not in i:
			i=i.replace('\n','')
			line2=i.split(',')
			#del line2[0]
			del line2[len(line2)-1]
			
			sum_line2=0
			for j in line2:
				sum_line2+=int(j)

			sum_file2+=sum_line2

	print 'File 1 => ',sum_file1
	print 'File 2 => ',sum_file2

	sum_file3=0
	for i in open(file_combind):
		if '.C.' not in i:
			i=i.replace('\n','')
			line3=i.split(',')
			#del line3[0]
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


def combine_feature_2file(file1,file2):
	#df1=pd.read_csv(file1, header=None)
	#df2=pd.read_csv(file2, header=None)

	df1=pd.read_csv(file1)
	df2=pd.read_csv(file2)
	print df1

	#df1=df1.astype('str')
	#df2=df2.astype('str')

	label1=df1.iloc[:,-1] #copy row of label in the last column
	label2=df2.iloc[:,-1]

	label_combine=pd.concat([label1,label2])
	#print label_combine

	#gene1=df1.iloc[:,0] #copy row of label in the first column
	#gene2=df2.iloc[:,0]

	#gene_combine=pd.concat([gene1,gene2])


	df1 = df1.iloc[:, :-1] #remove last column
	df2 = df2.iloc[:, :-1]

	#df1 = df1.iloc[:, 1:] #remove first column
	#df2 = df2.iloc[:, 1:]


	df1=df1.fillna(0).astype("int")
	df2=df2.fillna(0).astype("int")


	df_combine=df1.append(df2)

	df_combine=df_combine.fillna(0)
	df_combine=df_combine.astype("int")

	#df_combine=pd.concat([gene_combine, df_combine], axis=1)

	df_combine=pd.concat([df_combine, label_combine], axis=1)



	#print df_combine
	df_combine.to_csv('./comb.temp', index=False)
	
combine_feature_2file(sys.argv[1],sys.argv[2])
check_value_correct(sys.argv[1],sys.argv[2],'comb.temp')
