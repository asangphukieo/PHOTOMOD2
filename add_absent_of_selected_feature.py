import sys
import commands
import pandas as pd
#using first file as feature template, and using second file for adding absent feature
#first file is original full set feature file
#second file is for zero selected feature

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


def combine_feature_2file(file1,file2):
	#df1=pd.read_csv(file1, header=None)
	#df2=pd.read_csv(file2, header=None)

	df1=pd.read_csv(file1)
	df2=pd.read_csv(file2)

	df1 = df1.iloc[:, 1:] #remove first column (gene label)
	df2=df2.multiply(0)

	label1=df1.iloc[:,-1] #copy row of label in the last column

	df1 = df1.iloc[:, :-1] #remove last column
	df2 = df2.iloc[:, :-1]
	df2 = df2.head(len(df1.index))
	#print df2
	#print 'File1 total row',len(df1.index)

	df_combine= df1.combine_first(df2)
	df_combine=df_combine.fillna(0)

	df_combine=pd.concat([df_combine, label1], axis=1)

	#print df_combine
	#df_combine.sort_index(axis=0, inplace=True) #re-ordering column
	#print df_combine
	df_combine.to_csv('./comb.temp', index=False, float_format='%.0f')
	
combine_feature_2file(sys.argv[1],sys.argv[2])
#check_value_correct(sys.argv[1],sys.argv[2],'comb.temp')
#python add_absent_of_selected_feature.py test_combind_file1.csv test_combind_file2.csv
