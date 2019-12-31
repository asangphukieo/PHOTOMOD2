import sys
import commands
import pandas as pd
import numpy as np
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


def select_feature_from_file1(file1,file2,out_file):
	print 'Select feature column ...'
	#df1=pd.read_csv(file1, header=None)
	#df2=pd.read_csv(file2, header=None)

	df1=pd.read_csv(file1)
	df2=pd.read_csv(file2, low_memory=False)

	#df1 = df1.iloc[:, 1:] #remove first column (gene label)
	#df1=df1.multiply(np.nan)
	df1=df1[:0] #remove every row
	nan_table=pd.DataFrame(np.nan, index = np.arange(len(df2.index)), columns = df1.columns) #make the new table which the number of column equals to df1 columns and the number of rows equals to df2 row
 
	df1=df1.append(nan_table) #append the nan_table to  df1 table

	#label1=df1.iloc[:,-1] #copy row of label in the last column
	#df1 = df1.iloc[:, :-1] #remove last column
	#df_combine= df1.combine_first(df2) #combine features of 2 file but using the feature order and value from df1 and adding the new distinct feature from df2 and their value 

	df2_filtered=df1.fillna(df2) #replace NaN value with the vaue from df2
	df2_filtered=df2_filtered.fillna(0) #the rest of NaNs are replaced with 0

	#df_combine=pd.concat([df_combine, label1], axis=1)
	#df_combine.sort_index(axis=0, inplace=True) #re-ordering column

	#move 18 GO labels to the front
	GOs=df2_filtered.iloc[:,-18:] #copy last 18 columns of GO labels
	df2_filtered=df2_filtered.iloc[:, :-18] #remove last 18 columns of GO labels

	df2_GOs=pd.concat([GOs, df2_filtered], axis=1)
	
	df2_GOs.to_csv(out_file, index=False, float_format='%.0f')

	head_file3=commands.getoutput('head -n1 '+out_file)

	header3_list=(head_file3.replace(' ','')).split(',')
	print 'File combinded contains ',len(header3_list)-1,'features'

def select_feature_from_file1_forMEKA(file1,file2,out_file):
	print 'Select feature column ...'
	#df1=pd.read_csv(file1, header=None)
	#df2=pd.read_csv(file2, header=None)

	df1=pd.read_csv(file1)
	df2=pd.read_csv(file2, low_memory=False)

	#df1 = df1.iloc[:, 1:] #remove first column (gene label)
	#df1=df1.multiply(np.nan)
	df1=df1[:0] #remove every row
	nan_table=pd.DataFrame(np.nan, index = np.arange(len(df2.index)), columns = df1.columns) #make the new table which the number of column equals to df1 columns and the number of rows equals to df2 row
 
	df1=df1.append(nan_table) #append the nan_table to  df1 table

	#label1=df1.iloc[:,-1] #copy row of label in the last column
	#df1 = df1.iloc[:, :-1] #remove last column
	#df_combine= df1.combine_first(df2) #combine features of 2 file but using the feature order and value from df1 and adding the new distinct feature from df2 and their value 

	df2_filtered=df1.fillna(df2) #replace NaN value with the vaue from df2
	df2_filtered=df2_filtered.fillna(0) #the rest of NaNs are replaced with 0

	#df_combine=pd.concat([df_combine, label1], axis=1)
	#df_combine.sort_index(axis=0, inplace=True) #re-ordering column

	#move 18 GO labels to the front
	#GOs=df2_filtered.iloc[:,-18:] #copy last 18 columns of GO labels
	#df2_filtered=df2_filtered.iloc[:, :-18] #remove last 18 columns of GO labels

	#df2_GOs=pd.concat([GOs, df2_filtered], axis=1)
	
	df2_filtered.to_csv(out_file, index=False, float_format='%.0f')

	head_file3=commands.getoutput('head -n1 '+out_file)

	header3_list=(head_file3.replace(' ','')).split(',')
	print 'File combinded contains ',len(header3_list)-1,'features'

#combine_feature_2file(sys.argv[1],sys.argv[2],sys.argv[3])
#select_feature_from_file1('181GO_label_181nonphoto_removeuseless.csv','query_final_meka.csv','selected_feature_final_meka.csv') 
#check_value_correct(sys.argv[1],sys.argv[2],'comb.temp')
#python add_absent_of_selected_feature.py test_combind_file1.csv test_combind_file2.csv
