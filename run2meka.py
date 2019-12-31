import commands
import sys,os
#from select_feature_column import *
from select_feature_panda import select_feature_from_file1_forMEKA
import subprocess
import numpy as np

##History
# -0.1- add extract_label_with_prob_MEKA method


def extract_label_with_prob_MEKA(Type,MEKA_input,MEKA_raw_output,output,threshold):
	#usage
	#extract_label_with_prob_MEKA('allGO_multilabel_uniq_ready_FS.arff','BPNN_RF_5.mdout','MEKA_score.temp',0.5)
	print('extract probability score ...')
	
	#|==== PREDICTIONS (N=55.0) =====>
	#|    1 [ 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ] [ 1,0 0,4 0,6 0,3 0,0 1,0 0,0 0,4 0,4 0,0 0,0 0,0 0,3 0,0 0,3 0,0 0,3 0,3 1,0 0,0 0,4 0,3 0,4 0,4 ]
	#|    2 [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ] [ 0,0 0,2 0,0 0,0 0,0 0,1 0,0 0,2 0,2 0,1 0,0 0,0 0,0 0,1 1,0 0,1 0,1 0,0 0,1 0,1 0,2 0,1 0,2 0,2 ]
	#|==============================<
	w2file=open(output,'w')
	GO_term=[]
	for i in open(MEKA_input):
		i=i.rstrip()
		if '@attribute ' in i:
			if 'GO:' in i:
				col=i.split()
				GO_term.append(col[1])
	count=0
	for i in open(MEKA_raw_output):	
		if i != '':
			i=i.rstrip()
			#print(i)

			if '|==============================<' in i:
				count=0

			if count >= 1:
				if Type == 'Predicted':
					col=i.split('[')[2].split(']')[0]
					col=col.replace(',','.')
					#print(col)
					score=col.split()
					for j in range(0,len(GO_term)):
						if float(score[j]) >= threshold:
							w2file.write('protein_'+str(count)+'\t'+GO_term[j]+'\t'+score[j]+'\n')

					count+=1
				elif Type == 'Curated':
					col=i.split('[')[1].split(']')[0]
					#col=col.replace(',','.')
					#print(col)
					score=col.split()
					for j in range(0,len(GO_term)):
						if float(score[j]) == 1:
							w2file.write('protein_'+str(count)+'\t'+GO_term[j]+'\t'+score[j]+'\n')

					count+=1			
				else:
					print('Specify output type !!!')

			if '|==== PREDICTIONS' in i:
				count+=1
	print('Done!')
	w2file.close()

def convert_arff2csv(arff,csv):
	input_file=arff
	output_file=csv
	
	current_path=commands.getoutput('pwd')

	print 'Read input_file',input_file
	#input file is final_combination3.csv
	if os.path.isfile(input_file) == True:
		#use this flow to convert file to nominal and saving as arff
		file_kf=open('./convert2csv.kf','r')
		file_kf_read=file_kf.read()
		file_kf_read=file_kf_read.replace("csv_inputfile",current_path+'/'+input_file)
		file_kf_read=file_kf_read.replace("outputpath",current_path)
		file_kf_read=file_kf_read.replace("outputfilename",output_file)
		file_kf.close()

		file_kf_new=open('convert2csvRUN.kf','w')
		file_kf_new.write(file_kf_read)
		file_kf_new.close()
		num_c=1
		while True: #to make sure that converted_input.arff is created
			print 'Check output file , round ',num_c

			cmd='java -cp ./weka-3-9-1/weka.jar -Xmx10000m -Xss5000m weka.knowledgeflow.FlowRunner convert2csvRUN.kf'
			ps_con1=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
			print(ps_con1.stdout.read())

			num_c+=1
			if int(os.path.getsize(output_file)) > 0 :
				print 'Pass!! : outputfile is created'
				break

		os.system('rm convert2csvRUN.kf')

	else:
		print 'No input file =>  ',input_file

#select feature using first file as filter to select from the second file
#If some feature are absent in the second file, the zero column will be added.
  
def check_value_correct(file1,file2,file_combind): #check value in whole metrix
	print 'Start checking values in combined file ...'
	sum_file1=0
	for i in open(file1):
		if '#' not in i and 'e' not in i:
			i=i.replace('\n','')
			line2=i.split(',')
			del line2[0]
			del line2[len(line2)-1]
			
			sum_line2=0
			for j in line2:
				sum_line2+=int(j)

			sum_file1+=sum_line2

	sum_file2=0
	for i in open(file2):
		if '#' not in i and 'e' not in i:
			i=i.replace('\n','')
			line2=i.split(',')
			del line2[0]
			del line2[len(line2)-1]
			
			sum_line2=0
			for j in line2:
				if 'photo' not in j:
					sum_line2+=int(j)

			sum_file2+=sum_line2	

	sum_file3=0
	for i in open(file_combind):
		if '#' not in i and 'e' not in i:
			i=i.replace('\n','')
			line3=i.split(',')

			del line3[0]
			del line3[len(line3)-1]
			sum_line3=0
			for j in line3:
				sum_line3+=int(j)

			sum_file3+=sum_line3

	print "### Checking value (This step no need to be equal)###"
	print 'File 1 => ',sum_file1
	print 'File 2 => ',sum_file2
	print 'File 1 + File2 ',sum_file1+sum_file2
	print 'File combined => ',sum_file3


def convert_csv2arff(arff,csv):
	input_file=csv
	output_file=arff
	
	current_path=commands.getoutput('pwd')

	print 'Read input_file',input_file
	#input file is final_combination3.csv
	if os.path.isfile(input_file) == True:
		#use this flow to convert file to nominal and saving as arff
		file_kf=open('./convert2arff.kf','r')
		file_kf_read=file_kf.read()
		file_kf_read=file_kf_read.replace("csv_inputfile",current_path+'/'+input_file)
		file_kf_read=file_kf_read.replace("outputpath",current_path)
		file_kf_read=file_kf_read.replace("outputfilename",output_file)
		file_kf.close()

		file_kf_new=open('convert2arffRUN.kf','w')
		file_kf_new.write(file_kf_read)
		file_kf_new.close()
		num_c=1
		while True: #to make sure that converted_input.arff is created
			print 'Check output file , round ',num_c

			cmd='java -cp ./weka-3-9-1/weka.jar -Xmx10000m -Xss5000m weka.knowledgeflow.FlowRunner convert2arffRUN.kf'
			ps_con2=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
			print(ps_con2.stdout.read())

			num_c+=1
			if int(os.path.getsize(output_file)) > 0 :
				print 'Pass!! : outputfile is created'
				break

		os.system('rm convert2arffRUN.kf')
	else:
		print 'No input file =>  ',input_file


def replace_specific_header(modified_file,ref_file,numGOline,out_file): #modified_file ,ref_file = file for to be used to replace the modified file
	ref_head=commands.getoutput('head -n'+str(numGOline)+' '+ref_file)
	modified_tail=commands.getoutput('tail -n+'+str(numGOline+1)+' '+modified_file)

	w=open(out_file,'w')
	w.write(ref_head+'\n'+modified_tail)
	w.close()

def convert_meka_output(filename, label_n_col, outfile):
	#to convert binary output from arff file to GO term
	col=commands.getoutput('grep "@" -v '+filename+'| grep -v "^$" |cut -d"," -f 1-'+str(label_n_col)+' ')
	head=commands.getoutput('grep "@attribute" '+filename+'| head -n '+str(label_n_col)+' ')
	head=head.split('\n') 
	#head[1] is labeled GO term
	
	col=col.split('\n')
	list_GO='' #all GO terms

	for i in col:
		
		i=i.split(',')
		n=0 #column number
		No_GO=True
		for j in i:
			if j == '1':
				GO=head[n].split()
				GO[1]=GO[1].replace('_binarized','')
				list_GO+=GO[1]+'; '
				No_GO=False
			n+=1
		if No_GO==False:
			list_GO=list_GO[:-2]
			list_GO+='\n' #all GO terms	
		else:
			list_GO+='No_GO\n'
	list_GO=list_GO[:-1]
	w=open(outfile,'w')
	w.write(list_GO)
	w.close()

def count_go_predicted_with_prob(GO_file,outputGOfile):  ## -0.1- ## use GO_file =predicted_set_predicted.tsv

	#GO file#
	#protein_1	GO:0009579	1.0
	#protein_1	GO:0030089	1.0
	#protein_2	GO:0009579	1.0
	#protein_2	GO:0030089	1.0

	cmd='cut -f2 '+GO_file+'| sort -u'
	ps_grep=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	out_meka=ps_grep.stdout.read()

	count_GO=out_meka.split('\n')

	cmd2='cut -f2,3 '+GO_file
	ps_grep2=subprocess.Popen(cmd2,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	out_meka2=ps_grep2.stdout.read().split('\n') #GO:0030089	1.0

	w2file=open(outputGOfile,'w')
	no_go_found=0
	allgo_found=0

	for GOuniq in count_GO:
		collected_GO_prob=[]
		if GOuniq != '':
			if GOuniq == 'No_GO':
				no_go_found+=1

			else:
				allgo_found+=1	
				for search_GO in out_meka2:
					search_GO=search_GO.split('\t')
					if GOuniq == search_GO[0]:#GO:0030089	1.0
						collected_GO_prob.append(float(search_GO[1]))
				
				average_prob = np.mean(collected_GO_prob)
				w2file.write(GOuniq+'\t'+str(average_prob)+'\n')

	w2file.close()
	if allgo_found >= 1 :
		if no_go_found == allgo_found :
			return 'Photo-GO not found'
		else:
			return 'Photo-GO found'

def count_go_predicted(GO_file,outputGOfile):  ## -0.0- ##

	#GO file#
	#GO:0009579; GO:0030089; GO:0034357
	#GO:0009579; GO:0030089; GO:0034357
	#GO:0009579; GO:0030089; GO:0034357

	outfile='listGO.temp'
	out=open(outfile,'w')
	for line in open(GO_file):
		line=line.replace('\n','')
		col=line.split(';')
		for j in col:
			out.write(j.replace(' ','')+'\n')
	out.close()
	count_GO=commands.getoutput('sort '+outfile+'|uniq -c ')
	count_GO=count_GO.split('\n')

	total_prediction_instance=int(commands.getoutput('grep "" -c '+GO_file))

	w2file=open(outputGOfile,'w')
	no_go_found=0
	allgo_found=0
	for GOuniq in count_GO:
		if GOuniq != '':
			GOuniq=GOuniq.split()
			percentGO = (float(GOuniq[0])/total_prediction_instance)*100
			w2file.write(GOuniq[1]+'\t'+str(percentGO)+'\n')
			if GOuniq[1] == 'No_GO':
				no_go_found+=1

			allgo_found+=1	

	w2file.close()
	if allgo_found >= 1 :
		if no_go_found == allgo_found :
			return 'Photo-GO not found'
		else:
			return 'Photo-GO found'
	

def run_MEGO(file_label,csv_for_selected_feature,arff_for_selected_feature,input_from_weka):
	print 'Start running MEKA prediction ...'
	#run these steps
	#requirement file : allGO_multilabel_uniq.csv (trained file), query_final.arff (generated file from weka protocol), allGO_multilabel_uniq_sparse.arff (trained file in arff)
	print commands.getoutput('sed s/" \' "/" "/g '+input_from_weka+' |sed s/"\' "/" "/g |sed s/"\'"/""/g > inputfrom_weka.arff')
	convert_arff2csv('inputfrom_weka.arff','query_final_meka.csv')	#[1] arff , [2] csv


	#need to change file input
	select_feature_from_file1_forMEKA(csv_for_selected_feature,'query_final_meka.csv','selected_feature_final_meka.csv') #file1 => file contain selected header(allGO_multilabel_uniq.csv), file2 => unselected header file, [3] output

	check_value_correct(csv_for_selected_feature,'query_final_meka.csv','selected_feature_final_meka.csv') #[1] ref file , [2] unselected header file, result_file from previous step , [3] output

	convert_csv2arff('selected_feature_final_meka.arff','selected_feature_final_meka.csv') #[1] arff , [2] csv

	#need to change file input
	replace_specific_header('selected_feature_final_meka.arff',arff_for_selected_feature,27,'selected_feature_head.arff')#[1] modified_file =file to be replaced , [2] ref_file = file for to be used to replace the modified file, [3] number = row for replace, [4] final arff file

	#for MEKA prediction
	cmd='java -Xmx10000m -Xss10000m -cp "./meka-release-1.9.2-SNAPSHOT/lib/*" meka.classifiers.multilabel.MULAN -t '+arff_for_selected_feature+' -T selected_feature_head.arff -predictions prediction.out.arff -l ./update_MEKAmodel/MEKA_current.model -verbosity 6 -threshold 0.1 > raw_MEKA_predicted.out'
	ps_meka=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	out_meka=ps_meka.stdout.read()
	print(out_meka)
	#print commands.getoutput()

	print 'convert_meka_output ... '
	#convert_meka_output('prediction.out.arff', 24, 'MEGO_predicted.temp') ## -0.1- ##
	extract_label_with_prob_MEKA('Predicted','./update_MEKAmodel/allGO_multilabel_uniq_ready_FS_u.arff','raw_MEKA_predicted.out','predicted_set_predicted.tsv',0.1) ## -0.1- ##

	print 'done '

	print 'count_go_predicted ... '

	cmd='grep "" -c predicted_set_predicted.tsv' ## -0.1- ##
	ps_pred=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	out_pred=ps_pred.stdout.read()
	print(out_pred)

	if int(out_pred) <= 0:
		print 'No GO prediction result!!'
		found_go='No GO prediction result!!'
	else:
		#found_go=count_go_predicted('MEGO_predicted.temp',file_label+'_MEKA_go.go')
		found_go=count_go_predicted_with_prob('predicted_set_predicted.tsv',file_label+'_MEKA_go.go') ## -0.1- ##
	print 'done '

	commands.getoutput('rm prediction.out.arff predicted_set_predicted.tsv raw_MEKA_predicted.out selected_feature_final_meka.arff selected_feature_head.arff selected_feature_final_meka.csv query_final_meka.csv inputfrom_weka.arff')

	return found_go
	#required file for final version
	#1 181GO_label_181nonphoto.csv
	#2 final_18GO_label_181nonphoto.csv.arff

#run_MEGO('label_','181GO_label_181nonphoto_removeuseless.csv','final_18GO_label_181nonphoto.csv.arff','query_final.arff')
