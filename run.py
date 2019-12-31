import commands,os,sys
from conclude_prediction import *
from run2meka import *
import subprocess

path_to_blast='' #already installed
blast_eval=sys.argv[2]
blast_core=sys.argv[3]
ram_usage=str(int(sys.argv[4])*1000)
blast_hit=sys.argv[5]

print('PhotoMod2 v 1.0')
print('Photosynthetic protein classifier')
print('=====================================')
print('- Usage: python run.py [test_protein.fasta] [blast evalue] [cpu] [ram(GB)] [blast hit]')
print('- Output: (1) prediction_report.tab: containing photosynthetic protein classification result , (2) GO_prediction.tab: containing photosynthesis-related function prediction result')
print('=====================================')
print('============= Update ================')
print('- add PhotoMod2 multi-label model according to the second manuscript')
print('- update PhotoMod1 binary model according to the first manuscript')
print('=====================================')

def dict_24_GO_terms(csv_GO):
	dict_24GO={}
	for i in open(csv_GO):
		i=i.rstrip()
		j=i.split(',')
		dict_24GO[j[0]]=j[1]
	return dict_24GO

def check_header_feature(file1,file2): #chech identity of AFff header by sampling features
	row_f1=commands.getoutput('tail -n+3 '+file1+' ').split('\n')	
	row_f2=commands.getoutput('tail -n+3 '+file2+' ').split('\n')	

	if (row_f1[0] == row_f2[0]) and (row_f1[100] == row_f2[100]) and (row_f1[1450] == row_f2[1450]) and (row_f1[7000] == row_f2[7000]) and (row_f1[10000] == row_f2[10000]) and (row_f1[18000] == row_f2[18000]) :
		return True
	else:
		return False

def select_only_new_query(input_file,index_prot_file): #to separate between identified protein (have been collected) and the new proteins
	file_index=open(index_prot_file,'r')
	read_index=file_index.read()
	file_index.close()

	read_index=read_index.split('\n')
	file_match=open('matched_prot_file.query','w')
	file_unmatch=open('unmatched_prot_file.query','w')
	for i in open(input_file):
		i=i.replace('\n','')
		j=i.split()
		try:
			prot_match_index=read_index.index(j[1]) 
			file_match.write(i+'\t'+str(prot_match_index)+'\n')
		except ValueError:
			file_unmatch.write(i+'\n')
		
	file_match.close()
	file_unmatch.close()

def select_identified_query(index_prot_file,feature_file,header_file): #to separate between identified protein (have been collected) and the new proteins
	file_index=open(index_prot_file,'r') #index_prot_file is [query \t index] file 
	read_index=file_index.read()
	file_index.close()
	read_index=read_index.split('\n')

	#make list of file that contain feature
	feature_list=open(feature_file,'r')
	read_feature=feature_list.read()
	feature_list.close()
	read_feature=read_feature.split('\n')
	
	file_match=open('matched_prot_query.arff','w')

	for i in read_index:
		if i != '':
			i=i.split() #i[0] => query name , i[1] => index

			file_match.write(read_feature[int(i[8])]+'\n')
	
	file_match.close()
	commands.getoutput('cat '+header_file+' matched_prot_query.arff > matched_query.arff')

def blast_query_to_MyDB(list_query_fasta,DB_name,evalue_cutoff):#DB_name our photosynthetic genes in blast DB format, max_seq =max_target_seqs
	#blastp using query gene file
	
	blastp=path_to_blast+"blastp -db "+DB_name+" -query "+list_query_fasta+" -out eval_matched2DB.query -max_hsps 1 -max_target_seqs "+blast_hit+" -outfmt '6 qseqid sseqid evalue pident qstart qend sstart send' -num_threads "+str(blast_core)

	print "#blastp using query gene file"
	print blastp
	os.system(blastp)
	print commands.getoutput("awk '$3<="+str(evalue_cutoff)+"' eval_matched2DB.query > matched2DB.query") #filter by evalue
	print commands.getoutput("rm eval_matched2DB.query")
	print "blastp finished!!!"

inputfile=open(sys.argv[1],'r')
read_input=inputfile.read()
read_input=read_input.split('>')
inputfile.close()
dict_24GO = dict_24_GO_terms('./update_MEKAmodel/24selected_GO_Photosynthesis.csv')
collected_MEGO_out='NO.\tName\tGO\tDescription\tProbability\tMethod\n'
binaryout=open('prediction_report.tab','w')
binaryout.write('Query	totalBlastHits	PhotoHits(%)	nonphotoHits(%)	noPredictHits	prediction	Ave_Prob_prediction\n')
binaryout.close()
#############################################################################################
count=1
for seq in read_input:
	if seq != '':
		wfasta=open('1_query_file.fasta','w')
		wfasta.write('>'+seq)
		wfasta.close()
		print '>'+seq
		prot_name=seq.rstrip().split('\n')[0].split()[0]

		#1. Blast target proteins to our proteome 
		path2db = './OnlyGenome/all_154genomes.fasta'
		blast_query_to_MyDB("1_query_file.fasta",path2db,blast_eval) #change to evalue cutoff criteria

		#in the case of new genes
		if int(commands.getoutput("grep '' -c matched2DB.query")) > 0:
			

			path_DB='./collected_NB_feature/'
			select_only_new_query('matched2DB.query',path_DB+'index_prot.nb') #separate already identified prot and new prot

			#filter the new protein
			#######################################
			#get NBdata from DB if mached protein is found, or use the new prot profile to be a final query file if no matched prot
			print 'matched_prot_file.query => ',int(commands.getoutput('grep "" -c matched_prot_file.query' ))
			print 'unmatched => ', int(commands.getoutput("grep '' -c unmatched_prot_file.query"))

			if int(commands.getoutput("grep '' -c unmatched_prot_file.query")) > 0 :#check wether there are matched prot found

				os.system("python export_GO.py 24240 cds_MCL_cluster_evalue10")
				os.system("python make_metrix_conserve_makeprediction.py 10 './call_neighbor_3evalue/e10.C.uniq_nb_id.id' no_label 0.41079 2.61799")

				os.system("rm cluster_export.txt GO_neighbor_enrichment.db *.out *.out2")
				os.system("rm -r all_neighboring conservation_score")

				os.system("python export_GO.py 39068 cds_MCL_cluster_evalue50")
				os.system("python make_metrix_conserve_makeprediction.py 50 './call_neighbor_3evalue/e50.C.uniq_nb_id.id' no_label 0.41079 2.61799")

				os.system("rm cluster_export.txt GO_neighbor_enrichment.db *.out *.out2")
				os.system("rm -r all_neighboring conservation_score")

				os.system("python export_GO.py 39896 cds_MCL_cluster_evalue100")
				os.system("python make_metrix_conserve_makeprediction.py 100 './call_neighbor_3evalue/e100.C.uniq_nb_id.id' photo 0.41079 2.61799")

				os.system("rm cluster_export.txt GO_neighbor_enrichment.db *.out *.out2")
				os.system("rm -r all_neighboring conservation_score")

				os.system("python concat_feature.py metrix_cv_nb_predicted_e10.csv metrix_cv_nb_predicted_e50.csv no_label >metrix_cv_nb_e10_50.csv")
				os.system("python check_file_combind.py metrix_cv_nb_predicted_e10.csv metrix_cv_nb_predicted_e50.csv metrix_cv_nb_e10_50.csv")

				os.system("python concat_feature.py metrix_cv_nb_e10_50.csv metrix_cv_nb_predicted_e100.csv photo >metrix_cv_nb_3evalue.csv")
				os.system("python check_file_combind.py metrix_cv_nb_e10_50.csv metrix_cv_nb_predicted_e100.csv metrix_cv_nb_3evalue.csv")

				if int(commands.getoutput('grep "" -c metrix_cv_nb_3evalue.csv' )) >= 2: #to check whether the file contain neighbor profile
					os.system("python add_absent_of_selected_feature.py metrix_cv_nb_3evalue.csv comb.csv")

					#os.system("cp comb.temp "+str(count)+"_collect2DB.csv") #to copy for collectig into DB
			
					#in the case of already identified genes
					os.system("python file_convertor.py comb.temp 2_query")
					os.system("rm converted_input.arff")

					
					#check file identity between the new file and in DB file
					print commands.getoutput('sed s/" \' "/" "/g 2_query_converted_input.arff |sed s/"\' "/" "/g |sed s/"\'"/""/g > 3_query_converted_input.arff')#remove quate
					data_check = check_header_feature('3_query_converted_input.arff' , path_DB+'head_arff.nb')
					print 'Header checking ... ', data_check

					#add data to the DB
					if data_check==True:
						print commands.getoutput('tail -n+18856 3_query_converted_input.arff > new_arff_data.arff')
						print commands.getoutput('mv '+path_DB+'body_feature.nb '+path_DB+'body_feature.nb.previous')
						print commands.getoutput('cat '+path_DB+'body_feature.nb.previous new_arff_data.arff > '+path_DB+'body_feature.nb')
						#update index proteins
						print commands.getoutput('mv '+path_DB+'index_prot.nb '+path_DB+'index_prot.nb.previous')
						print commands.getoutput('cut -f1 -d"," metrix_cv_nb_3evalue.csv | tail -n+2 > updated_proteinindex.txt')
						print commands.getoutput('cat '+path_DB+'index_prot.nb.previous updated_proteinindex.txt > '+path_DB+'index_prot.nb')
						print commands.getoutput('rm updated_proteinindex.txt ')
					


			if int(commands.getoutput('grep "" -c matched_prot_file.query' )) > 0 :#check wether there are matched prot found

				select_identified_query('matched_prot_file.query',path_DB+'body_feature.nb',path_DB+'head_arff.nb')
				
				if int(commands.getoutput("grep '' -c unmatched_prot_file.query")) > 0: #to check whether the file contain neighbor profile
					if os.path.exists('./metrix_cv_nb_3evalue.csv') == True :
						if int(commands.getoutput('grep "" -c metrix_cv_nb_3evalue.csv' )) >= 2 :

							print commands.getoutput('cat matched_query.arff new_arff_data.arff > query_final.arff')

				else:
					print commands.getoutput('mv matched_query.arff query_final.arff')
			
			elif int(commands.getoutput("grep '' -c unmatched_prot_file.query")) > 0: #to check whether the file contain neighbor profile
					if int(commands.getoutput('grep "" -c metrix_cv_nb_3evalue.csv' )) >= 2 :
						print commands.getoutput('cat '+path_DB+'head_arff.nb new_arff_data.arff > query_final.arff')

			else:
				print ' NO both matched protein and unmatched proteins..... !!!!! '
			#prediction
			print('Run PhotoMod1...')
			cmd="java -cp ./weka-3-9-1/weka.jar -Xmx"+str(ram_usage)+"m -Xss10000m weka.Run weka.classifiers.meta.FilteredClassifier -T query_final.arff -o -l current_model.model -classifications weka.classifiers.evaluation.output.prediction.CSV > prediction.csv"
			ps_weka=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
			out_weka=ps_weka.stdout.read()

			conclude_result()

			########################### Run MEKA model ###################
			print('Run PhotoMod2...')
			status_MEGO=run_MEGO('query_'+str(count),'./update_MEKAmodel/allGO_multilabel_uniq_ready_FS_u.arff.csv','./update_MEKAmodel/allGO_multilabel_uniq_ready_FS_u.arff','query_final.arff') #label file as input

			if status_MEGO=='Photo-GO found':
				print 'Photo-GO found'
				for h in open('query_'+str(count)+'_MEKA_go.go'):
					h=h.rstrip()
					if h != '':
						h=h.split('\t')
						collected_MEGO_out+='query_'+str(count)+'\t'+prot_name+'\t'+h[0]+'\t'+dict_24GO[h[0]]+'\t'+h[1]+'\t'+'\tphotomod2\n'


				os.system('rm query_'+str(count)+'_MEKA_go.go')

			else:
				collected_MEGO_out+='query_'+str(count)+'\t'+prot_name+'\tNo_GO_predicted\t0.0\tphotomod\n'
				

			#os.system("rm prediction.csv")
			if os.path.exists("metrix_cv_nb_3evalue.csv") == True:
				os.system("mv metrix_cv_nb_3evalue.csv metric_gene"+str(count)+".csv")
				os.system("rm metrix_cv_nb_e10_50.csv")
				os.system("rm comb.temp matched_query.arff new_arff_data.arff")
				os.system("rm 2_query_converted_input.arff 3_query_converted_input.arff ")

		else:
			print 'No blast matched!!'
		os.system("rm matched2DB.query")

		os.system("rm matched_prot_file.query")
		os.system("rm unmatched_prot_file.query")
		count+=1

wfile=open('GO_prediction.tab','w')
wfile.write(collected_MEGO_out)
wfile.close()	
		
#need to change number of feature if the number of feature changed at line "tail -n+18856 2_query_converted_input.arff > new_arff_data.arff" 
