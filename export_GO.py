import commands
import sqlite3
import os,sys
import scipy.stats as ss
import numpy
import itertools
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from build_tree import *

#input_blast_file='blast_photo_cutoff.out'
cluster_conserverd=sys.argv[1] #number of cluster that conserved more than 2 members
DB_mcl_tablename=sys.argv[2] #'cds_MCL_cluster_evalue100'

#grep '' blast_photo_cutoff.out |cut -f2|cut -f2 -d'|'|wc -l

def call_cluster_id(gen_query):
	conn = sqlite3.connect("Photosyn.db")
	c = conn.cursor()
	#convert gene query to cluster id
	c.execute("SELECT cluster_id FROM "+DB_mcl_tablename+" WHERE seq_id = '"+gen_query+"'")
	result1 = c.fetchall()
	conn.close()

	if len(result1) != 0:
		return result1[0][0],gen_query
	else:
		return "Not cluster found: unique gene"

def multiple_call_cluster_id(file_gen_query,output):
	seqs=commands.getoutput("cut -f1 "+file_gen_query)
	seqs=seqs.split('\n')
	w=open(output,'w')
	for ind in seqs:
		out=call_cluster_id(ind)
		if out != 'Not cluster found: unique gene':
			#print str(out[0])+'\t'+out[1]
			w.write(str(out[0])+'\t'+out[1]+'\n')
	w.close()

def call_cluster_id_from_neighor_file(file_gen_query,clust_id): #file_gen_query=neighboring file name , clust_id = output file label
	list_prot=commands.getoutput("awk '! /#/' "+file_gen_query+"|cut -f1,2")
	line_list=list_prot.split('\n')

	if len(line_list) >1:#To discard the file that no neighboring gene found
		conn = sqlite3.connect("Photosyn.db")
		c = conn.cursor()
		out_write=open('cluster_'+str(clust_id)+'.nout','w')

		for prot_query in line_list:
			prot_query=prot_query.split()
			query_gene=prot_query[0]
			neighbor_gene=prot_query[1]	
			c.execute("SELECT cluster_id FROM "+DB_mcl_tablename+" WHERE seq_id = '"+neighbor_gene+"' ")#Cluster id higher than cluster_conserverd is unique cluster
			result = c.fetchall()
			if len(result) == 0:
				#To confirm whether proteins exist in our set
				c.execute("SELECT cds_ID FROM cds WHERE cds_ID = '"+neighbor_gene+"'")
				result2 = c.fetchall()

				if len(result2) == 0:
					print '#Unknown protein (not present in our data set) of cluster ',i
				else:
					out_write.write('uniq'+'\t'+str(neighbor_gene)+'\t'+str(query_gene)+'\n') #uniq mean that protein might no have any relationship with our set of protein and are filtered out in the blast all vs all step
			else:	
				#Neighboring_cluster_id Neighboring_gene_name Query_gene_name
				out_write.write(str(result[0][0])+'\t'+str(neighbor_gene)+'\t'+str(query_gene)+'\n')
		out_write.close()
		conn.close()
		return True
	else:
		return False
		

def neighbor_GO(cluster_id):
	conn = sqlite3.connect("Photosyn.db")
	c = conn.cursor()
	#extract GO term that are enriched in such cluster
	cluster_query=str(cluster_id)
	c.execute("SELECT GO_id FROM GO_enrichment WHERE cluster_id = "+cluster_query+" AND pvalue_adjust <= 0.05")
	result2 = c.fetchall()
	neighbor_members=len(result2)
	#print 'cluster_id',cluster_id,'found',len(result2),'GO_terms'
	GO_list=[]
	for i in result2:
		line=''
		for j in i:		
			GO_list.append(j)

	conn.close()	
	return GO_list,neighbor_members

def neighbor_GO_random(cluster_id): #random GO term with the random times is the same with the number of GO terms of each cluster_id  
	conn = sqlite3.connect("Photosyn.db")
	c = conn.cursor()
	#extract GO term that are enriched in such cluster
	cluster_query=str(cluster_id)
	c.execute("SELECT GO_id FROM GO_enrichment WHERE cluster_id = "+cluster_query+" AND pvalue_adjust <= 0.05")
	result1 = c.fetchall()
	neighbor_members=len(result1)
	
	c.execute("SELECT GO_id FROM GO_enrichment WHERE pvalue_adjust <= 0.05 ORDER BY RANDOM() LIMIT "+str(neighbor_members))
	result2=c.fetchall()
	#print 'Random set contains ',len(result2)
	#print result2

	GO_list=[]
	for i in result2:
		line=''
		for j in i:		
			GO_list.append(j)

	conn.close()
	return GO_list,neighbor_members

def call_neighboring_gene(list_cluster_id):
	conn = sqlite3.connect("Photosyn.db")
	c = conn.cursor()
	path_to_neighboring='./all_neighboring/'
	cluster_id=commands.getoutput("cut -f1 "+list_cluster_id)
	cluster_id=cluster_id.split('\n')
	#print cluster_id

	file_p=open('neighbor_enrichment.out','w')
	file_p.write("Query_cluster_id\tneighbor_cluster_id\ttotal_neighbor_cluster_id\tx\tM\tK\tN\tPvalue\tadjusted_Pvalue\n")
	file_p.close()	
	for i in cluster_id:
		if int(i) <= cluster_conserverd:# cluster id that more than cluster_conserverd is unique (not conserved)	
			#we analyse only gene that is conserved among these genome
			file_neighboring_cluster='query_prot_id_'+str(i)+'.ext.neighbor'
			N_found=call_cluster_id_from_neighor_file(path_to_neighboring+file_neighboring_cluster,i)
			if N_found == True:
				total_query_neighbor=commands.getoutput('grep "" -c cluster_'+str(i)+'.nout')
				total_neighbor_cluster_id=commands.getoutput('cut -f1 cluster_'+str(i)+'.nout | sort|uniq -c|grep "" -c')
				#print 'total',total_neighbor_cluster_id
				freq_neighbor=commands.getoutput('cut -f1 cluster_'+str(i)+'.nout | sort|uniq -c')
				freq_neighbor=freq_neighbor.split('\n')
				for line in freq_neighbor:
					col=line.split()
					#if len(col)<=1:
						#print 'Found empty line'
					if col[1] != 'uniq':
						freq_cluster = col[0] #while col[1] is cluster_id

						c.execute("SELECT COUNT(seq_id) FROM "+DB_mcl_tablename+" WHERE cluster_id="+col[1])
						count_all_seq_in_cluster = c.fetchall()
						#print count_all_seq_in_cluster[0][0]

					
						enrichment(str(i),col[1],total_neighbor_cluster_id,freq_cluster,604857,count_all_seq_in_cluster[0][0],total_query_neighbor) #604857 is total number of sequences
					else:
						uniq_seq=col[0] #count neighboring gene that is unique sequence
					

				os.system('rm cluster_'+str(i)+'.nout')
	conn.close()

def enrichment(query_cluster_name,neighbor_cluster_id,total_neighbor_cluster_id,x,M,K,N):
	file_p=open('neighbor_enrichment.out','a')
	total_neighbor_cluster_id=int(total_neighbor_cluster_id)
	x=int(x)
	M=int(M)
	K=int(K)
	N=int(N)
	pvalue=1-ss.hypergeom.cdf(x-1,M,K,N)
	adj_pvalue=pvalue*int(total_neighbor_cluster_id) #bonferroni correction

	file_p.write(str(query_cluster_name)+'\t'+str(neighbor_cluster_id)+'\t'+str(total_neighbor_cluster_id)+'\t'+str(x)+'\t'+str(M)+'\t'+str(K)+'\t'+str(N)+'\t'+str(pvalue)+'\t'+str(adj_pvalue)+'\n')
	#print pvalue

	file_p.close()

def write_GO_match_output(file_name,match_method,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11):#matching_method pvalue C_score D_score CD_score phylo_conservation_score precision pre_sd recall re_sd cov cov_not_nan
	#SS=sematic similarity, ANC=ancestor_GO_match, GO=exact_GO_match

	if match_method == 'SS':
		header='#matching_method\tpvalue\tC_score\tD_score\tCD_score\tphylo_conservation_score\tave_GOsemsim\tsd\tcov\tcov_not_nan'
		out_file=open(file_name,'a')
		out_file.write(header+'\n')
		out_file.write(match_method+'\t'+str(col1)+'\t'+str(col2)+'\t'+str(col3)+'\t'+str(col4)+'\t'+str(col5)+'\t'+str(col6)+'\t'+str(col7)+'\t'+str(col8)+'\t'+str(col9)+'\n')
		out_file.close()
	else:
		header='#matching_method\tpvalue\tC_score\tD_score\tCD_score\tphylo_conservation_score\tprecision\tpre_sd\trecall\tre_sd\tcov\tcov_not_nan'
		out_file=open(file_name,'a')
		out_file.write(header+'\n')
		out_file.write(match_method+'\t'+str(col1)+'\t'+str(col2)+'\t'+str(col3)+'\t'+str(col4)+'\t'+str(col5)+'\t'+str(col6)+'\t'+str(col7)+'\t'+str(col8)+'\t'+str(col9)+'\t'+str(col10)+'\t'+str(col11)+'\n')
		out_file.close()
		

	
def call_neighboring_gene_conservation(list_cluster_id): 
	conn = sqlite3.connect("Photosyn.db")
	c = conn.cursor()
	path_to_neighboring='./all_neighboring/'
	cluster_id=commands.getoutput("cut -f1 "+list_cluster_id)
	cluster_id=cluster_id.split('\n')

	for i in cluster_id:
		if int(i) <= cluster_conserverd:# cluster id that more than cluster_conserverd is unique (not conserved)	
			#we analyse only gene that is conserved among these genome
			out_conservation_file=open('CD_score_Qid_'+str(i)+'.out','w')
			#out_conservation_file.write('#query_cluster_id neighbor_cluster_id neighbor_found qeury_found number_of_genome(fix_id)_found C-score D-score CD-score\n')
			out_conservation_file.close()							
			file_neighboring_cluster='query_prot_id_'+str(i)+'.ext.neighbor'
			file_query_cluster='query_prot_id_'+str(i)+'.ext'

			N_found=call_cluster_id_from_neighor_file(path_to_neighboring+file_neighboring_cluster,i)
			Total_query_gene=commands.getoutput("grep '' -c "+path_to_neighboring+file_query_cluster)

			if N_found == True:
				temp_file=open('neighbor_id_'+str(i)+'.temp','w')

				for line_neighbor_ID in open('cluster_'+str(i)+'.nout'):
					line_neighbor_ID=line_neighbor_ID.replace('\n','')
					col_seqID=line_neighbor_ID.split()
					
					if col_seqID[0] != 'uniq':
						#print col_seqID
						c.execute("SELECT Genome_accession.FIX_ID FROM cds,Genome_accession WHERE cds_ID='"+str(col_seqID[1])+"' AND cds.genome_acc=Genome_accession.genome_acc")
						fix_id= c.fetchall()
						#neighboring_cluster_id Fix_id neighboring_gene_name query_gene_name (as output file)
						temp_file.write(str(col_seqID[0])+'\t'+str(fix_id[0][0])+'\t'+str(col_seqID[1])+'\t'+str(col_seqID[2])+'\n')#col_seqID[0] is cluster id of neighbors
				temp_file.close()
				#[query_gene_name neighboring_cluster_id Fix_id ] > sort by column 1 and 2
				os.system("awk '{print $4,$1,$2}' neighbor_id_"+str(i)+".temp|sort -k1,2 -u > neighbor_id_"+str(i)+".temp2")
				os.system("awk -F' ' '{print $2}' neighbor_id_"+str(i)+".temp2| sort | uniq -c > neighbor_id_"+str(i)+".temp3")
				
				
				for freq in open("neighbor_id_"+str(i)+".temp3"):
					freq=freq.split()
					#[frequency neighboring_cluster_ID frequency/total_query_gene]
					#print freq[0],freq[1],float(freq[0])/float(Total_query_gene)#C-score
					conservation_score = float(freq[0])/float(Total_query_gene)	#C-score
					
					#Calculate D-score 	
					get_fixID=commands.getoutput("awk '$2=="+str(freq[1])+" {print $3}' neighbor_id_"+str(i)+".temp2 |sort -u")
					get_fixID=get_fixID.split('\n')
				
					combination='1genome'
					if len(get_fixID) > 1:	#To consider only neighbors that occured more than 1 genomes. If neighbor that occures in 1 genomes (fix id),it will not calculated in D-score, but it will be added to C-score.
						combination=list(itertools.combinations(get_fixID,2))
						
					if combination != '1genome' :
						list_Gdistance=[]
						for fix_g in combination:
							#print fix_g[0],fix_g[1]
							c.execute("SELECT feature_scaling FROM GenomeDistance WHERE (FIX_ID_sp1='"+str(fix_g[0])+"' AND FIX_ID_sp2="+str(fix_g[1])+") OR (FIX_ID_sp1='"+str(fix_g[1])+"' AND FIX_ID_sp2="+str(fix_g[0])+")")
							result_fix_g= c.fetchall()						
						
							list_Gdistance.append(result_fix_g[0][0])
						ave_pairwise_dist = sum(list_Gdistance)/len(list_Gdistance) #D-score
					else:
						ave_pairwise_dist=0.0 #Neighbors are conserved in only 1 genome , so No Distance score 

					CD_score = (conservation_score + ave_pairwise_dist)/2.0 #CD-score
					#CD_score = 2*(conservation_score * ave_pairwise_dist)/(conservation_score + ave_pairwise_dist) #CD-score harmonic mean

					#query_cluster_id neighbor_cluster_id neighbor_found qeury_found number_of_genome(fix_id)_found C-score D-score CD-score					
					out_conservation_file=open('CD_score_Qid_'+str(i)+'.out','a')
					out_conservation_file.write(str(i)+'\t'+str(freq[1])+'\t'+str(freq[0])+'\t'+str(Total_query_gene)+'\t'+str(len(get_fixID))+'\t'+str(conservation_score)+'\t'+str(ave_pairwise_dist)+'\t'+str(CD_score)+'\n' )
					out_conservation_file.close()
				os.system('rm cluster_'+str(i)+'.nout')
				os.system('rm neighbor_id_'+str(i)+'.temp')
				os.system('rm neighbor_id_'+str(i)+'.temp2')
				os.system('rm neighbor_id_'+str(i)+'.temp3')
	conn.close()	

from select_gene_for_blastall import *
def blast_2_ref_genome(query,db):
	
	#os.system("makeblastdb -in "+db+" -parse_seqids -dbtype prot")

	os.system("blastp -db "+db+" -query "+query+" -out "+query+".blastout -max_hsps 1 -outfmt '7 qseqid sseqid evalue pident qstart qend qlen slen length sstart send' -num_threads 1 -max_target_seqs 1")

	#select_gene(query+".blastout",80.0, 80.0, 80.0, 1e-10)

def insert_to_DB(input_file,index):

	conn = sqlite3.connect("GO_neighbor_enrichment.db")
	c = conn.cursor()
	#c.execute("DROP TABLE neighbor_conservation_phylo")
	c.execute("CREATE TABLE neighbor_conservation_phylo (Query_cluster_id INTEGER NOT NULL, neighbor_cluster_id INTEGER NOT NULL, neighbor_found INTEGER NOT NULL, qeury_found INTEGER NOT NULL, number_of_genome_found INTEGER NOT NULL, Check1 Text NOT NULL, Check2 Text NOT NULL, conservation_score DOUBLE NOT NULL , FOREIGN KEY(Query_cluster_id) REFERENCES GO_neighbor_enrichment (Query_cluster_id), FOREIGN KEY(neighbor_cluster_id) REFERENCES GO_neighbor_enrichment (neighbor_cluster_id), PRIMARY KEY (Query_cluster_id, neighbor_cluster_id) )")


	db_list = []
	with open(input_file) as infile:
		
	    for line in infile:
		if '#' not in line:
			data = line.replace('\n','')
			data = data.split()
			#print data
			db_list.append(data)

	c.executemany("INSERT INTO neighbor_conservation_phylo VALUES (?, ?, ?, ?, ?, ?, ?, ?)", db_list)
	c.execute("CREATE INDEX id"+str(index)+" ON neighbor_conservation_phylo(Query_cluster_id, neighbor_cluster_id)")
	conn.commit()
	conn.close()	

def call_neighboring_gene_conservation_phylo(list_cluster_id): 
	conn = sqlite3.connect("Photosyn.db")
	c = conn.cursor()
	path_to_neighboring='./all_neighboring/'
	cluster_id=commands.getoutput("cut -f1 "+list_cluster_id)
	cluster_id=cluster_id.split('\n')

	for i in cluster_id:
		if int(i) <= cluster_conserverd:# cluster id that more than cluster_conserverd is unique (not conserved)	
			#we analyse only gene that is conserved among these genome

			out_conservation_file=open('Conservation_score_Qid_'+str(i)+'.out','w')
			#out_conservation_file.write('#query_cluster_id\tneighbor_cluster_id\tneighbor_gene_found\t qeury_gene_found\tnumber_of_genome(fix_id)_found\tcheck1(genome_pairwise_distance)\tcheck2(totalGenome_found)\ttotal_branch_lenngth(Phylogeny)')
			out_conservation_file.close()
						
			file_neighboring_cluster='query_prot_id_'+str(i)+'.ext.neighbor'
			file_query_cluster='query_prot_id_'+str(i)+'.ext'

			N_found=call_cluster_id_from_neighor_file(path_to_neighboring+file_neighboring_cluster,i) #OUTPUT .nout file
			Total_query_gene=commands.getoutput("grep '' -c "+path_to_neighboring+file_query_cluster)

			if N_found == True:
				temp_file=open('neighbor_id_'+str(i)+'.temp','w')

				for line_neighbor_ID in open('cluster_'+str(i)+'.nout'):
					line_neighbor_ID=line_neighbor_ID.replace('\n','')
					col_seqID=line_neighbor_ID.split()
					
					if col_seqID[0] != 'uniq':
						#print col_seqID
						c.execute("SELECT Genome_accession.FIX_ID FROM cds,Genome_accession WHERE cds_ID='"+str(col_seqID[1])+"' AND cds.genome_acc=Genome_accession.genome_acc")
						fix_id= c.fetchall()
						#neighboring_cluster_id Fix_id neighboring_gene_name query_gene_name (as output file)
						temp_file.write(str(col_seqID[0])+'\t'+str(fix_id[0][0])+'\t'+str(col_seqID[1])+'\t'+str(col_seqID[2])+'\n')#col_seqID[0] is cluster id of neighbors
				temp_file.close()
				#[query_gene_name neighboring_cluster_id Fix_id ] > sort by column 1 and 2
				os.system("awk '{print $4,$1,$2}' neighbor_id_"+str(i)+".temp|sort -k1,2 -u > neighbor_id_"+str(i)+".temp2")
				os.system("awk -F' ' '{print $2}' neighbor_id_"+str(i)+".temp2| sort | uniq -c > neighbor_id_"+str(i)+".temp3")
				
				
				for freq in open("neighbor_id_"+str(i)+".temp3"):
					freq=freq.split()
					#[frequency neighboring_cluster_ID frequency/total_query_gene]
					#print freq[0],freq[1],float(freq[0])/float(Total_query_gene)#C-score
					conservation_score = float(freq[0])/float(Total_query_gene)	#C-score
					
					#Calculate D-score 	
					get_fixID=commands.getoutput("awk '$2=="+str(freq[1])+" {print $3}' neighbor_id_"+str(i)+".temp2 |sort -u")
					get_fixID=get_fixID.split('\n')
				
					combination='1genome'

					check1='No_status' #check the number of pair distance found in build_tree.py
					check2='No_status' #No status means , the process cannot calculate because no value provided
					if len(get_fixID) >= 3:	#To consider only neighbors that occured more than 3 genomes. If neighbor that occures in 1 genomes (fix id),it will not calculated in C score.
						combination=list(itertools.combinations(get_fixID,2))
						#print combination
						
					if combination != '1genome' : #to check wheter gene is conserved
						file_name_pair_dist="pairwise_distance_"+str(i)+".temp"
						list_Gdistance=open(file_name_pair_dist,'w')
						for fix_g in combination:
							#print fix_g[0],fix_g[1]
							c.execute("SELECT orthologous_distance FROM OrthoDistance WHERE (FIX_ID_sp1='"+str(fix_g[0])+"' AND FIX_ID_sp2="+str(fix_g[1])+") OR (FIX_ID_sp1='"+str(fix_g[1])+"' AND FIX_ID_sp2="+str(fix_g[0])+")") #Change DB to genome pairwise distance (bidirectional-blast best)
							result_fix_g= c.fetchall()						
							list_Gdistance.write(str(fix_g[0])+'\t'+str(fix_g[1])+'\t'+str(result_fix_g[0][0])+'\n') #write to be input of building phylogenetic tree 
						list_Gdistance.close()
						build_tree = build_phylo_tree(file_name_pair_dist) #Total_genomes, Total_pair_distance, check1, total_branch_length
						total_branch_len = build_tree[3] #Conservation score - total branch length


						check1=build_tree[2] #check the number of pair distance found in build_tree.py
						if len(get_fixID) == build_tree[0] :
							check2=True #check the number of genome found between pairwise_distance_i.temp and neighbor_id_i.temp2

						os.system('rm '+file_name_pair_dist)

					else:
						total_branch_len =0.0 #Neighbors are conserved in only 1 genome , so No Distance score 


					#query_cluster_id neighbor_cluster_id neighbor_gene_found qeury_gene_found number_of_genome(fix_id)_found	total_branch_lenngth(Phylogeny)				
					out_conservation_file=open('Conservation_score_Qid_'+str(i)+'.out','a')
					out_conservation_file.write(str(i)+'\t'+str(freq[1])+'\t'+str(freq[0])+'\t'+str(Total_query_gene)+'\t'+str(len(get_fixID))+'\t'+str(check1)+'\t'+str(check2)+'\t'+str(total_branch_len)+'\n' )
					out_conservation_file.close()
				
				os.system('rm cluster_'+str(i)+'.nout')
				os.system('rm neighbor_id_'+str(i)+'.temp')
				os.system('rm neighbor_id_'+str(i)+'.temp2')
				os.system('rm neighbor_id_'+str(i)+'.temp3')
	conn.close()	

def blast_query_to_MyDB(list_query_fasta,DB_name,evalue_cutoff):#DB_name our photosynthetic genes in blast DB format, max_seq =max_target_seqs

	#blastp using query gene file
	blastp="blastp -db "+DB_name+" -query "+list_query_fasta+" -out matched2DB.query -max_hsps 1 -evalue "+str(evalue_cutoff)+" -outfmt '6 qseqid sseqid evalue pident qstart qend sstart send' -num_threads 4"
	print "#blastp using query gene file"
	print blastp
	os.system(blastp)
	print "blastp finished!!!"

def get_cluster_ID(seq_name_file,cds_MCL_cluster): #edited

	import commands
	import sqlite3
	conn = sqlite3.connect("Photosyn.db")
	c = conn.cursor()

	w=open('query2myClusterID.out','w')
	for prot in open(seq_name_file):
		prot=prot.replace('\n','')
		prot=prot.split()

		c.execute("SELECT cluster_id FROM "+cds_MCL_cluster+" WHERE seq_id = '"+prot[1]+"'")
		result = c.fetchall()
		for i in result:
			w.write(str(prot[0])+'\t'+str(prot[1])+'\t'+str(i[0])+'\n')
	w.close()

	conn.close()

def copyneighbor_from_MyDB(query2myClusterID,path_to_all_neighbor,conserved_cluster): #to copy only gene neighbors for single query

	for i in open('query2myClusterID.out'):
		i=i.replace('\n','')
		i=i.split()

		if int(i[2]) <= int(conserved_cluster) : #24240
			
			nb_list=commands.getoutput('grep "^'+i[1]+'" '+path_to_all_neighbor+'query_prot_id_'+str(i[2])+'.ext.neighbor')
			if nb_list =='':
				print 'No neighbor FOUND for ', i[0], 'using ', i[1] , 'as virtual query gene.'
			else:
				w=open('query_prot_id_'+str(i[2])+'.ext.neighbor','w') ######
				w.write(nb_list)
				w.close()

def copyneighborcluster_from_MyDB(query2myClusterID,path_to_all_neighbor,conserved_cluster): #to copy the entire cluster neighbor

	for i in open('query2myClusterID.out'):
		i=i.replace('\n','')
		i=i.split()

		if int(i[2]) <= int(conserved_cluster) : #24240
			if os.path.exists(path_to_all_neighbor+'query_prot_id_'+str(i[2])+'.ext.neighbor'):
				commands.getoutput('cp '+path_to_all_neighbor+'query_prot_id_'+str(i[2])+'.ext.neighbor ./')
				commands.getoutput('cp '+path_to_all_neighbor+'query_prot_id_'+str(i[2])+'.ext ./')
			else:
				print '++++++++++++ No neighbor for cluster ',i[2]



cluster_file='cluster_export'
if int(commands.getoutput("grep '' -c unmatched_prot_file.query")) > 0:

	get_cluster_ID("unmatched_prot_file.query",DB_mcl_tablename)
	#commands.getoutput("rm matched2DB.query")

	#case1: call neighbor gene from already called genes
	evalue_foldernum=sys.argv[2].replace('cds_MCL_cluster_evalue','') #get all neighboring folder number
	path_to_all_neighbor='./call_neighbor_3evalue/e'+str(evalue_foldernum)+'/'
	copyneighborcluster_from_MyDB('query2myClusterID.out',path_to_all_neighbor,sys.argv[1])
	#####################################################################################################

	#case2: call new neighbor genes (remove only ### to run )
	#input file
	#os.system('cut -f1 '+list_prodigal_photo+' >get_photo_gene.out2')
	os.system('cut -f2 query2myClusterID.out >get_photo_gene.out2')
	#1 make input file to neghbor calling, 2 column file, 1)cluster_id 2)gene_query , namely (cluster_export_Ref5.txt)
	multiple_call_cluster_id('get_photo_gene.out2',cluster_file+'.txt') 
	#2 use n.py to call neighbor genes using cluster_id as input
	###os.system('python n.py "'+cluster_file+'.txt" 250 '+DB_mcl_tablename+' ') #specify input file and gene neighbor distance(bp) and database name of Cluster ID

	###os.system('mv *.cyto *.ext *.ext.neighbor ./all_neighboring')
	#####################################################################################################

	os.system('mkdir ./all_neighboring')
	os.system('mv *.ext *.neighbor ./all_neighboring')

	#3 calling cluster_id of neighboring genes, 2 column file, 1)cluster_id 2)gene_neighboring_query ()
	os.system('cut -f1 '+cluster_file+'.txt|sort -u >'+cluster_file+'_uniq.txt')
	call_neighboring_gene(cluster_file+'_uniq.txt') 

	#os.system('python get_photo_uniprot_id.py >cds_uniprot.txt')#
	#os.system('python merge2file.py '+cluster_file+'.txt cds_uniprot.txt >combined_cluster_id_uniprot.txt')#

	os.system('tail -n +2 neighbor_enrichment.out > neighbor_enrichment_no_header.out')
	os.system('python make_temp_database.py')

	call_neighboring_gene_conservation_phylo(cluster_file+'_uniq.txt')
	os.system('mkdir conservation_score')
	os.system('mv Conservation_score_Qid_*.out ./conservation_score')
	os.system('cat ./conservation_score/Conservation_score_Qid_*.out > ./Conservation_score_Qid_all.out')
	insert_to_DB('Conservation_score_Qid_all.out',3)
	
	
else:
	print 'No blast matched'


