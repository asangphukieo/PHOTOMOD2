import os
import re
import commands
import sys

print 'input file >>',sys.argv[1]
print 'Intergenic distance condition(bp) >>',sys.argv[2]
print 'Database table name >>',sys.argv[3]
#update : use for only select neighborhood (skip blast process)

##############################################################################################################	
##############################################################################################################	
#Version7 fix bug, loss some diverging neighobor gene (replace < with <= and remove method add method split_gene_position specific for divergent transciption)
#to select the neighboring gene around the query gene
#add input from user
#Version6 connect to photosyn db 
#Version5 modify output report , remove > sign from query name and neighbor, add Diverging_neighbor type
#Version4 Include neighboring gene from diverging transcription
#Version3 complete neighboring genes selection in condition of in_range X bp, overlaped frame, circular genome, read both strands 

def create_arrange_genome(genome_id):
	#original file is arrange_orf.py
	#To arrange orf file, separate 2 different strand
	#strand_1 for 5' to 3'
	#strand_2 for 3' to 5'
	file_name="./OnlyGenome/"+genome_id+"_prodigal.orf.fsa" # need to edit !!!
	os.system("awk '/>/' "+file_name+" > out.orf_name")
	strand_1=''
	strand_2=''
	for i in open('out.orf_name'):
		i=i.replace('\n','')
		col=i.split('_')
		cds=col[len(col)-1]
		cds_strand=cds.split('-')
		if int(cds_strand[1])-int(cds_strand[0]) >0:
			strand_1=strand_1+i+'\n'			
		else:
			strand_2=strand_2+i+'\n'

	write_to_file=open("strand_1.orf_name",'w') 
	write_to_file.write(strand_1)
	write_to_file.close()
	write_to_file=open("strand_2.orf_name",'w') 
	write_to_file.write(strand_2)
	write_to_file.close()

	#for remove intermediate files
	os.system("rm out.orf_name")	


def select_neighborhood(query_first_input,input_name,chr_id,dist,circle, chr_len, out_cyto, min_operonic, operonic_dist):
	Cytoscape=""
	if os.path.exists(query_first_input+'.neighbor') == False:
		out_neighbor_list=open(query_first_input+'.neighbor','w')
		out_neighbor_list.close()
	out_neighbor_list=open(query_first_input+'.neighbor','a')	
	
	query=input_name
	genome_id=query.split('_')
	
	#To arrange orf file
	create_arrange_genome(chr_id)	
	#query is AE006470_AE01278_CDS_203664-202288

	num_gene=query
	tab=query.split('_')
	gene_position_str_end=tab[len(tab)-1].split('-')
	query_str = gene_position_str_end[0]
	query_end = gene_position_str_end[1]

	out_neighbor_list.write("#Find gene neighborhood for query "+num_gene+'\n')
	#To check which strand
	if int(query_end)-int(query_str) > 0:	

		out_neighbor_list.write("#read orf file strand 1\n")
		
		neighboring_gene=neighboring(dist,circle, chr_len,query, "plus", "strand_1.orf_name", query_str, query_end, num_gene,'neighbor',query, '000')
		out_neighbor_list.write(neighboring_gene[0])
		Cytoscape = Cytoscape+neighboring_gene[1]

		if len(neighboring_gene[0]) ==0:
			out_neighbor_list.write("#No neighboring gene found\n")

			out_neighbor_list.write("#==============================================================================\n")

		else:	
			first_cds = split_position(neighboring_gene[4], 'plus')
			first_cds_str = first_cds[0]
			first_cds_end = first_cds[1]
		
			#search_gene_in_range(str_position, operonic_distance, oposit_strand, first_cds_str)
			divergent=search_gene_in_range(first_cds_str-min_operonic, operonic_dist-min_operonic, "minus", circle, chr_len, first_cds_str)
			if divergent[0] == None:
				out_neighbor_list.write("#No diverging transcription found\n")
				out_neighbor_list.write("#==============================================================================\n") 
			else:
				divergent_position = split_gene_position(divergent[0]) #method (neighboring) will invert str and end position again
				diverging_str = divergent_position[0]
				diverging_end = divergent_position[1]

				neighboring_gene_diverging = neighboring(dist,circle, chr_len, divergent[0], "minus", "strand_2.orf_name", diverging_str, diverging_end, num_gene+' Diverging neighborhood','Diverging_neighbor',query,divergent[1])

				Cytoscape = Cytoscape+neighboring_gene_diverging[1]
				out_neighbor_list.write('#Found Diverging trascription near '+neighboring_gene[4]+'\n')
				out_neighbor_list.write('#Diverging trascription\n'+neighboring_gene_diverging[0])
				out_neighbor_list.write("#==============================================================================\n")

	elif int(query_end)-int(query_str) < 0:
		out_neighbor_list.write("#read orf file strand 2\n")
		neighboring_gene=neighboring(dist,circle, chr_len,query, 'minus', "strand_2.orf_name", query_str, query_end, num_gene,'neighbor',query,'000')
		out_neighbor_list.write(neighboring_gene[0])
		Cytoscape=Cytoscape+neighboring_gene[1]
	

		if len(neighboring_gene[0]) ==0:
			out_neighbor_list.write("#No neighboring gene found\n")
			out_neighbor_list.write( "#==============================================================================\n")
		else:	
			first_cds = split_position(neighboring_gene[4], 'minus')
			first_cds_str = first_cds[0]
			first_cds_end = first_cds[1]
		
			#search_gene_in_range(str_position, operonic_distance, oposit_strand, first_cds_str)
			divergent=search_gene_in_range(first_cds_str+min_operonic, operonic_dist-min_operonic, "plus", circle, chr_len, first_cds_str)#if min_operonic = 200, operonic_dist should be operonic_dist-min_operonic or if operonic dist = 1000 then 1000-200=800
		
			if divergent[0] == None:
				out_neighbor_list.write("#No diverging transcription found\n")
				out_neighbor_list.write( "#==============================================================================\n")
			else:
				divergent_position = split_gene_position(divergent[0])  #method (neighboring) will invert str and end position again
				diverging_str = divergent_position[0]
				diverging_end = divergent_position[1]
				neighboring_gene_diverging = neighboring(dist,circle, chr_len, divergent[0], "plus", "strand_1.orf_name", diverging_str, diverging_end, num_gene+' Diverging neighborhood','Diverging_neighbor',query,divergent[1])

				Cytoscape =  Cytoscape+neighboring_gene_diverging[1]
				out_neighbor_list.write('#Found Diverging trascription near '+neighboring_gene[4]+'\n')

				out_neighbor_list.write('#Diverging trascription\n'+neighboring_gene_diverging[0])
				out_neighbor_list.write("#==============================================================================\n")
	
	else:
		print "Query input is incorrect"

	if out_cyto == True:
		write_to_file=open(query_first_input+'.cyto','w')
		write_to_file.write(Cytoscape)		
		write_to_file.close()	
		#print Cytoscape
	out_neighbor_list.close()
	
############################################
#neighboring(dist,circle, chr_len,col (query full name), "minus", "strand_1.orf_name", query_str, query_end, num_gene_of_query, original_query_gene, operonic_distance)
def neighboring(dist,circle, chr_len,col, strand, file_input, query_str, query_end, num_gene, type, original_query,operonic_dist):
	col=col.replace('>','')	
	neighbor=""
	Cytoscape=""
	cds_str=[]
	cds_end=[]
	first_cds=""
	query_str_in_def = int(query_str)
	query_end_in_def = int(query_end)
	strand_sign='+'
	if strand == 'minus':
		strand_sign='-'
		query_str_in_def = int(query_end)
		query_end_in_def = int(query_str)

	file_orf1=open(file_input,'r') 
	file_orf1_out=file_orf1.read()
	orf1_out=file_orf1_out.split('\n')
	file_orf1.close()
	del orf1_out[len(orf1_out)-1]
	
	#for report start gene of diverging transcription
	if type == 'Diverging_neighbor':
		original_query=original_query.replace('>','')
		neighbor=neighbor+original_query+'\t'+col+'\t'+strand_sign+'\t'+"Separated_"+str(operonic_dist)+'_bp\t'+'*'+'\t'+type+'\t'+'\n'	
		col=original_query	###
	
	k=0	# k is position of pivot 
	while k <= len(orf1_out)-1:
		
		db_match = split_position(orf1_out[k],strand)
		db_str = db_match[0]
		db_end = db_match[1]

		#To match between query and orf file
		#print query_str_in_def, db_str, query_end_in_def , db_end
		if query_str_in_def == db_str and query_end_in_def == db_end:

			first_cds= orf1_out[k]
			print "#Query gene "+str(num_gene)+" Found"
									
			cds_str.append(int(db_str))
			cds_end.append(int(db_end))
			
			#Right direction >> 
			l=1	# l is position of neighborhood according to the arranged genes in file strand.orf_name  
			#k postion of matched gene
			if k == len(orf1_out)-1:
				l=0
			while k+l <= len(orf1_out)-1:
				#k+l-1 first next position
				db = split_position(orf1_out[k+l-1],strand)
				db_str = db[0]
				db_end = db[1]

				#k+l second next position
				db2 = split_position(orf1_out[k+l],strand)
				db2_str = db2[0]
				db2_end = db2[1]
				
				db3 = split_position(orf1_out[0],strand)
				db3_str = db3[0]
				db3_end = db3[1]

				if int(db_str)<= int(db2_str) <= (int(db_end)+int(dist)) :	####fix	
					
					intergenic_dist=int(db2_str) - int(db_end) #+1 for first bp counting ####fix	
					
					#Checking neighborhood conditions
					overlap='Separated_'+str(intergenic_dist)+'_bp'
					
					#to check wether overlaped gene [3]
					if int(db_str)<= int(db2_str) <= int(db_end) :
						overlapped_dist= int(db_end) - int(db2_str) +1
						overlap='Overlapped_'+str(overlapped_dist)+'_bp'
						
					########### OUTPUT UPDATE ###########	
					neighbor_name2=orf1_out[k+l].replace('>','')
					neighbor=neighbor+col+'\t'+neighbor_name2+'\t'+strand_sign+'\t'+overlap+'\t'+'>'+'\t'+type+'\n'
					Cytoscape=Cytoscape+orf1_out[k+l-1]+'\t'+orf1_out[k+l]+'\n'
					if strand == 'minus': 	#To save the first cds of operon
						first_cds= orf1_out[k+l]
					db2=split_position(orf1_out[k+l],strand)
					cds_str.append(int(db2[0]))
					cds_end.append(int(db2[1]))	
					#Checking circular genome neighborhood
					if circle == True and k+l == len(orf1_out)-1 and int(db3_str) <= (int(db2_end) + dist ) - chr_len :			
						#need to separate condition because k+l == len(orf1_out)-1 is not in first loop	
						cir_intergenic_dist = ( chr_len - int(db2_end) ) + int(db3_str)

						########### OUTPUT UPDATE ###########
						neighbor_name=orf1_out[0].replace('>','')
						neighbor=neighbor+col+'\t'+neighbor_name+'\t'+strand_sign+'\t'+"Separated_"+str(cir_intergenic_dist)+'_bp'+'\t'+'>'+'\t'+type+'\t'+"Circular_neighbor"+'\n'
						Cytoscape=Cytoscape+orf1_out[-1]+'\t'+orf1_out[0]+'\n'
						cds_str.append(int(db3_str))
						cds_end.append(int(db3_end))
						
						if strand == 'minus':
							first_cds= orf1_out[0]

						l= -k+1	# to reset l for covering circular genome
						
					l=l+1
							
				else :
					break
			
					

			#Left direction <<
			l= 0	# l is position of neighborhood according to the arranged genes in file strand.orf_name  
			while k+l <= k :#
				db = split_position(orf1_out[k+l],strand)
				db_str = db[0]
				db_end = db[1]

				db2 = split_position(orf1_out[k+l-1],strand)
				db2_str = db2[0]
				db2_end = db2[1]	
				
				ex_len= int(db_str) - dist #excessive chromosome
				ex_len_chr= chr_len + ex_len
				
				if int(dist) < int(db_str) and (k+l-1 != 0) and (int(db_str)-int(dist)) <= int(db2_end) <= int(db_end): ####fix #Checking whether first position and circular chlomosome  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

					intergenic_dist=  int(db_str) - int(db2_end) ####fix
					overlap='Separated_'+str(intergenic_dist)+'_bp'
					
					#to check wether overlapped gene
					if int(db_str)<= int(db2_end) <= int(db_end):
						overlapped_dist = int(db2_end) - int(db_str) +1
						overlap='overlapped_'+str(overlapped_dist)+'_bp'
						
					########### OUTPUT UPDATE ###########	
					neighbor_name3=orf1_out[k+l-1].replace('>','')
					neighbor=neighbor+col+'\t'+neighbor_name3+'\t'+strand_sign+'\t'+overlap+'\t'+'<'+'\t'+type+'\n'
					Cytoscape=Cytoscape+orf1_out[k+l]+'\t'+orf1_out[k+l-1]+'\n' 
					if strand == 'plus': 	#To save the first cds of operon
						first_cds= orf1_out[k+l-1]	
				
					cds_str.append(int(db2_str))
					cds_end.append(int(db2_end))
					
					l=l-1	
				
				elif circle == True and int(dist) >= int(db_str)  and k+l == 0 and int(ex_len_chr) <= int(db2_end) <= int(chr_len): #Checking whether first position and circular chlomosome
					#because array position can be a minus sign e.g. array[-1] >> the last array , no need to reset 'l' like above
					cir_intergenic_dist = int(db_str) + ( chr_len - int(db2_end) )
					
					########### OUTPUT UPDATE ###########
					neighbor_name4=orf1_out[k+l-1].replace('>','')
					neighbor=neighbor+col+'\t'+neighbor_name4+'\t'+strand_sign+'\t'+"Separated_"+str(cir_intergenic_dist )+'_bp'+'\t'+'<'+'\t'+type+'\t'+"Circular_neighbor"+'\n'
					Cytoscape=Cytoscape+orf1_out[k+l]+'\t'+orf1_out[k+l-1]+'\n'

					if strand == 'plus': 	#To save the first cds of operon
						first_cds= orf1_out[k+l-1]

					cds_str.append(int(db2_str))
					cds_end.append(int(db2_end))
					
					l=l-1	
				else:
					break

		#else:
			#print "No gene match"
		k=k+1
	#neighbor; list of neighboring gene and details, Cytoscape; list of neighboring gene in cytoscape relationship format, cds_position; position of genes
	return neighbor, Cytoscape, cds_str, cds_end, first_cds

#to convert name to only location
def split_position(name,strand):
	tab2=name.split('_')
	db2_position_str_end=tab2[len(tab2)-1].split('-')
	db2_str = int(db2_position_str_end[0])
	db2_end = int(db2_position_str_end[1])
	if strand == 'minus':
		db2_str = int(db2_position_str_end[1])
		db2_end = int(db2_position_str_end[0])
	return db2_str, db2_end

def split_gene_position(name):
	tab2=name.split('_')
	db2_position_str_end=tab2[len(tab2)-1].split('-')
	db2_str = int(db2_position_str_end[0])
	db2_end = int(db2_position_str_end[1])
	return db2_str, db2_end


def search_gene_in_range(db_position, distance, strand, circular, chr_len,first_cds_str): #first_cds_str= position of query gene for searching
	accual_operonic_dist=None
	list_diverging=None
	filename="strand_1.orf_name"
	if strand == "minus":
		filename="strand_2.orf_name"
	open_arrange_file=open(filename,'r')
	arrange_file=open_arrange_file.read()
	arrange_file=arrange_file.split('\n')
	open_arrange_file.close()
	del arrange_file[len(arrange_file)-1]
	#db_position is strat position of query which run on different strand

	index_count = 0
	while index_count <= len(arrange_file)-1 and index_count >= -(len(arrange_file)-1):
		position=split_position(arrange_file[index_count],strand)
				
		if strand == 'minus':
			#run on minus strand
			if circular == True and (db_position-distance) < 0: #for checking circular chr.
				if int(position[1]) in range(chr_len+(db_position-distance), chr_len): #circular genome
					list_diverging = arrange_file[index_count] #return only nearest gene of the query
					accual_operonic_dist =  chr_len - int(position[1]) + int(first_cds_str) +1
					break
			else:
				if int(position[1]) in range(db_position-distance, db_position):
					list_diverging = arrange_file[index_count] #return only nearest gene of the query
					accual_operonic_dist = int(first_cds_str) -  int(position[1]) +1
					break
			index_count = index_count-1
		else:
			#run on plus strand
			if circular == True and (db_position+distance) > chr_len: #distance for plus strand is positive value
				if int(position[0]) in range(0, (db_position+distance)-chr_len): #circular genome
					list_diverging = arrange_file[index_count] #return only nearest gene of the query
					accual_operonic_dist =  chr_len - int(first_cds_str) + int(position[0]) +1
					break
			else:
				if int(position[0]) in range(db_position, db_position+distance):
					list_diverging = arrange_file[index_count] #return only nearest gene of the query
					accual_operonic_dist = int(position[0]) - int(first_cds_str) +1
					break
			index_count = index_count+1
	return list_diverging, accual_operonic_dist
##############################################################################################################	
##############################################################################################################	

##########################
#write from list to file
def write2file(list_1,out_1):
	file_out=open(out_1,'w')
	for i in list_1:
		line=''
		for j in i:		
			line=line+str(j)+'\t'
		file_out.write(line+'\n') 
	file_out.close()
#########################

import sqlite3
conn = sqlite3.connect("Photosyn.db")
c = conn.cursor()
DB_mcl_tablename=sys.argv[3]
#blast_query_to_genome("list_genomes.txt", "list_query.txt")
#report_file(0,0,1e-10) #set criteria for qiden, qcov, evalue
#select_neighborhood(250, True , 4986812, True, 200, 1000)  #select_neighborhood(query_list_input,query_name,FIX_id,intergenic_dist,circle_chromosome,chromosome_length, Cytoscape_output, interoperonic_min_distance, interoperonic_max_distance)


#Make quey protein id file =>' default is query_prot_id '
#os.system('ls ./query_cluster/ > list_query_prot_id')
file_query_cluster=sys.argv[1]

os.system('cut -f1 '+file_query_cluster+' |sort -u > query_cluster_id.ext')

for cluster_id in open('query_cluster_id.ext'):
	cluster_id=cluster_id.replace('\n','')
	
	#os.system('cp ./query_cluster/'+cluster_id+' ./query_prot_id')#optional
	#in this case , I extract cluster id from database to be queries
	print 'Extract protein sequences from cluster ID '+str(cluster_id)
	c.execute("SELECT Genome_accession.FIX_ID,cds.genome_acc,Chr_lenght,Chr_circulation,"+DB_mcl_tablename+".seq_id  FROM "+DB_mcl_tablename+",cds,Genome_accession WHERE cluster_id="+str(cluster_id)+" AND cds.cds_ID="+DB_mcl_tablename+".seq_id AND Genome_accession.genome_acc=cds.genome_acc")
	result = c.fetchall()
	
	query_first_input='query_prot_id_'+str(cluster_id)+'.ext'
	write2file(result,query_first_input)

	for genomes in open(query_first_input):		
		genomes=genomes.replace('\n','')		
		col=genomes.split()

		if col[3] == 'circular':
			cir=True
		else:
			cir=False
		
		gen=int(sys.argv[2])
		op_min=200
		op_max=1000
		
		print 'Read genome	: '+col[1]+' (FIX_ID '+str(col[0])+')'
		print 'Chromosome length	: '+str(col[2])
		print 'Circular chromosome	: '+str(cir)
		print 'Intergenic distance cutoff	: '+str(gen)
		print 'Interoperonic minimum distance	: '+str(op_min)
		print 'Interoperonic maximum distance	: '+str(op_max)
		
		select_neighborhood(query_first_input,col[4],str(col[0]),gen, cir , int(col[2]), True, op_min, op_max)
		#select_neighborhood(query_list_input,query_name,FIX_id,intergenic_dist,circle_chromosome,chromosome_length, Cytoscape_output, interoperonic_min_distance, interoperonic_max_distance)
		

conn.close()

