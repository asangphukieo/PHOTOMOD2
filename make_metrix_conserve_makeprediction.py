import commands
import sqlite3
import numpy,sys

def range_quartile(Q1,Q3,input):
	if input < Q1:
		return 1
	elif input >= Q1 and input < Q3:
		return 2
	elif input >= Q3 :
		return 3

def make_metrix_conservation(cluster_export,evalue,model_nb_idfile,add_label,Q1,Q3):
	
	Q1=float(Q1)
	Q3=float(Q3)
	label='e'+str(evalue)+'.C.'
	conn = sqlite3.connect("GO_neighbor_enrichment.db")
	c = conn.cursor()

	header_nb_id=''
	list_uniq_nb=[]
	for k in open(model_nb_idfile):
		k=k.replace('\n','')
		list_uniq_nb.append(k)
		header_nb_id+=str(label+str(k))+', '
	
	header_nb_id=header_nb_id[:-2]
	w2file=open('metrix_cv_nb_predicted_e'+str(evalue)+'.csv','w')
	
	if add_label == 'no_label':
		w2file.write('#Query, '+header_nb_id+'\n')
		w2file.close()
	else:
		w2file.write('#Query, '+header_nb_id+', label'+'\n')
		w2file.close()

	#c.execute("SELECT conservation_score FROM neighbor_conservation_phylo WHERE conservation_score>0 ")
	#cutoff_Q = c.fetchall()
	#Q1=numpy.percentile(cutoff_Q, 25) #Quantile 1
	#Q3=numpy.percentile(cutoff_Q, 75) #Quantile 3

	print 'Quantile 1 of conservation score = ',Q1,', Quantile 3 of conservation score = ',Q3

	c.execute("SELECT neighbor_cluster_id FROM neighbor_conservation_phylo WHERE conservation_score>0 ")
	list_neighbor_id = c.fetchall()

	for query in open(cluster_export):	#read prodigal prot and its cluster id in the file 'cluster_export.txt'
		query=query.replace('\n','')
		query=query.split() #query[0]=>cluster ID, query[1]=>prot name
		#print query[1]
		c.execute("SELECT neighbor_cluster_id, conservation_score FROM neighbor_conservation_phylo WHERE Query_cluster_id = "+query[0]+" AND conservation_score>0 ")
		result = c.fetchall()
		
		list_nb_row = [0]*len(list_uniq_nb) #make list of zero for replacing by conservation score
		#print result
		for i in result:
			try:
				# i[0]=> nb cluster ID , i[1] => conservation_score
				index_match=list_uniq_nb.index(str(i[0]))
				level_cv=range_quartile(Q1,Q3,float(i[1])) ############ set quartile range 1 = low, 2= med , 3 =high conservation
				list_nb_row[index_match] = level_cv
			except ValueError:
				print "Not include cluster ID : "+str(i[0])

		list_nb_row=str(list_nb_row)
		list_nb_row=list_nb_row.replace('[','')
		list_nb_row=list_nb_row.replace(']','')
		w2file=open('metrix_cv_nb_predicted_e'+str(evalue)+'.csv','a')

		if add_label == 'no_label':
			w2file.write(query[1]+', '+list_nb_row+'\n')
			w2file.close()
		else:
			w2file.write(query[1]+', '+list_nb_row+', '+str(add_label)+'\n')
			w2file.close()
	conn.close()

make_metrix_conservation('cluster_export.txt',sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]) ###### 1) evalue, 2) neighbor id , 3) label no_label = no add label, 4) quartile 1, 5) quaretile 3
