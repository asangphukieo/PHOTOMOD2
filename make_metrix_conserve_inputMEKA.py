import commands
import sqlite3
import sys
#edited to ouput for MEKA
def range_quartile(Q1,Q3,input):
	if input < Q1:
		return 0
	elif input >= Q1 and input < Q3:
		return 1
	else :
		return 2

def make_metrix_conservation(cluster_export,DB_conserve_score,cutoff_score,model_nb_idfile,label):
	
	conn = sqlite3.connect("GO_neighbor_enrichment.db")
	c = conn.cursor()
	header_nb_id=''
	model_nb=[]
	for k in open(model_nb_idfile):
		k=k.replace('\n','')
		model_nb.append(k)
		header_nb_id+=str(label+str(k))+', '
	
	header_nb_id=header_nb_id[:-2]
	w2file=open('metrix_cv_nb_predicted.csv','w')

	w2file.write('#Query, '+header_nb_id+'\n')
	w2file.close()
	for query in open(cluster_export):	#read prodigal prot and its cluster id in the file 'cluster_export.txt'
		query=query.replace('\n','')
		query=query.split() #query[0]=>cluster ID, query[1]=>prot name
		#print query[1]

		list_nb_row = [0]*len(model_nb) #make list of zero for replacing by conservation score

		for model_nb_id in model_nb:
			
			c.execute("SELECT neighbor_cluster_id, conservation_score FROM neighbor_conservation_phylo WHERE neighbor_cluster_id = "+str(model_nb_id)+" AND Query_cluster_id = "+query[0]+" AND conservation_score>"+str(cutoff_score)+" ")
			result = c.fetchall()
			
			if len(result) >= 1:
				for i in result:
					try:
						# i[0]=> nb cluster ID , i[1] => conservation_score
						if sys.argv[2] == 'e10.C.':
							q1=0.267
							q3=2.306
						elif sys.argv[2] == 'e50.C.':
							q1=0.352
							q3=2.252
						elif sys.argv[2] == 'e100.C.':
							q1=0.252
							q3=1.642
						else:
							print 'Please check feature label'

						level_cv=range_quartile(q1,q3,float(i[1]))

						index_match=model_nb.index(str(i[0]))
						list_nb_row[index_match] = level_cv
					except ValueError:
						print "Not include cluster ID : "+str(i[0])
			
							
		list_nb_row=str(list_nb_row)
		list_nb_row=list_nb_row.replace('[','')
		list_nb_row=list_nb_row.replace(']','')
		w2file=open('metrix_cv_nb_predicted.csv','a')
		w2file.write(query[1]+', '+list_nb_row+'\n')
		w2file.close()


	conn.close()

make_metrix_conservation('cluster_export.txt','GO_neighbor_enrichment.db',sys.argv[3],sys.argv[1],sys.argv[2]) #'e10.C.uniq_nb_id.id'
