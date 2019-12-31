#to select the output gene above the criteria
#v2 add scov criteria cutoff

import os
import re

#Tempolary file 
#out.blast  >> raw output from blast with option -outfmt "7 qseqid sseqid evalue pident qstart qend qlen slen length" -evalue 0.001 
#out1.blast >> select only gene_name and its score
#out2.blast >> select only gene that have evalue lower than the criteria
#out3.blast >> select only gene that have identity more than the criteria
#out4.blast >> select only gene that have coverage more than the criteria

def select_gene(file_name,qiden, qcov, scov, evalue):
	os.system("awk '! /#/' "+file_name+" > out1.blast") #Set input file name
	outfile="blast_photo_cutoff.out" #Set output file name
	
	if os.stat("out1.blast").st_size==0:
		print "No gene pass the E-value criteria"

		
	else :
		
		os.system("awk '$3<="+str(evalue)+"' out1.blast > out2.blast")
		os.system("awk '$4>="+str(qiden)+"' out2.blast > out3.blast")
		if os.stat("out3.blast").st_size==0:
			print "No gene pass the identity criteria"
		else:
			
			output_file=open("out3.blast",'rt')
			output_file_open=output_file.read()
			output_file.close()
			output_file_read = output_file_open.split('\n')

			del output_file_read[len(output_file_read)-1]
		
			output_1=''
			output_2=''
			write_to_file=open(outfile,'w') 
			write_to_file.close()
			
			for line in open("out3.blast"): #read file as object for very large file avoiding memory error
				fields = re.split(r'\t+', line.strip())
				q_start, q_end, q_len = map(float, (fields[4], fields[5], fields[6]))
				q_cov = 100 * ((q_end - q_start)+ 1) / q_len

				s_start, s_end, s_len = map(float, (fields[9], fields[10], fields[7]))
				s_cov = 100 * ((s_end - s_start) + 1) / s_len

				if q_cov >= qcov and s_cov >= scov:
					#number_genome=genome.replace('prodigal.orf.fsa.blastout','')
					#Sseq = fields[12].replace('-','')
					#Qseq = fields[11].replace('-','') #for report sequence

					
					#output_1= output_1+'>FX'+number_genome+fields[1]+'\n'+Sseq+'\n'

					#report cytoscape file relation
					# Fields: 0 query id, 1 subject id, 2 evalue, 3 % identity, 4 q. start, 5 q. end, 6 query length, 7 subject length, 8 alignment length, 9 s. start, 10 s. end, 11 query seq, 12 subject seq
					output_2= fields[0]+'\t'+fields[1]+'\t'+fields[2]+'\t'+fields[3]+'\t'+str(q_cov)+'\t'+str(s_cov)+'\t'+fields[6]+'\t'+fields[7]+'\n'
					write_to_file=open(outfile,'a') 
					write_to_file.write(output_2)
					write_to_file.close()

			#write_to_file=open("blastall.fasta",'w') 
			#write_to_file.write(output_1)
			#write_to_file.close()
			
			
			
			#if os.stat("blastall.fasta").st_size==0:
			#	print "No gene pass the coverage criteria"

	#for remove intermediate files	
	#os.system("rm out.blast")
	#os.system("rm out2.blast")
	#os.system("rm out3.blast")
	
#select_gene(file_name,qiden, qcov, scov, evalue)

#if user do not want to set any value, just set to the lowest criteri
#such as in the case of; user want to display only sequences that have identity more than 50%
#just set => select_gene(50.0, 0.0, 1)


#select_gene('genome_5.blastout',80.0, 80.0, 80.0, 1e-10)
	
