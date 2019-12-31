from Bio.Phylo.TreeConstruction import _Matrix,DistanceTreeConstructor,_DistanceMatrix
from Bio import Phylo
import commands,sys
#Ref:http://biopython.org/wiki/Phylo
#Require Biopython package

def build_phylo_tree(file_name_pair_dist): #to return total branch lenght
	file_pair_dist=file_name_pair_dist	#all.reciprocal
	pair_dist = commands.getoutput("cut -f1,2,3 "+file_pair_dist+" ") #need to delete header of the table
	pair_dist = pair_dist.split('\n')
	list_genome = commands.getoutput("awk '{print $1}{print $2}' "+file_pair_dist+" |sort -g|uniq")
	list_genome = list_genome.split('\n')

	#print 'Total genomes >> ',len(list_genome)
	Total_genomes= len(list_genome)
	#print 'Total pair distance >> ',(len(list_genome)*(len(list_genome)-1))/2
	Total_pair_distance=(len(list_genome)*(len(list_genome)-1))/2
	check1=False
	if len(pair_dist) == (len(list_genome)*(len(list_genome)-1))/2 : #To check the number of pair distance
		#print 'Pass check 1 : total pair distances are correctly found '
		check1=True

	matrix =[]
	for n in range(1,len(list_genome)+1 ):
		matrix.append([0]*n)
	#print matrix
	Total_metrix=len(matrix)
	#print 'Total matrix >>',len(matrix)

	Ds=_DistanceMatrix(list_genome,matrix)

	for pair in pair_dist:
		i=pair.split()
		#print i
		Ds[i[0],i[1]] = float(i[2])
	#print Ds
	constructor = DistanceTreeConstructor()
	tree = constructor.nj(Ds)
	#print (tree) #for visualization
	#Phylo.draw(tree, branch_labels=lambda c: c.branch_length)

	print 'total_branch_length >> ', tree.total_branch_length()
	return Total_genomes, Total_pair_distance, check1, tree.total_branch_length()

#build_phylo_tree('scaled_16sdistance_154_forDB.csv')
