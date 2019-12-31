import commands
import sys
#fix label != 'no_label'


def concat_feature(file1,file2,label): #concat the row of each file that have the same label in the first collumn,1)first fie 2) second file ,3 additional label (default = no_label)
	#print file1,file2
	headfile2_cid=commands.getoutput('head -n1 '+file2+'| cut -d"," -f2-')
	headfile2=headfile2_cid.replace('\n','')
	headfile2=headfile2.split(',')
	num_feature2=len(headfile2)
	#print num_feature2
	list_zero=[0]*num_feature2
	list_zero=str(list_zero)
	list_zero=list_zero.replace('[','')
	list_zero=list_zero.replace(']','')
	if label != 'no_label':
		list_zero=list_zero[:-1]
		list_zero=list_zero+label

	for i in open(file1):
		i=i.replace('\n','')
		row=i.split(',')
		if '#' in row[0]:
			print i+','+headfile2_cid

		else:
			matchfile2=commands.getoutput('grep -m 1 -w "^'+row[0]+'" '+file2+'|cut -d"," -f2-')
			matchfile2=matchfile2.replace('\n','')

			if matchfile2 != '':
				print i+','+matchfile2

			else:
				print i+', '+list_zero			

concat_feature(str(sys.argv[1]) , str(sys.argv[2]), str(sys.argv[3]))

