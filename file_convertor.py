import commands
import sys,os

input_file=sys.argv[1]
current_path=commands.getoutput('pwd')

print 'Read input_file',input_file
#input file is final_combination3.csv
if os.path.isfile(input_file) == True:
	#use this flow to convert file to nominal and saving as arff
	file_kf=open('./convert2arff.kf','r')
	file_kf_read=file_kf.read()
	file_kf_read=file_kf_read.replace("csv_inputfile",current_path+'/'+input_file)
	file_kf_read=file_kf_read.replace("outputpath",current_path)
	file_kf_read=file_kf_read.replace("outputfilename",'converted_input.arff')
	file_kf.close()

	file_kf_new=open('convert2arffRUN.kf','w')
	file_kf_new.write(file_kf_read)
	file_kf_new.close()
	num_c=1
	while True: #to make sure that converted_input.arff is created
		print 'Check converted_input.arff , round ',num_c
		print commands.getoutput('java -cp ./weka-3-9-1/weka.jar -Xmx24000m -Xss5000m weka.knowledgeflow.FlowRunner convert2arffRUN.kf')
		num_c+=1
		if int(os.path.getsize("converted_input.arff")) > 0 :
			print 'Pass!! : converted_input.arff is created'
			break

	commands.getoutput("sed \"s/@attribute \' label\' {\' photo\'}/@attribute \' label\' {\' photo\',\' non_photo\'}/g\" converted_input.arff > "+sys.argv[2]+"_converted_input.arff")
else:
	print 'No input file =>  ',input_file
