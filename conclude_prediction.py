import commands
import os
def conclude_result():
	totalBlastHits=commands.getoutput("grep '' matched2DB.query -c")
	count_photo=int(commands.getoutput("tail -n+6 prediction.csv |cut -d',' -f3|grep '1:photo' -c"))
	count_nonphoto=int(commands.getoutput("tail -n+6 prediction.csv |cut -d',' -f3|grep '2:non_photo' -c"))

	total_predict=count_photo+count_nonphoto

	percent_photo=(float(count_photo)/total_predict)*100
	percent_nonphoto=(float(count_nonphoto)/total_predict)*100
	prot=commands.getoutput("head -n1 matched2DB.query |cut -f1")

	predict='non_photo'
	if percent_photo >=50.0:
		predict='photo'

	if predict=='non_photo': #report only confident score of the prediction class
		ave_prob=commands.getoutput("grep '2:non_photo' prediction.csv|awk -F',' '{ total += $5; count++ } END { print total/count }'")
	else:
		ave_prob=commands.getoutput("grep '1:photo' prediction.csv|awk -F',' '{ total += $5; count++ } END { print total/count }'")


	print os.path.exists('prediction_report.tab')
	if os.path.exists('prediction_report.tab') == True:
		print 'write report...'
		fileout=open('prediction_report.tab','a')
		fileout.write(prot+'\t'+totalBlastHits+'\t'+str(count_photo)+' ('+str(percent_photo)+')\t'+str(count_nonphoto)+' ('+str(percent_nonphoto)+')\t'+str(int(totalBlastHits)-int(total_predict))+'\t'+predict+'\t'+ave_prob+'\n')
		fileout.close()
		
	else:
		print 'Make report file...'
		j=open('prediction_report.tab','w')
		j.write('Query	totalBlastHits	PhotoHits(%)	nonphotoHits(%)	noPredictHits	prediction	Ave_Prob_prediction\n')
		j.write(prot+'\t'+totalBlastHits+'\t'+str(count_photo)+' ('+str(percent_photo)+')\t'+str(count_nonphoto)+' ('+str(percent_nonphoto)+')\t'+str(int(totalBlastHits)-int(total_predict))+'\t'+predict+'\t'+ave_prob+'\n')
		j.close()
	
