# PhotoModGO or PhotoMod2
Multi-label classifier for predicting photosynthetic protein via genome neighborhood feature<br>

========= install ==============<br>
- download scripts from github (PhotoMod2-master)
- download tools from https://drive.google.com/file/d/1q3C3i8EDdJ7B-HOejMd7l5z8DRhFAECe/view?usp=sharing
- extract tools and put them in the same folder of PhotoMod2-master 
- make conda environment using file photomod2.yml

=========== run ==============<br>
python run.py [input fasta] [blast evalue] [cpu core] [ram(GB)] [blast hit] <br>
(e.g. python run.py test_protein.fasta 1e-25 8 2 1)

=========== citation ==============<br>
PhotoModPlus: A webserver for photosynthetic protein prediction from a genome neighborhood feature
Apiwat Sangphukieo, Teeraphan Laomettachit, Marasri Ruengjitchatchawalya
bioRxiv 2020.05.10.087635; doi: https://doi.org/10.1101/2020.05.10.087635
