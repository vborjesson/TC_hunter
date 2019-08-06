#!/usr/bin/python

import sys 
import pandas as pd
import subprocess 
import os
os.path.dirname(os.path.abspath(__file__))
#import pysam

# run:  	 python ../Code/ExtractConstruct.py Wip1.sam test.bam Wip1

reads = str(sys.argv[1]) # Wip1.sam
construct_name = str(sys.argv[3]) # construct name that will be used as output name to txt-file
bam = str(sys.argv[2]) # bam file to extract reads from
#bam_file = pysam.AlignmentFile(bam, "rb")
#new_sam = pysam.AlignmentFile("allpaired.sam", "wb", template=bam_file)
output_name = '{}{}'.format(construct_name, '_IDlist.txt') 
running_dir = os.getcwd()
output_sam = '{}{}'.format(construct_name, '_sup_reads.sam')

def uniq_ID (reads):
	df = pd.read_csv(reads, sep= '\t', header = None)
	with open (output_name, 'w') as f_out:
		id_list = list(df[0])
		#print(id_list) 
		uniq_id = []
		for item in id_list:
			item = str(item)
			#print(item)
			ex = id_list.count(item)
			#print (item, ex)
			if ex == 1:
				uniq_id.append(item)
				f_out.write('{}{}'.format(item, '\n'))
	#print(uniq_id)
	return uniq_id


def create_sam (ids, bam_file, output_sam):
	with open(output_sam, 'w') as outfile:
		subprocess.call('samtools view ' + bam_file + ' |grep -f ' + output_name, shell=True, stdout=outfile)


#def create_sam (ids, bam_file):
#	for read in bam_file.fetch():
#		#print (read.query_name)
#		if str(read.query_name) in ids:
#			new_sam.write(read)			
#			#print(read, 'is found in listand printed to file')

#	bam_file.close()
#	new_sam.close()

print('Find all IDs with reads mapped to construct and host genome')
ids = uniq_ID(reads)
n_id = len(ids)
print(n_id , 'sequences were found with supplementary alignment in host genome')
print('Create SAM-file based on this ids')	
create_sam(ids, bam, output_sam)
print ('new SAM-file created with reads mapped to both host genome and construct')

#print(ids)
print ('------------------------------------------')
print ('Done!')
print ('------------------------------------------')
