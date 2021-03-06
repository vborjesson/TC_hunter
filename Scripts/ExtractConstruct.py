#!/usr/bin/python

import sys 
import pandas as pd
import os

# run:  	 python ../Code/ExtractConstruct.py Wip1.sam 
sam = sys.argv[1]
out = sys.argv[2]

# Check if sam file is empty, and create an Empty file
if os.path.getsize(sam) == 0:
	f = open(out, "w")
	f.write("")
	f.close()	


#Create a dataframe with data 
elif os.path.getsize(sam) != 0: 
	df = pd.read_csv(sam, sep= "\t", header = None)
	
	# Change column name 
	df = df.rename(columns={0:'ID', 5:'cigar', 2:'chr1', 3:'pos1', 6:'chr2', 7:'pos2'})

	#print (df) 
	# rm rows with softclips (these are treated in FindLinks instead)
	df = df[~df.cigar.str.contains('S')]
	df = df[~df.cigar.str.contains('H')]

	# drop all column that we don't need 
	df_keep = df[['chr1', 'pos1', 'pos1', 'chr2', 'pos2', 'pos2']]
	df_keep.to_csv(out, index=False, header=False, sep=' ')



