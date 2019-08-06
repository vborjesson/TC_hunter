#!/usr/bin/python

'''
Part of TC-hunter
Description: Takes an sam file with softclipped reads (created by runSoftClipExtraction.sh) as input and extract 
reads with map quality above 40. The script creates an txt file containing candidate breakpoint positions.   
Date: 2019-01-22
Author: Vanja Boerjesson
Usage: python cand_insite.py --sam 
'''

import sys
import os 
import argparse
import subprocess

########################################### ARGPARSER #############################################

usage = '''cand_insite.py takes a sam file created by runSoftClipExtraction.sh as input and returns a txt-file with possible insertion sites''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--sam', dest='sam', help = 'Add path to sam-file created by runSoftClipExtraction.sh', required=True)
parser.add_argument('--mapq', dest='q', help = 'Threshold to use for map quality filtering', default=60, required=False)
#parser.add_argument('--output_name', dest='o', help = 'name of output file', required=True)


args = parser.parse_args()

sam = args.sam
out = 'links.txt'
q = int(args.q)


########################################## Calculate breakpoint position ###################################


def find_bp (pos, cig):

	for i in cig:
		if i == 'S':
			return pos
			break

		if i == 'M':	
			n_M = cig.split('M')
			bp = int(pos) + int(n_M[0])
			return bp


########################################## Filter based on Quality and create txt-file ###################################

def create_txt (sam_file):
	print ('Filter sam-file based on map quality and calculate breakpoint positions')

	with open (sam_file, 'r') as f_in, open (out, 'w') as f_out:
		for line in f_in:
			sam_tab = line.split('\t')

			# sort out reads having map quality score below 40
			if int(sam_tab[4]) < q:
				continue 

			sa = sam_tab[14]	
			sa_sep_com = sa.split(',') 
			
			if int(sa_sep_com[4]) < q:
				continue
			sa_sep_col = sa_sep_com[0].split(':')
			#print (sa_sep_com[4])
			chrom1 = sam_tab[2]
			chrom2 = sa_sep_col[2]

			# Mapping position 1 and 2
			map_pos1 = sam_tab[3]
			map_pos2 = sa_sep_com[1]

			# Cigar for both pos
			cigar1 = sam_tab[5]
			cigar2 = sa_sep_com[3]

			bp1 = find_bp(map_pos1, cigar1)
			bp2 = find_bp(map_pos2, cigar2)

			#print (chrom1, chrom2, bp1, bp2)
			#print ('{}\t{}\t{}\t{}'.format(chrom1, bp1, chrom2, bp2))
			f_out.write('{} {} {} {} {} {}\n'.format(chrom1, bp1, bp1, chrom2, bp2, bp2))

			#print line 


create_txt(sam)




























