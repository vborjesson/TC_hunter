#!/usr/bin/python

'''
Part of TC-hunter
Description: Takes an sam file with softclipped reads (created by runSoftClipExtraction.sh) as input and extract 
reads with map quality equal to or above 60. The script creates an txt file containing candidate breakpoint positions.   
Date: 2019-01-25
Author: Vanja Boerjesson
Usage: python karyotype.py --links links.txt 
'''

import sys
import os 
import argparse
import subprocess

########################################### ARGPARSER #############################################

usage = '''karyptype.py takes a txt-file with possible insertion sites generated by FindLinks.py, and returns a karyptype.txt file''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--links', dest='links', help = 'links.txt', required=True)
parser.add_argument('--construct_length', dest='length', help = 'Length of the construct', required=True)
parser.add_argument('--threshold', dest='thres', help = 'The number of links that most exist for one region to be reported', default=  2, required=False)

args = parser.parse_args()

links = args.links
length = args.length
threshold = int(args.thres)

# Create genome list 

genome = ['x', 'y', 'X', 'Y']
for i in range (1,24):
	genome.append(str(i))

#print(genome)

##########################################    Create karyotype    ##########################################

def create_karyotype (links, length):

	karyo_dict = dict()
	construct_dict = dict()

	with open (links, 'r') as f_in, open ('karyotype.txt', 'w') as f_out:
		
		# Extract unique chromosomes and positions 
		for line in f_in:

			list_line = line.split(' ')
			#print(list_line)
			chr1 = list_line[0]
			chr2 = list_line[3]
			pos1 = list_line[1]
			pos2 = list_line[4]

			# see if chr1 or 2 are the construct	
			if chr1 in genome: 
				#Add chr and bp position to dict
				if chr1 in karyo_dict:
					karyo_dict[chr1].append(pos1)
				else:
					karyo_dict[chr1] = [pos1]
				
				#Add construct bp to dict	
				if chr2 in construct_dict:
					construct_dict[chr2].append(pos2)
				else:
					construct_dict[chr2] = [pos2]	

			else: 
				#Add chr and bp position to dict
				if chr2 in karyo_dict:
					karyo_dict[chr2].append(pos2)
				else:
					karyo_dict[chr2] = [pos2]

				#Add construct bp to dict	
				if chr1 in construct_dict:
					construct_dict[chr1].append(pos1)
				else:
					construct_dict[chr1] = [pos1]

		# number of chrosomes: 
		n_karyo = len (karyo_dict)	
	
		# counter = 0
		for chrom in karyo_dict:

			if len(karyo_dict[chrom]) < threshold:
				continue

			bp_max = int(max(karyo_dict[chrom])) + 5000
			bp_min = int(min(karyo_dict[chrom])) - 5000
			chr_name = '{}{}'.format('chr', str(chrom))

			f_out.write('{} {} {} {} {} {} {}\n'.format('chr', '-', chrom, chrom, bp_min, bp_max, "255,128,0"))

		# Print construct to file
		
		for chrom in construct_dict:
			construct_name = chrom	
			f_out.write('{} {} {} {} {} {} {}\n'.format('chr', '-', construct_name, construct_name, 0, length, construct_name))				

		return (karyo_dict, construct_dict)		


create_karyotype(links, length)	


















