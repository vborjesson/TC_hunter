#!/usr/bin/python

'''
Part of TC-hunter
Description: Takes karyotype.txt as input and execute samtools depth for these regions. The output will be saved in a depth.tsv file and 
is used as input to create a hist.txt file.   
Date: 2019-01-29
Author: Vanja Boerjesson
Usage: python createHistogram.py --karyo karyotype.txt 
'''

import sys
import os 
import argparse
import subprocess

########################################### ARGPARSER #############################################

usage = '''Takes karyotype.txt as input and execute samtools depth for these regions. The output will be saved in a depth.tsv file and 
is used as input to create a hist.txt file''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--karyo', dest='karyo', help = 'karyotype.txt', required=True)
parser.add_argument('--bam', dest='bam', help = 'sorted and indexted bam file', required=True)

args = parser.parse_args()

kar = args.karyo
bam = args.bam


########################################### Create hist.txt function #############################################

def create_depth (kar, bam):

	with open (kar, 'r') as f_in, open ('depth.tsv', 'w') as f_out:
		
		script = 'samtools depth -r '
		output = ' >> depth.tsv'

		for line in f_in:
			l = line.split(' ')
			#print(l)
			chrom = l[2]
			start_pos = l[4]
			end_pos = l[5]

			# run samtools to find depth 
			print (script + chrom + ':' + start_pos + '-' + end_pos + ' ' + bam + output)
			subprocess.call(script + chrom + ':' + start_pos + '-' + end_pos + ' ' + bam + output, shell=True)


def create_hist (file):
	with open (file, 'r') as f_in, open ('hist.txt', 'w') as f_out:
		# print list
		l_list = []
		for line in f_in:
			l=line.split('\t')
			# first rond
			if len(l_list) == 0:
				l_list = [l[0], l[1], l[1], l[2]] 

			else: 
				if l[0] == l_list[0] and l[2] == l_list[3]: # Make region larger
					l_list[2] = l[1]
				else:
					f_out.write('{} {} {} {}'.format(l_list[0], l_list[1], l_list[2], l_list[3]))
					l_list = [l[0], l[1], l[1], l[2]]	




# Run functions
create_depth(kar, bam)			
create_hist('depth.tsv')