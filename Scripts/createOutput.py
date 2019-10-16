#!/usr/bin/env python3

'''
Part of TC-hunter
Description: this script takes links, karyotype and histogram files created so far as input. Statistics is performed in order 
to rank the possible insertin sites and rank them. This script outputs file with possible sites together with two 
figures (pdf).    
Date: 2019-01-29
Author: Vanja Boerjesson
Usage: python createOutput.py --hist hist.txt --links links.txt --sup_links sup_links.txt --karyo karyotype.txt --construct construct.csv 
'''

import sys
import os 
import argparse
import subprocess
import pandas as pd

########################################### ARGPARSER #############################################

usage = '''This script takes links, karyotype and histogram files created so far as input. Statistics is performed in order 
to rank the possible insertin sites and rank them. This script outputs file with possible sites together with two 
figures (pdf)''' 

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('--hist', dest='hist', help = 'hist.txt', required=True)
parser.add_argument('--links', dest='links', help = 'links.txt', required=True)
parser.add_argument('--sup_links', dest='sup_links', help = 'sup_links.txt', required=True)
parser.add_argument('--karyo', dest='karyo', help = 'karyotype.txt', required=True)
parser.add_argument('--construct', dest='construct', help = 'construct.txt', required=True)
parser.add_argument('--WorkDir', dest='WD', help = 'WorkingDir', required=True)
parser.add_argument('--tchunter', dest='TC', help = 'Path to TC_hunter', required=True)
parser.add_argument('--bam', dest='bam', help = 'path to bam file', required=True)
parser.add_argument('--ref', dest='ref', help = 'reference fasta file for host', required=True)

args = parser.parse_args()

hist = args.hist
links = args.links
sup_links = args.sup_links
karyo = args.karyo
construct = args.construct
WD = args.WD
TC = args.TC
bam = args.bam
ref = args.ref

########################################### ranking function ###########################################

# Read in file if more then 1 insertion site: count for how many supporting links they have and rank 
def rank_sites (kar, links, sup_links):

	with open (kar, 'r') as f_in:
		site_list = []
		for line in f_in:
			line = line.split(' ')
			site_list.append(line[0])
			

		if len(site_list) == 2:
			subprocess.call ('cat ' + str(kar) + ' > tmp_karyotype.txt', shell = True)
			first_line = True
			with open (kar, 'r') as f_in:
				for line in f_in:
					if first_line:
						print('this is the first line')
						line = line.split(' ')
						biomart_string = '{}:{}:{}'.format(line[0], line[1], line[2])
						igv_position_1 = '{}:{}-{}'.format(line[0], int(line[1])+4900, int(line[2])-4900)
						igv_position_2 = '{}:{}-{}'.format(line[0], int(line[1])+4600, int(line[2])-4600)
						first_line == False
						break	

				Gene_annotation_R(biomart_string)
				out_name_R = str('circlize.pdf')
				out_name_igv_1 = str ('igv_zoom.png')
				out_name_igv_2 = str ('igv.png')
				makePlot_R(out_name_R)
				makePlot_igv(igv_position_1, out_name_igv_1)
				makePlot_igv(igv_position_2, out_name_igv_2)

		else: 
			link_df = pd.read_csv(links, sep=' ', names=['chrom1', 'start1','start2', 'chrom2', 'end1', 'end2'])
			sup_link_df = pd.read_csv(sup_links, sep=' ', names=['chrom1', 'start1','start2', 'chrom2', 'end1', 'end2'])

			karyo_df_all = pd.read_csv(kar, sep=' ', names=['chrom', 'start', 'end'])

			# Remove last row = the construct 
			karyo_df = karyo_df_all[:-1]
			construct_df = karyo_df_all[-1:]

			# Add new csolumn for score 
			karyo_df = karyo_df.assign(score = '')

			# Convert sup_link_df chromosomes to string instead of integer
			sup_link_df['chrom1'] = sup_link_df['chrom1'].astype(str)
			sup_link_df['chrom2'] = sup_link_df['chrom2'].astype(str)
			
			#print(karyo_df)
			
			# count witch site is the most suported
			for i in range(len(karyo_df)):
				site = karyo_df['chrom'].iloc[i]
				
				score_softlink = int((link_df['chrom1'] == site).sum())
				score_softlink += int((link_df['chrom2'] == site).sum())
				
				score_suplink = int((sup_link_df['chrom1'] == site).sum())
				score_suplink += int((sup_link_df['chrom2'] == site).sum())
				#print (site, score_softlink, score_suplink)		

				total_score = score_softlink * 2
				total_score += score_suplink
				karyo_df.loc[i, 'score'] = total_score

			# sort dataframe based on score 
			karyo_df.sort_values(by = ['score'])
			# remove score row 
			karyo_sorted = karyo_df.drop(columns= 'score')
			#print (karyo_sorted)

			for i in range(len(karyo_sorted)):
				# create karyotype files for plotting
				R_karyo_construct = construct_df.iloc[0]
				R_karyo = karyo_sorted.iloc[i, :]
				R_karyo_string = '{} {} {}'.format(R_karyo[0], R_karyo[1], R_karyo[2])
				R_karyo_string_construct = '{} {} {}'.format(R_karyo_construct[0], R_karyo_construct[1], R_karyo_construct[2]) 

				# Make file
				subprocess.call('echo chr start end > ' + WD + '/tmp_karyotype.txt', shell = True)
				subprocess.call('echo ' + R_karyo_string + ' >> ' + WD + '/tmp_karyotype.txt', shell = True)  
				subprocess.call('echo ' + R_karyo_string_construct + ' >> ' + WD + '/tmp_karyotype.txt', shell = True)

				# make position string for biomart gene annotation
				biomart_karyo_string = '{}:{}:{}'.format(R_karyo[0], R_karyo[1], R_karyo[2]) 
				Gene_annotation_R(biomart_karyo_string)

				# Make position string for igv 
				igv_position_1 = '{}:{}-{}'.format(R_karyo[0], int(R_karyo[1])+4950, int(R_karyo[2])-4950) 
				igv_position_2 = '{}:{}-{}'.format(R_karyo[0], int(R_karyo[1])+4600, int(R_karyo[2])-4600) 
				out_name_R = '{}{}'.format((i+1),'_circlize.pdf') 
				out_name_igv_1 = '{}{}'.format((i+1),'_igv_zoom.png') 
				out_name_igv_2 = '{}{}'.format((i+1),'_igv.png') 
				makePlot_R(out_name_R)

				makePlot_igv(igv_position_1, out_name_igv_1)
				makePlot_igv(igv_position_2, out_name_igv_2)

			#print (karyo_df)	
			#print (link_df)
			#print (sup_link_df)

 

			# Look for regions 

########################################### create gene annotation function ###########################################
def Gene_annotation_R (position):
	print('Rscrip ' + TC + '/Scripts/biomart.R ' + position)
	subprocess.call('Rscript ' + TC + '/Scripts/biomart.R ' + position, shell = True)



########################################### Plot function ###########################################
def makePlot_R (out_name):
	print('Rscript ' + TC + '/Scripts/circlize.R ', out_name + ' ' + construct)
	subprocess.call('Rscript ' + TC + '/Scripts/circlize.R ' + out_name + ' ' + construct, shell = True)


#def makePlot_python ():	

########################################### IGV function ###########################################
def makePlot_igv (position, out_name):

	# Create batch file to use as input in igv 
	subprocess.call ('echo new > ' + WD + '/igv.bat', shell = True)
	subprocess.call ('echo load ' + bam + ' >> ' + WD + '/igv.bat', shell= True)
	subprocess.call ('echo snapshotDirectory ' + WD +  ' >> ' + WD + '/igv.bat', shell = True)
	subprocess.call ('echo genome ' + ref +  ' >> ' + WD + '/igv.bat', shell = True)
	subprocess.call ('echo goto ' + position +  ' >> ' + WD + '/igv.bat', shell = True)
	subprocess.call ('echo sort base >> ' + WD + '/igv.bat', shell = True)
	subprocess.call ('echo collapse  >> ' + WD + '/igv.bat', shell = True)
	subprocess.call ('echo snapshot ' + out_name + ' >> ' + WD + '/igv.bat', shell = True)	
	subprocess.call ('echo exit >> ' + WD + '/igv.bat', shell = True)

	print('igv.bat file is done!')
	print('igv.sh -b ' + WD + '/igv.bat')

	#subprocess.call ('igv.sh -b ' + WD + '/igv.bat', shell = True)

	subprocess.call ('igv.sh -b ' + WD + '/igv.bat', shell = True)


########################################### Create html output ###########################################



########################################### run functions #############################################


n_sites = rank_sites(karyo, links, sup_links)



