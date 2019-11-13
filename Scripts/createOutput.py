#!/usr/bin/env python3

'''
Part of TC-hunter
Description: this script takes links, karyotype and histogram files created so far as input. Statistics is performed in order 
to rank the possible insertin sites and rank them. This script outputs file with possible sites together with two 
figures (pdf).    
Date: 2019-01-29
Author: Vanja Boerjesson
Usage: python createOutput.py --hist hist.txt --links links.txt --sup_links sup_links.txt --karyo karyotype.txt --construct construct.csv 


test run on simulated data.
	
#!/bin/bash -ue
python /jumbo/WorkingDir/B19-001/TC_hunter/Scripts/createOutput.py --hist hist.txt --links links.txt --sup_links 
sup_links.txt --karyo karyotype.txt --construct /jumbo/WorkingDir/B19-001/Intermediate/tc_wgsim/construct.txt 
--WorkDir /jumbo/WorkingDir/B19-001/Intermediate/tc_wgsim_igv --tchunter /jumbo/WorkingDir/B19-001/TC_hunter 
--bam /jumbo/WorkingDir/B19-001/Intermediate/wgsim/PPM1D_M42_BWA_sorted.bam.bam 
--ref /jumbo/WorkingDir/B19-001/Data/Raw/Ref/JointRefGenome.fasta



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
			
		# If only one insertion site	
		if len(site_list) == 2:
			subprocess.call ('cp ' + str(kar) + ' tmp_karyotype.txt', shell = True)
			first_line = True
			with open (kar, 'r') as f_in:
				for line in f_in:
					print (line)
					if first_line:
						print ('first line', line)

						R_karyo = line.split(' ')

						#print(line[0], line[1], line[2])
						
						biomart_string = '{}:{}:{}'.format(R_karyo[0], R_karyo[1], R_karyo[2]) # string with pos to search for genes 
						igv_position_1 = '{}:{}-{}'.format(R_karyo[0], int(R_karyo[1])+4900, int(R_karyo[2])-4900) # igv span
						igv_position_2 = '{}:{}-{}'.format(R_karyo[0], int(R_karyo[1])+4600, int(R_karyo[2])-4600) # igv span zoomed

						html_host_bp = '{} - {}'.format(int(R_karyo[1])+5000, int(R_karyo[2])-5000) # bp position to be reported
						
						first_line = False
						continue
							
					else: # Construct line

						R_karyo_construct = line.split(' ')

						const_bp1, const_bp2 = find_construct_bp(links, R_karyo[0], int(R_karyo[1])+5000, int(R_karyo[2])-5000)
						html_cons_bp = '{} - {}'.format(const_bp1, const_bp2) # bp position to be reported
							

				Gene_annotation_R(biomart_string)
				out_name_R = str('circlize.pdf')
				out_name_igv_1 = str ('igv_zoom.png')
				out_name_igv_2 = str ('igv.png')
				out_name_html = str('best_output.html')
				makePlot_R(out_name_R)
				makePlot_igv(igv_position_1, out_name_igv_1)
				makePlot_igv(igv_position_2, out_name_igv_2)
				makeHTML_output(R_karyo[0], html_host_bp, R_karyo_construct[0], html_cons_bp, out_name_R, out_name_igv_1, out_name_igv_2, out_name_html) # Create output file 


		else: 
			link_df = pd.read_csv(links, sep=' ', names=['chrom1', 'start1','start2', 'chrom2', 'end1', 'end2'])
			sup_link_df = pd.read_csv(sup_links, sep=' ', names=['chrom1', 'start1','start2', 'chrom2', 'end1', 'end2'])

			karyo_df_all = pd.read_csv(kar, sep=' ', names=['chrom', 'start', 'end'])

			# Remove last row = the construct 
			karyo_df = karyo_df_all[:-1]
			construct_df = karyo_df_all[-1:]

			# Add new column for score 
			karyo_df = karyo_df.assign(score = '')

			# Convert sup_link_df chromosomes to string instead of integer
			sup_link_df['chrom1'] = sup_link_df['chrom1'].astype(str)
			sup_link_df['chrom2'] = sup_link_df['chrom2'].astype(str)
			
			#print(karyo_df)
			
			# count which site is the most suported
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
				igv_position_1 = '{}:{}-{}'.format(R_karyo[0], int(R_karyo[1])+4950, int(R_karyo[2])-4950) # long region
				igv_position_2 = '{}:{}-{}'.format(R_karyo[0], int(R_karyo[1])+4600, int(R_karyo[2])-4600) # short region (4600 orig)

				html_host_bp = '{} - {}'.format(R_karyo[1], R_karyo[2])
				
				const_bp1, const_bp2 = find_construct_bp(links, R_karyo[0], int(R_karyo[1])+5000, int(R_karyo[2])-5000)
				html_cons_bp = '{} - {}'.format(const_bp1, const_bp2) # bp position to be reported
				#html_cons_bp = '{} - {}'.format(R_karyo_construct[1], R_karyo_construct[2])

				out_name_R = '{}{}'.format((i+1),'_circlize.pdf') 
				out_name_igv_1 = '{}{}'.format((i+1),'_igv_zoom.png') 
				out_name_igv_2 = '{}{}'.format((i+1),'_igv.png') 
				out_name_html = '{}{}'.format((i+1),'_output.html')
				
				makePlot_R(out_name_R) # circlize plot

				makePlot_igv(igv_position_1, out_name_igv_1) # igv plot
				makePlot_igv(igv_position_2, out_name_igv_2) # zoomed igv plot

				print('igv plot are finished!')
				print('\n')
				print('Create output files')

				makeHTML_output(R_karyo[0], html_host_bp, R_karyo_construct[0], html_cons_bp, out_name_R, out_name_igv_1, out_name_igv_2, out_name_html) # Create output file 

			#print (karyo_df)	
			#print (link_df)
			#print (sup_link_df)

 

			# Look for regions 

########################################### create gene annotation function ###########################################
def Gene_annotation_R (position):
	print('Rscript ' + TC + '/Scripts/biomart.R ' + position)
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

def makeHTML_output (host_chr, host_bp, cons_chr, cons_bp, circ_name, igv1_name, igv2_name, html_name):
	#subprocess.call ('cp ' + os.path.join(TC,"template/output.html") + ' ' + WD, shell=True)
	#subprocess.call ('mv ' + WD + '/output.html ' + WD + '/' + html_name, shell=True)
	output_file = '{}/{}'.format(WD, html_name)
	with open(os.path.join(TC,"template/output.html"), 'r') as myfile:
		template=myfile.read()

		template= template.replace("£££££", "{}\t{}".format(host_chr, host_bp)) # host bp position
		template= template.replace("ˇˇˇˇˇ", "{}\t{}".format(cons_chr, cons_bp)) # construct bp position

		template= template.replace("ªªªªª", circ_name) # Circlize.png
		template= template.replace("ßßßßß", igv1_name) # igv.png
		template= template.replace("®®®®®", igv2_name) # igv_zoom.png

		f= open(output_file, "w")
		f.write(template)
		f.close()

# What we need for output HTML; template path WD/Template/TC_hunter_out.html
# 1) rank number 2) bp position 3) plot_nams * 2


########################################### Find construct bp ###########################################

def find_construct_bp (links, host_chr, host_bp1, host_bp2):
	with open (links, 'r') as links_in:
		
		construct_bp1 = []
		construct_bp2 = []

		for line in links_in:
			line = line.split(' ')

			print(line)
			print(str(host_chr), str(host_bp1), str(host_bp2))

			if (line[0] == host_chr and str(line[1]) == str(host_bp1)):
				print('hej')
				construct_bp1.append(line[4])
				print(construct_bp1)

			elif line[3] == host_chr and str(line[4]) == str(host_bp1):	
				construct_bp1.append(line[1])

			elif line[0] == host_chr and str(line[1]) == str(host_bp2):
				construct_bp2.append(line[4])

			elif line[3] == host_chr and str(line[4]) == str(host_bp2):
				construct_bp2.append(line[1])
	
			#print(max(set(numbers), key=numbers.count))	

		#print(construct_bp1)
		bp_1 = most_frequent(construct_bp1)
		#print(construct_bp2)
		bp_2 = most_frequent(construct_bp2)

	return (bp_1, bp_2)


########################################### Find most frequent number bp ###########################################

def most_frequent(List): 
    counter = 0
    num = List[0] 
      
    for i in List: 
        curr_frequency = List.count(i) 
        if(curr_frequency> counter): 
            counter = curr_frequency 
            num = i 
  
    return num 



########################################### run functions #############################################


n_sites = rank_sites(karyo, links, sup_links)



