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
--WorkDir /jumbo/WorkingDir/B19-001/Intermediate/wgsim_2insite/tc_hunter --tchunter /jumbo/WorkingDir/B19-001/TC_hunter 
--bam /jumbo/WorkingDir/B19-001/Intermediate/wgsim_2insite/tc_hunter/BWA_sorted.bam 
--ref /jumbo/WorkingDir/B19-001/Intermediate/wgsim_2insite/tc_hunter/JointRefGenome.fasta



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

						html_host_bp = '{} - {}'.format(int(R_karyo[1])+5000+1, int(R_karyo[2])-5000) # bp position to be reported
						
						first_line = False
						continue
							
					else: # Construct line

						R_karyo_construct = line.split(' ')

						const_bp1, const_bp2 = find_construct_bp(links, R_karyo[0], int(R_karyo[1])+5000+1, int(R_karyo[2])-5000)
						html_cons_bp = '{} - {}'.format(const_bp1, const_bp2) # bp position to be reported
							

				Gene_annotation_R(biomart_string)
				out_name_R = str('circlize.pdf')
				out_name_igv_1 = str ('igv_zoom.png')
				out_name_igv_2 = str ('igv.png')
				out_name_html = str('best_output.html')
				makePlot_R(out_name_R)
				#makePlot_igv(igv_position_1, out_name_igv_1)
				#makePlot_igv(igv_position_2, out_name_igv_2)
				makeHTML_pre(R_karyo[0], html_host_bp, R_karyo_construct[0], html_cons_bp, out_name_R, out_name_igv_1, out_name_igv_2, out_name_html) # Create output file 
				makeHTML(out_name_html)


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

				total_score = score_softlink * 10
				total_score += score_suplink
				karyo_df.loc[i, 'score'] = total_score

			print(karyo_df)

			# sort dataframe based on score 
			karyo_df = karyo_df.sort_values(by = ['score'], ascending=False)
			#print(karyo_df)
			# remove score row 
			karyo_sorted = karyo_df.drop(columns= 'score')
			print (karyo_sorted)

			for i in range(len(karyo_sorted)):
				# create karyotype files for plotting
				R_karyo_construct = construct_df.iloc[0]
				R_karyo = karyo_sorted.iloc[i, :]
				R_karyo_string = '{} {} {}'.format(R_karyo[0], R_karyo[1], R_karyo[2])
				R_karyo_string_construct = '{} {} {}'.format(R_karyo_construct[0], R_karyo_construct[1], R_karyo_construct[2]) 

				# Make file
				subprocess.call('echo chr start end > tmp_karyotype.txt', shell = True)
				subprocess.call('echo ' + R_karyo_string + ' >> tmp_karyotype.txt', shell = True)  
				subprocess.call('echo ' + R_karyo_string_construct + ' >> tmp_karyotype.txt', shell = True)

				# make position string for biomart gene annotation
				biomart_karyo_string = '{}:{}:{}'.format(R_karyo[0], R_karyo[1], R_karyo[2]) 
				Gene_annotation_R(biomart_karyo_string)

				# Make position string for igv 
				igv_position_1 = '{}:{}-{}'.format(R_karyo[0], int(R_karyo[1])+4950, int(R_karyo[2])-4950) # long region
				igv_position_2 = '{}:{}-{}'.format(R_karyo[0], int(R_karyo[1])+4600, int(R_karyo[2])-4600) # short region (4600 orig)

				html_host_bp = '{} - {}'.format(R_karyo[1]+5001, R_karyo[2]-5000)

				# Find construct break point				
				const_bp1, const_bp2 = find_construct_bp(links, R_karyo[0], int(R_karyo[1])+5001, int(R_karyo[2])-5000)
				html_cons_bp = '{} - {}'.format(const_bp1, const_bp2) # bp position to be reported
				#html_cons_bp = '{} - {}'.format(R_karyo_construct[1], R_karyo_construct[2])

				out_name_R = '{}{}'.format((i+1),'_circlize.pdf') 
				out_name_igv_1 = '{}{}'.format((i+1),'_igv_zoom.png') 
				out_name_igv_2 = '{}{}'.format((i+1),'_igv.png') 
				out_name_html = '{}{}'.format((i+1),'_output.html')
				
				makePlot_R(out_name_R) # circlize plot

				#makePlot_igv(igv_position_1, out_name_igv_1) # igv plot
				#makePlot_igv(igv_position_2, out_name_igv_2) # zoomed igv plot

				print('igv plot are finished!')
				print('\n')
				print('Create output files')

				makeHTML_pre(R_karyo[0], html_host_bp, R_karyo_construct[0], html_cons_bp, out_name_R, out_name_igv_1, out_name_igv_2, out_name_html, str(i+1)) # Create output file 

	makeHTML('TC_hunter_report.html')

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

def makeHTML_pre (host_chr, host_bp, cons_chr, cons_bp, circ_name, igv1_name, igv2_name, html_name, rank_n):
	#subprocess.call ('cp ' + os.path.join(TC,"template/output.html") + ' ' + WD, shell=True)
	#subprocess.call ('mv ' + WD + '/output.html ' + WD + '/' + html_name, shell=True)

	subprocess.call ('echo "<hr size=4 align=left color=darkgrey>" >> ranking.txt', shell=True)
	subprocess.call ('echo "<h1>' + rank_n + ' </h1>" >> ranking.txt', shell=True)
	
	subprocess.call ('echo "<hr size=4 align=left color=darkgrey>" >> bp_position.txt', shell=True)	
	subprocess.call ('echo "<h1>' + host_chr + ' ' + host_bp + '</h1>" >> bp_position.txt', shell=True)

	subprocess.call ('echo "<hr size=4 align=left color=darkgrey>" >> figures.txt', shell=True)
	subprocess.call ('echo "<a href=' + circ_name + '><img src=' + circ_name + ' title="Circlize figure" width=40 height=40 /></a> &nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;" >> figures.txt', shell=True)
	subprocess.call ('echo "<a href=' + igv1_name + '><img src=' + igv1_name + ' title="IGV figure" width=40 height=40 /></a> &nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;" >> figures.txt', shell=True)
	subprocess.call ('echo "<a href=' + igv2_name + '><img src=' + igv2_name + ' title="IGV figure" width=40 height=40 /></a> &nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;" >> figures.txt', shell=True)


def makeHTML (name):
	output_file = '{}/{}'.format(WD, name) 
	template1 = '{}/{}'.format(TC, "template/out_part1.txt")
	template2 = '{}/{}'.format(TC, "template/out_part2.txt")

	subprocess.call ('cat ' + template1 + ' > ' + output_file, shell=True)
	subprocess.call ('cat ranking.txt >> ' + output_file, shell=True)
	subprocess.call ('echo "</nav>" >> ' + output_file, shell=True)
	subprocess.call ('echo "<article>" >> ' + output_file, shell=True)
	subprocess.call ('echo "<h2>Predicted breakpoint position in host genome </h2>" >> ' + output_file, shell=True)
	subprocess.call ('cat bp_position.txt >> ' + output_file, shell=True)
	subprocess.call ('echo "</article>" >> ' + output_file, shell=True)
	subprocess.call ('echo "<figures>" >> ' + output_file, shell=True)
	subprocess.call ('echo "<h2>figures</h2>" >> ' + output_file, shell=True)
	subprocess.call ('cat figures.txt >> ' + output_file, shell=True)
	subprocess.call ('echo "</figures>" >> ' + output_file, shell=True)
	subprocess.call ('cat ' + template2 + ' >> ' + output_file, shell=True)

	subprocess.call ('rm ranking.txt', shell=True)
	subprocess.call ('rm bp_position.txt', shell=True)
	subprocess.call ('rm figures.txt', shell=True)

#	output_file = '{}/{}'.format(WD, html_name)
#	with open(os.path.join(TC,"template/output.html"), 'r') as myfile:
#		template=myfile.read()#

#		template= template.replace("£££££", "{}\t{}".format(host_chr, host_bp)) # host bp position
#		template= template.replace("ˇˇˇˇˇ", "{}\t{}".format(cons_chr, cons_bp)) # construct bp position#

#		template= template.replace("ªªªªª", circ_name) # Circlize.png
#		template= template.replace("ßßßßß", igv1_name) # igv.png
#		template= template.replace("®®®®®", igv2_name) # igv_zoom.png#

#		f= open(output_file, "w")
#		f.write(template)
#		f.close()

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

			# Check construct bp that matches the host bp
			if (line[0] == host_chr and str(line[1]) == str(host_bp1)):
				print('hej')
				construct_bp1.append(line[4])
				print(construct_bp1)

			if line[3] == host_chr and str(line[4]) == str(host_bp1):	
				print('new')
				construct_bp1.append(line[1])

			if line[0] == host_chr and str(line[1]) == str(host_bp2):
				print('bp2')
				construct_bp2.append(line[4])

			if line[3] == host_chr and str(line[4]) == str(host_bp2):
				print('bp4')
				construct_bp2.append(line[1])
	
			#print(max(set(numbers), key=numbers.count))	
		print(construct_bp1, construct_bp2)	
		#print(construct_bp1)
		if len(construct_bp1)==0:
			bp_1 = 'Unknown'

		elif len(construct_bp1)!=0:	
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



