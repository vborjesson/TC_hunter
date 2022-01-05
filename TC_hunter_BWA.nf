#!/usr/bin/env nextflow

workingdirectory = params.workingDir
construct_file = file(params.construct_file)

params.folder = ""
params.sample = ""


//----------------------------- HELP MESSAGE ------------------------------------

if (params.help) {
    helpMessage()
    exit 0
}

def helpMessage() {
  log.info """

Usage:
It's recomended to use the config file to add all arguments, but you can also add it in the command line. 
nextflow run TC_hunter_BWA.nf -c tc_hunter.config [other Nextflow arguments]

Mandatory arguments
--workingDir 				Working directory path
--tc_hunter_path			TC_hunter path
--construct_file 			construct.txt path
--construct_length 			Length of construct, a number
--construct_name			construct name, must be identical to identifer in construct ref fasta file 
--host_reference 			host reference fasta file
--construct_ref 			construct reference fasta file
(--folder)   				if running several samples: path to folder with sample folders
(--sample)					if running one sample: path to folder with R1.fa + R2.fa

Optional arguments:
--mapq      				mapping quality threshold for soft clipped links [30]
--n_threads 				number of threads to use for mapping [8]
--supporting_reads 			Number of supporting chimeric reads needed to call a TIS [1]
--help                      This usage statement

other optional Nextflow arguments
-resume						Execute the script using the cached results, useful to continue executions that was stopped by an error
-with-report				Create processes execution html report

"""
}


//--------------------- TC-HUNTER MAIN ------------------------------------


// Channel all samples from folder in order to parallelize longranger
if (params.folder) {
	String character = "/*";
	String folder_path = params.folder;
	sample_path =folder_path+character;
	fastq_path = Channel.fromPath(sample_path, type: 'dir').map {path -> tuple(path.name, path)}
}

// If just one sample; no channels are needed

if (params.sample) {
	String character = "/*";
	String folder_path = params.sample;
	sample_path =folder_path;
	//System.out.println(otherString);
	fastq_path = Channel.fromPath(sample_path, type: 'dir').map {path -> tuple(path.name, path)}
}




//----------------------Run BWA mem (mapping)-------------------------------

process bwa_mem {
	publishDir workingdirectory, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

//	module 'bwa/0.7.5a:samtools/0.1.19'

	input:
		set ID, path from fastq_path

	output:
		set ID, "${ID}_indexed.bam", "${ID}_indexed.bam.bai" into bwa_mem_out	
		file "JointRefGenome.fasta" into bwa_mem_out_ref

	script:

	"""
		cat ${params.host_reference} > JointRefGenome.fasta
		cat ${params.construct_ref} >> JointRefGenome.fasta
		bwa index JointRefGenome.fasta
		bwa mem -t ${params.bwa_threads} JointRefGenome.fasta ${path}/*R1* ${path}/*R2* | samtools view -Sb - >  ${ID}.bam	
		samtools sort -o ${ID}_indexed.bam ${ID}.bam
		samtools index ${ID}_indexed.bam
		samtools flagstat ${ID}_indexed.bam > ${params.workingDir}/${ID}_indexed.flagstat
		samtools idxstats ${ID}_indexed.bam > ${params.workingDir}/${ID}_indexed.idxstats		
	"""		
}

// outputs can only be used once as input in a new process, therefor we copy them into several identical outputs. 
bwa_mem_out.into {
  	bwa_mem_out_extractReads
  	bwa_mem_out_soft
	bwa_mem_out_links
	bwa_mem_out_hist
	bwa_mem_out_plots
}


//----------------------Extract reads with soft clips in construct = BWA mem-------------------------------

/***
process extract_reads_bwa {
	publishDir workingdirectory, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

//	module 'samtools/1.9'

	input:
		set ID, bam, bai from bwa_mem_out_extractReads

	output:
		set ID, "${ID}_softclipped.sam" into softclipped_out	

	script:

	"""
		bash ${params.tc_hunter_path}/Scripts/runSoftClipExtraction.sh ${bam} ${ID}_softclipped.sam ${params.construct_name}	
	"""		
}
***/


//----------------------Extract reads with supplementary alignments in construct-------------------------------

process create_links_sup {
	publishDir workingdirectory, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

//	module 'samtools/1.9'

	input:
		set ID, bam, bai from bwa_mem_out_links

	output:
		set ID, "${ID}_sup_links.txt" into sup_links	

	script:
		"""
		bash ${params.tc_hunter_path}/Scripts/ExtractConstruct.sh ${bam} ${params.construct_name}
		python ${params.tc_hunter_path}/Scripts/ExtractConstruct.py ${params.construct_name}.sam ${ID}_sup_links.txt
		"""

}




//----------------------Extract positions and create links.txt-------------------------------

process create_links_soft {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	input:
		set ID, bam, bai from bwa_mem_out_soft

	output:
		set ID, "${ID}_links.txt" into links_out

	script:
	"""
		bash ${params.tc_hunter_path}/Scripts/runSoftClipExtraction.sh ${bam} ${ID}_softclipped.sam ${params.construct_name}
		python ${params.tc_hunter_path}/Scripts/FindLinks.py --sam ${ID}_softclipped.sam --mapq ${params.mapq}
		mv links.txt ${ID}_links.txt 
	"""	

}


// outputs can only be used once as input in a new process, therefor we copy them into several identical outputs. 
links_out.into {
  	links_out_karyo
	links_out_circos
	links_out_plot
}


//----------------------Create karyotype file-------------------------------

process create_karyotype {
	publishDir params.workingDir, mode: 'copy'
	errorStrategy 'ignore'

	input:
		set ID, file(links) from links_out_karyo

	output:
		set ID, file("${ID}_karyotype.txt") into karyotype_out

	script:
	"""
		python ${params.tc_hunter_path}/Scripts/createKaryotype.py --links ${links} --construct_length ${params.construct_length} --construct_name ${params.construct_name} --threshold ${params.supporting_reads}
		mv karyotype.txt ${ID}_karyotype.txt 
	"""	
}

karyotype_out.into {
  	karyotype_out_hist
	karyotype_out_circos
}


// merge karyotype and bam output
//karyotype_out_hist.cross(bwa_mem_out_hist).map{it ->  [it[0][0],it[1][1],it[1][2],it[1][3]]}
//karyotype_out_hist.cross(bwa_mem_out_hist).subscribe { println it }

//----------------------Create histogram file-------------------------------

combined_bam_karyotype = bwa_mem_out_hist.cross(karyotype_out_hist).map{
        it ->  [it[0][0],it[0][1],it[0][2],it[1][1]]
}

process create_histogram {
	publishDir params.workingDir, mode: 'copy'
	errorStrategy 'ignore'

//	module 'samtools/1.9'

	input:
		set ID, file(bam), file(bai), file(karyo) from combined_bam_karyotype 
	
	output:
		set ID, "${ID}_hist.txt" into hist_out	

	script:

	"""
		python ${params.tc_hunter_path}/Scripts/createHistogram.py --karyo ${ID}_karyotype.txt --bam ${bam} 
	 	mv hist.txt ${ID}_hist.txt
	"""	

}

//----------------------Cretae Plots-------------------------------

combined_bam_karyotype2 = bwa_mem_out_plots.cross(karyotype_out_circos).map{
        it ->  [it[0][0],it[0][1],it[0][2],it[1][1]]
}
combined_2 = combined_bam_karyotype2.cross(links_out_plot).map{
        it ->  [it[0][0],it[0][1],it[0][2],it[0][3],it[1][1]]
}
combined_3 = combined_2.cross(hist_out).map{
        it ->  [it[0][0],it[0][1],it[0][2],it[0][3],it[0][4],it[1][1]]
}
combined_all = combined_3.cross(sup_links).map{
        it ->  [it[0][0],it[0][1],it[0][2],it[0][3],it[0][4],it[0][5],it[1][1]]
}


jointref_path = workingdirectory+'/JointRefGenome.fasta'

process create_plots {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

//	module "R/3.5"

	input:
		set ID, bam, bai, karyo, links, hist, sup_links from combined_all

	output:
		file "${ID}_output.html" into plots_out  

	script:
	"""	
		cp ${links} .
		cp ${karyo} .
		cp ${hist} .
		cp ${sup_links} .
		python ${params.tc_hunter_path}/Scripts/createOutput.py --hist ${hist} --links ${links} --sup_links ${sup_links} --karyo ${karyo} --construct $construct_file --WorkDir ${params.workingDir} --tchunter ${params.tc_hunter_path} --bam ${bam} --ref $jointref_path --name ${ID} 
		mkdir ${params.workingDir}/meta_data_${ID} 
		mv ${params.workingDir}/${ID}_* ${params.workingDir}/meta_data_${ID}			
		cp *pdf ${params.workingDir} || :
		cp *png ${params.workingDir} || :
	"""				
}





//----------------------Cretae final html summary report -------------------------------

process create_html {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	input:
		file "*_output.html" from plots_out.collect()

	output:
		file 'output_summary.html' into html_out  

	script:
	"""
		cat ${params.tc_hunter_path}/template/header.html > output_summary.html
		cat *_output.html >> output_summary.html
		cat ${params.tc_hunter_path}/template/tail.txt >> output_summary.html		
	"""				
}








