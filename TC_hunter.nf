#!/usr/bin/env nextflow

params.bam = ''
params.construct_length = ''
params.workingDir = ''
params.tc_hunter_path = ''
params.dry_run = ''
params.dry_run_softclipped = ''

workingdirectory = params.workingDir

//sequences = Channel
//                .fromPath(params.bam)
//                .map { file -> tuple(file.baseName, file) }


//----------------------Extract reds with soft clips in construct-------------------------------

process extract_reads {
	publishDir workingdirectory, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	module 'samtools/1.9'

//	input:
//		set file_ID, file from sequences

	output:
		file "softclipped.sam" into softclipped_out	

	script:

	if (params.dry_run){

	"""
		echo 'hej'
		cp ${params.dry_run_softclipped} softclipped.sam 	
	"""		

	}else{

	"""
		bash ${params.tc_hunter_path}/Scripts/runSoftClipExtraction.sh ${params.bam} softclipped.sam 	
	"""	
	}
}


//----------------------Extract positions and create links.txt-------------------------------

process create_links {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	input:
		file sam from softclipped_out

	output:
		file 'links.txt' into links_out

	script:
	"""
		python ${params.tc_hunter_path}/Scripts/FindLinks.py --sam ${sam}  
	"""	

}



//----------------------Create karyotype file-------------------------------

process create_karyotype {
	publishDir params.workingDir, mode: 'copy'
	errorStrategy 'ignore'

	input:
		file links from links_out

	output:
		file 'karyotype.txt' into karyotype_out

	script:
	"""
		python ${params.tc_hunter_path}/Scripts/createKaryotype.py --links ${links} --construct_length ${params.construct_length}  
	"""	
}


//----------------------Create histogram file-------------------------------


process create_histogram {
	publishDir params.workingDir, mode: 'copy'
	errorStrategy 'ignore'

	input:
		file karyo_file from karyotype_out

	output:
		file 'hist.txt' into hist_out	

	script:

	"""
		python ${params.tc_hunter_path}/Scripts/createHistogram.py --karyo ${karyo_file} --bam ${params.bam} 
	"""	

}

