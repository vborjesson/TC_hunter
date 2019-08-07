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


// outputs can only be used once as input in a new process, therefor we copy them into several identical outputs. 
links_out.into {
  	links_out_karyo
	links_out_circos
}


//----------------------Create karyotype file-------------------------------

process create_karyotype {
	publishDir params.workingDir, mode: 'copy'
	errorStrategy 'ignore'

	input:
		file links from links_out_karyo

	output:
		file 'karyotype.txt' into karyotype_out

	script:
	"""
		python ${params.tc_hunter_path}/Scripts/createKaryotype.py --links ${links} --construct_length ${params.construct_length}  
	"""	
}

karyotype_out.into {
  	karyotype_out_hist
	karyotype_out_circos
}


//----------------------Create histogram file-------------------------------


process create_histogram {
	publishDir params.workingDir, mode: 'copy'
	errorStrategy 'ignore'

	input:
		file karyo_file from karyotype_out_hist

	output:
		file 'hist.txt' into hist_out	

	script:

	"""
		python ${params.tc_hunter_path}/Scripts/createHistogram.py --karyo ${karyo_file} --bam ${params.bam} 
	"""	

}

//----------------------Run Circos-------------------------------

process create_plots {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	module "circos/0.69"

	input:
		file 'links.txt' from links_out_circos
		file 'karyotype.txt' from karyotype_out_hist
		file 'hist.txt' into hist_out 

	output:
		set 'circos.png', 'circos.tif' into circos_out 

	script:
	"""
		circos -conf ${params.tc_hunter_path}/Circos/circos.conf
	"""				

}




