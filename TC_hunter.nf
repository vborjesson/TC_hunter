#!/usr/bin/env nextflow

workingdirectory = params.workingDir
construct_file = file(params.construct_file)

//sequences = Channel
//                .fromPath(params.bam)
//                .map { file -> tuple(file.baseName, file) }


//----------------------Extract reads with soft clips in construct-------------------------------

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
		cp ${params.dry_run_softclipped} softclipped.sam 	
	"""		

	}else{

	"""
		bash ${params.tc_hunter_path}/Scripts/runSoftClipExtraction.sh ${params.bam} softclipped.sam 	
	"""	
	}
}


//----------------------Extract reads with supplementary alignments in construct-------------------------------

process create_links_sup {
	publishDir workingdirectory, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	module 'samtools/1.9'

	output:
		file "sup_links.txt" into sup_links	

	script:
		"""
		bash ${params.tc_hunter_path}/Scripts/ExtractConstruct.sh ${params.bam} ${params.construct_name}
		python ${params.tc_hunter_path}/Scripts/ExtractConstruct.py ${params.construct_name}.sam sup_links.txt
		"""

}




//----------------------Extract positions and create links.txt-------------------------------

process create_links_soft {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	input:
		file sam from softclipped_out

	output:
		file 'links.txt' into links_out

	script:
	"""
		python ${params.tc_hunter_path}/Scripts/FindLinks.py --sam ${sam} --mapq ${params.mapq} 
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

	module 'samtools/1.9'

	input:
		file karyo_file from karyotype_out_hist

	output:
		file 'hist.txt' into hist_out	

	script:

	"""
		python ${params.tc_hunter_path}/Scripts/createHistogram.py --karyo ${karyo_file} --bam ${params.bam} 
	 	
	"""	

}

//----------------------Run Outputs-------------------------------

process create_plots {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	module "R/3.5.1"
	module "igv"

	input:
		file 'links.txt' from links_out_circos
		file 'karyotype.txt' from karyotype_out_circos
		file 'hist.txt' from hist_out 
		file "sup_links.txt" from sup_links

	output:
		file '*.pdf' into circos_out 

	script:
	"""
		python ${params.tc_hunter_path}/Scripts/createOutput.py --hist hist.txt --links links.txt --sup_links sup_links.txt --karyo karyotype.txt --construct $construct_file --WorkDir ${params.workingDir} --tchunter ${params.tc_hunter_path} --bam ${params.bam} --ref ${params.reference}
	"""				

}