#!/usr/bin/env nextflow

workingdirectory = params.workingDir
construct_file = file(params.construct_file)


sequences = Channel
                .fromPath(params.bam)
                .map { file -> tuple(file.baseName, file) }  


//----------------------Samtools sort and index-------------------------------

process samtools_index {
	publishDir workingdirectory, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	module 'samtools/0.1.19'

	input:
		set ID, bam from sequences 

	output:
		set ID, "${ID}_indexed.bam", "${ID}_indexed.bam.bai" into bwa_mem_out	

	script:

	"""  
		ln -s $bam ${ID}_indexed.bam	
		samtools index ${ID}_indexed.bam
	"""	
	
}


// outputs can only be used once as input in a new process, therefor we copy them into several identical outputs. 
bwa_mem_out.into {
  	bwa_mem_out_extractReads
	bwa_mem_out_links
	bwa_mem_out_hist
	bwa_mem_out_plots
}


//----------------------Extract reads with soft clips in construct = BWA mem-------------------------------

process extract_reads_bwa {
	publishDir workingdirectory, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	module 'samtools/1.9'

	input:
		set ID, bam, bai from bwa_mem_out_extractReads

	output:
		set ID, "${ID}_softclipped.sam" into softclipped_out	

	script:	

	"""
		bash ${params.tc_hunter_path}/Scripts/runSoftClipExtraction.sh ${bam} ${ID}_softclipped.sam ${params.construct_name}	
	"""	
	
}



//----------------------Extract reads with supplementary alignments in construct-------------------------------

process create_links_sup {
	publishDir workingdirectory, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	module 'samtools/1.9'

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
		set ID, sam from softclipped_out

	output:
		set ID, "${ID}_links.txt" into links_out

	script:
	"""
		python ${params.tc_hunter_path}/Scripts/FindLinks.py --sam ${sam} --mapq ${params.mapq}
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
		python ${params.tc_hunter_path}/Scripts/createKaryotype.py --links ${links} --construct_length ${params.construct_length}  
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

	module 'samtools/1.9'

	input:
		set ID, file(bam), file(bai), karyo from combined_bam_karyotype 
		//set ID, file(bam), file(bai) from bwa_mem_out_hist		
		//file "${ID}_karyotype.txt" from karyotype_out_hist
		//set file(karyo_file), ID, file(bam), file(bai) from karyotype_out_hist

	output:
		set ID, "${ID}_hist.txt" into hist_out	

	script:

	"""
		python ${params.tc_hunter_path}/Scripts/createHistogram.py --karyo ${karyo} --bam ${bam} 
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

process create_plots {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	module "igv"
	module "R/3.5.1"

	input:

		set ID, bam, bai, karyo, links, hist, sup_links from combined_all
		//set name, links from links_out_circos
		//set ID, bam, bai from bwa_mem_out_plots
		//file "${ID}_links.txt" from links_plot
		//file "${ID}_karyotype.txt" from karyotype_out_circos
		//file "${ID}_hist.txt" from hist_out 
		//file "${ID}_sup_links.txt" from sup_links
		//file jointRef from bwa_mem_out_ref

	output:
		file "${ID}_output.html" into plots_out  

	script:
	"""	
		cp ${links} .
		cp ${karyo} .
		cp ${hist} .
		cp ${sup_links} .
		python ${params.tc_hunter_path}/Scripts/createOutput.py --hist ${hist} --links ${links} --sup_links ${sup_links} --karyo ${karyo} --construct $construct_file --WorkDir ${params.workingDir} --tchunter ${params.tc_hunter_path} --bam ${bam} --ref ${params.reference} --name ${ID}
		cp *pdf ${params.workingDir} || :
		cp *png ${params.workingDir} || :
	"""				
}





//----------------------Cretae final html summary report -------------------------------

process create_html {
	publishDir params.workingDir, mode: 'copy', overwrite: true
	errorStrategy 'ignore'

	input:
		file html from plots_out.collect()

	output:
		file 'output_summary.html' into html_out  

	script:
	"""
		cat ${params.tc_hunter_path}/template/header.html > output_summary.html
		cat ${html} >> output_summary.html
		cat ${params.tc_hunter_path}/template/tail.txt >> output_summary.html
	"""				
}



