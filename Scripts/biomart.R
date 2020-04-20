#!/usr/bin/R

# Create list of genes found within region 
# usage: Rscript biomart.R chrom:pos_start:pos_end

args = commandArgs(trailingOnly=TRUE)
loc <- args[1] # region to look for genes within

#host <- 'mouse'

#library(biomart)

#print('start ')#

##human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#if (host == 'human') {	
#	try(ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "useast"))
#	try(all.genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"), filters = "chromosomal_region", values = loc, mart = ensembl))#

#}  else if (host == 'mouse') {
#	try(ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", mirror = "useast"))
#	try(all.genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"), filters = "chromosomal_region", values = loc, mart = ensembl))#

#} else {
#	host == 'Other' 
#	all.genes <- ''
#}


#ensembl = userMart("ensembl") 
#ensembl = useDataset("mmusculus_gene_ensembl",mart = ensembl)

all.genes = ''

#all.genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"), filters = "chromosomal_region", values = loc, mart = ensembl)
#all.genes
write.csv(all.genes, file= 'genes.csv')


