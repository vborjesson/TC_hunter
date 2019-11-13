#!/usr/bin/R

# Create list of genes found within region 
# usage: Rscript biomart.R chrom:pos_start:pos_end

args = commandArgs(trailingOnly=TRUE)
loc <- args[1] # region to look for genes within

library(biomaRt)

#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

all.genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"), filters = "chromosomal_region", values = loc, mart = ensembl)

write.csv(all.genes, file= 'genes.csv')

