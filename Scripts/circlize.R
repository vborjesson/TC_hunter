#!/usr/bin/env Rscript

# install.packages("circlize")

##########################################

library(circlize)
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)
log(-1) 
#option_list = list(
#  make_option(c("-f", "--links"), type="character", default='links.txt', 
#              help="links file", metavar="character"),
#	make_option(c("-o", "--sup_links"), type="character", default="sup_links.txt", 
#              help="sup_links file", metavar="character")
#  make_option(c("-f", "--karyotype"), type="character", default='tmp_karyotype,txt', 
#              help="karyotype file", metavar="character"),
#	make_option(c("-o", "--genes"), type="character", default="genes.txt", 
#              help="gene annotation file", metavar="character")
#  make_option(c("-f", "--construct"), type="character", default="constract.csv", 
#              help="file with construct info", metavar="character"),
#	make_option(c("-o", "--hist"), type="character", default="hist.txt", 
#              help="histogram file", metavar="character")              	
#); 
 
#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);

# Load data
links <- data.table::fread("links.txt", data.table = F)
sup_links <- data.table::fread("sup_links.txt", data.table = F)
karyo <- data.table::fread("1_tmp_karyotype.txt", data.table = F)
genes <- data.table::fread('genes.csv', data.table = F) 
construct <- data.table::fread('construct.csv', data.table = F)
hist <- data.table::fread("hist.txt", data.table = F)
out <- args[1]
out <- 'out_test2.pdf'

chrom_name = karyo[1,1]
construct_name = karyo[2,1]

hist_split <- split(hist, hist$V1)

hist_chrom <- hist_split[[chrom_name]]
hist_construct <- hist_split[[construct_name]]

# Data wrangling
hist_to_long <- function(d){
	a <- integer()
	b <- integer()
	for (i in 1:nrow(d))
	{
	   a <- c(a, d[i, 2], d[i, 3])
	   b <- c(b, rep(d[i, 4], 2)) 
	}
	return(data.frame(position = a, read_depth = b))
} 

chrom_table <- hist_to_long(hist_chrom)
chrom_table$id <- chrom_name

construct_table <- hist_to_long(hist_construct)
construct_table$id <- construct_name

df <- rbind(chrom_table, construct_table)
#head(df)

max_h <- max(chrom_table$read_depth)
max_h
#circos.initializeWithIdeogram(plotType = NULL)

pdf(out,width=10,height=10,paper='special')

# Initiate plot 
circos.genomicInitialize(karyo, tickLabelsStartFromZero = FALSE, axis.labels.cex = 1, labels.cex = 1.5, major.by = 1000)

# add track 
circos.track(factors = df$id, ylim = c(0, 2), bg.col = c("#BDBDBD", "#8A0808"), track.height=0.02, cell.padding=c(0,0,0,0))
#circos.track(factors = genes$chromosome_name, ylim = c(0, 2), bg.col = c("#D9681D", "#15576D"), track.height=0.02, cell.padding=c(0,0,0,0))

# Initiate construct tracks
circos.track(df$id, ylim=c(0,1), track.height=0.015, cell.padding=c(0,0,0,0), bg.border = 'white')
circos.track(df$id, ylim=c(0,1), track.height=0.015, cell.padding=c(0,0,0,0), bg.border = 'white')

# Plot construct
for (pos in 1:nrow(construct)){
	circos.rect(construct[pos,2], 0, construct[pos,3], 1, sector.index = construct_name, col='black', track.index=3)
	circos.text(((construct[pos,2]+construct[pos,3])/2), 0.5, labels= construct[pos, 1], cex= 1.2, track.index=4, bg.border = 'white')
}

# Plot annotation genes
#circos.track(df$id, ylim=c(0,1), track.height=0.015, cell.padding=c(0,0,0,0), bg.border = 'white')

# restrict genes to be within region 

for (pos in 1:nrow(genes)){
	# restrict genes to be within region  	
	if (genes[pos, 5] < karyo[1,2]){
		genes[pos, 5] <- karyo[1,2]
	}	
	if (genes[pos, 6] > karyo[1,3]){
		genes[pos, 6] <- karyo[1,3]
	}
	# Create gene-ines for annotations 
	circos.lines(c(genes[pos,5], genes[pos,6]), c(0.5, 0.5), sector.index = chrom_name, col='#2E2E2E', track.index=3)
}

a_start_pos <- genes$start_position
a_end_pos <- genes$start_position + 80
a_box <- data.frame(a_start_pos, a_end_pos)
#a_box

a_start_pos <- genes$end_position
a_end_pos <- genes$end_position + 20
b_box <- data.frame(a_start_pos, a_end_pos)
#b_box

box <- rbind(a_box, b_box)
#box

circos.genomicRect(box, ytop=1, ybottom=0, col='black', sector.index=chrom_name, track.index=3)

# Add text to genes 
#for (pos in 1:nrow(genes)){
#	circos.text(((genes[pos,5]+genes[pos,6])/2), 0.5, labels= genes[pos, 2], cex= 1.2, sector.index = chrom_name, track.index=4, bg.border = 'white')
#}

# Histogram - read coverage
circos.track(factors = df$id, ylim = c(0, max_h), bg.col = 'gray96', bg.border = 'white')

circos.lines(construct_table$position, construct_table$read_depth, sector.index = construct_name, type = "s", area = TRUE, col = "#A4A4A4")
circos.lines(chrom_table$position, chrom_table$read_depth, sector.index = chrom_name, type = "s", area = TRUE, col = "#A4A4A4")



# links

for (i in 1:nrow(sup_links)){
	tryCatch({
	#print(sup_links[i,6])
	circos.link(sup_links[i,1], c(sup_links[i,2]-20, sup_links[i,3]+20), sup_links[i,4], c(sup_links[i,5]-20, sup_links[i,6]+20), col = "#B40404")}, 
        error=function(error_message) {
            message("This is my custom message.")
            message("And below is the error message from R:")
            message(error_message)
            return(NA)
        }
    )
}


for (i in 1:nrow(links)){
	#print(sup_links[i,6])
	circos.link(links[i,1], c(links[i,2]-30, links[i,3]+30), links[i,4], c(links[i,5]-30, links[i,6]+30), col = "#1C1C1C")
}

dev.off()
circos.clear()


