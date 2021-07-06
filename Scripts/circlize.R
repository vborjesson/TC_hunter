#!/usr/bin/env Rscript

if (!require('circlize')) install.packages('circlize'); library('circlize')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('data.table')) install.packages('data.table'); library('data.table')


##########################################

args = commandArgs(trailingOnly=TRUE)
options(warn=-1)
#log(-1) 
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


sample_name <- args[3]


# Load data
links <- data.table::fread(paste(sample_name, "_links.txt", sep=""), data.table = F)
sup_links <- data.table::fread(paste(sample_name, "_sup_links.txt", sep=""), data.table = F)
karyo <- data.table::fread("tmp_karyotype.txt", data.table = F)
genes <- data.table::fread("genes.csv", data.table = F) 
#construct <- data.table::fread('construct.txt', data.table = F)
hist <- data.table::fread(paste(sample_name, "_hist.txt", sep=""), data.table = F)
out <- args[1]
construct_path <- args[2]
construct <- data.table::fread(construct_path, data.table = F)

sample_name = 'ins2_100'
out = '1_ins2_100_circlize.pdf'
construct_path = '/jumbo/WorkingDir/B19-001/Intermediate/validation/sim_data/fly_2ins_new/construct.txt'


####################################Only for test#####################
#sample_name <- 'M411_BWA_sorted'
#out <- 'circ_test.pdf'
#out <- 'out.pdf'
#construct_path = '/jumbo/WorkingDir/B19-001/Intermediate/FINAL_M42/construct.txt'
#sample_name = "sample1"

#links <- data.table::fread("PPM1D_M45_links.txt", data.table = F)
#sup_links <- data.table::fread("PPM1D_M45_sup_links.txt", data.table = F)
#karyo <- data.table::fread("tmp_karyotype.txt", data.table = F)
#genes <- data.table::fread("genes.csv", data.table = F) 
#construct <- data.table::fread('/jumbo/WorkingDir/B19-001/Intermediate/FINAL_M41/construct.txt', data.table = F)
#hist <- data.table::fread("PPM1D_M45_hist.txt", data.table = F)
#out <- '3_PPM1D_M45_circlize.pdf'
#construct_path <- '/jumbo/WorkingDir/B19-001/Intermediate/FINAL_M41/construct.txt'
#construct <- data.table::fread(construct_path, data.table = F)
#sample_name <- 'PPM1D_M45'
####################################Only for test#####################


chrom_name = karyo[1,1]
construct_pos = length(karyo$V1)
construct_name = karyo[2,1]

hist_split <- split(hist, hist$V1)

hist_chrom <- hist_split[[chrom_name]]
hist_construct <- hist_split[[construct_name]]

#------- Function: covert hist tab-df to circlize friendly format----
hist_to_long2 <- function(d){
	position <- as.matrix(d[,2:3])
	position <- as.vector(t(position))
	read_depth <- as.matrix(d[, 4])
	read_depth <- cbind(read_depth, read_depth)
	read_depth <- as.vector(t(read_depth))
	
	return(data.frame(position = position, read_depth = read_depth))	
}

# Run hist to long for chromosome
chrom_table_dup <- hist_to_long2(hist_chrom)
chrom_table <- chrom_table_dup %>% distinct()
chrom_table$id <- chrom_name

#print (chrom_table)

# Try to shorten the dataframe if it's more than > 5000 rows 
space <- nrow(chrom_table) / 5000
#print (dim(chrom_table))



space <- round(space, digits = -1)


if (space > 10){
	chrom_table_short <- chrom_table %>% filter(row_number() %% space == 1)
	chrom_table <- chrom_table_short
} else if (space > 0) {
	chrom_table_short <- chrom_table %>% filter(row_number() %% space == 1)
	chrom_table <- chrom_table_short
} else {
	chrom_table <- chrom_table
}



# Run hist to long for construct
construct_table_dup <- hist_to_long2(hist_construct)
construct_table <- construct_table_dup %>% distinct()
construct_table$id <- construct_name


df <- rbind(chrom_table, construct_table)

#head(chrom_table)
#head(df)
#max_h <- max(chrom_table$read_depth)
max_h <- max(df$read_depth)
#print(max_h)

#circos.initializeWithIdeogram(plotType = NULL)


# Calculate what major (tick space) to use for labeling
start_chrom = max(chrom_table[1])
end_chrom = min(chrom_table[1])
label_major = (start_chrom - end_chrom) / 10
label_major = round(label_major, digits = -3)



#-----------------------------------------------------------------------------------
############################### Plot circle plot ###################################
#-----------------------------------------------------------------------------------

pdf(out,width=10,height=10,paper='special')
#png(out, width=1000, height=1000, res=300)

# Initiate plot 
circos.genomicInitialize(karyo, tickLabelsStartFromZero = FALSE, axis.labels.cex = 1, labels.cex = 1.5, major.by = label_major, sector.width=1)
# add track 
circos.track(factors = df$id, ylim = c(0, 2), bg.col = c("#BDBDBD", "#8A0808"), track.height=0.02, cell.padding=c(0,0,0,0))
#circos.track(factors = genes$chromosome_name, ylim = c(0, 2), bg.col = c("#D9681D", "#15576D"), track.height=0.02, cell.padding=c(0,0,0,0))



# Initiate construct tracks
circos.track(df$id, ylim=c(0,1), track.height=0.015, cell.padding=c(0,0,0,0), bg.border = 'white')
circos.track(df$id, ylim=c(0,1), track.height=0.015, cell.padding=c(0,0,0,0), bg.border = 'white')



# Plot construct
for (pos in 1:nrow(construct)){
	print(pos)
	circos.rect(construct[pos,2], 0, construct[pos,3]-30, 1, sector.index = construct_name, col='black', track.index=3)
	circos.text(((construct[pos,2]+construct[pos,3])/2), 0.5, labels= construct[pos, 1], cex= 1.2, track.index=4, bg.border = 'white', facing= 'bending.inside')
}


# Plot annotation genes
#circos.track(df$id, ylim=c(0,1), track.height=0.015, cell.padding=c(0,0,0,0), bg.border = 'white')

# restrict genes to be within region 
if (nrow(genes) > 1){
	for (pos in 1:nrow(genes)){
		# restrict genes to be within region  	
		if (genes[pos, 5] < karyo[1,2]){
			genes[pos, 5] <- karyo[1,2]
		}	
		if (genes[pos, 6] > karyo[1,3]){
			genes[pos, 6] <- karyo[1,3]
		}
		# Create gene-lines for annotations 
		circos.lines(c(genes[pos,5], genes[pos,6]), c(0.5, 0.5), sector.index = chrom_name, col='#2E2E2E', track.index=3)
		#circos.lines(c(23250000, 23270000), c(0.5, 0.5), sector.index = chrom_name, col='#2E2E2E', track.index=3)
	}


	a_start_pos <- genes$start_position
	a_end_pos <- genes$start_position + 80
	a_box <- data.frame(a_start_pos, a_end_pos)
	#a_box <- data.frame(23250000, 23250000+80)

	#a_box

	a_start_pos <- genes$end_position
	a_end_pos <- genes$end_position + 20
	b_box <- data.frame(a_start_pos, a_end_pos)
	#b_box <- data.frame(23257000-20, 23257000)
	#b_box

	box <- rbind(a_box, b_box)
	#box

	circos.genomicRect(box, ytop=1, ybottom=0, col='black', sector.index=chrom_name, track.index=3)

	# Add text to genes 
	for (pos in 1:nrow(genes)){
		circos.text(((genes[pos,5]+genes[pos,6])/2), 0.5, labels= genes[pos, 2], cex= 1.2, sector.index = chrom_name, track.index=4, bg.border = 'white', facing= 'bending.inside')
	}

}



## Add text to genes 
#for (pos in 1:nrow(genes)){
#	circos.text(((genes[pos,5]+genes[pos,6])/2), 0.5, labels= genes[pos, 2], cex= 1.2, sector.index = chrom_name, track.index=4, bg.border = 'white', facing= 'bending.inside')
#}

# Histogram - read coverage
circos.track(factors = df$id, ylim = c(0, max_h), bg.col = 'gray96', bg.border = 'white')



circos.lines(construct_table$position, construct_table$read_depth, sector.index = construct_name, type = "s", area = TRUE, col = "#A4A4A4")
circos.lines(chrom_table$position, chrom_table$read_depth, sector.index = chrom_name, type = "s", area = TRUE, col = "#A4A4A4")



# links
#head(sup_links)
#head(links)

for (i in 1:nrow(sup_links)){
	tryCatch({
	#print(sup_links[i,6])
	circos.link(sup_links[i,1], c(sup_links[i,2]-20, sup_links[i,3]+20), sup_links[i,4], c(sup_links[i,5]-20, sup_links[i,6]+20), col = "#B40404")}, 
        error=function(error_message) {
            #message("This is my custom message.")
            #message("And below is the error message from R:")
            message(error_message)
            return(NA)
        }
    )
}

for (i in 1:nrow(links)){
	tryCatch({
	#print(sup_links[i,6])
	circos.link(links[i,1], c(links[i,2]-30, links[i,3]+30), links[i,4], c(links[i,5]-30, links[i,6]+30), col = "#1C1C1C")},
        error=function(error_message) {
            #message("This is my custom message.")
            #message("And below is the error message from R:")
            message(error_message)
            return(NA)
		}            
	)	
}


text(-1,-1, "-", adj = 0, col = 'red', cex=3)
text(-0.92,-1, "Anchor read", adj = 0)
text(-1,-1.05, "-", adj = 0, cex=3)
text(-0.92,-1.05, "Soft clipped (chimeric) read", adj = 0)



dev.off()
circos.clear() 



#-----------------------------------------------------------------------------------
############################### Plot line plot ###################################
#-----------------------------------------------------------------------------------
#

#head(chrom_table)#

#library(ggplot2)
#pdf('test_lineplot.pdf',width=20,height=5,paper='special')
#p2 <- ggplot() + geom_line(aes(x=position, y= read_depth), data = chrom_table, stat="identity") 
#p2 + geom_rect(aes(xmin = 69401000, ymin = -Inf, 
#                 xmax = 69401500, ymax = Inf),
#             fill = "steelblue") 
#p2
#dev.off()#
#

#pdf('test_lineplot.pdf',width=20,height=5,paper='special')
#p2 <- ggplot()+
#	geom_rect(aes(xmin = 69401000, ymin = -Inf, 
#                 xmax = 69401500, ymax = Inf),
#             fill = "steelblue", alpha = 0.5) +	
#	geom_line(aes(x=position, y= read_depth), data = chrom_table, stat="identity") 
#	theme_bw()
#p2
#dev.off()




