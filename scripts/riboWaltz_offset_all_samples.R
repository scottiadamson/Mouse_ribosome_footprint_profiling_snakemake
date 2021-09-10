library(devtools)
library(data.table)
library(riboWaltz)

args <- commandArgs(trailingOnly = TRUE)

gtf_file <- args[1]
annotation <- create_annotation(gtf_file, dataSource = "gencode", organism = "Mus musculus")
fasta_file <- args[2]
merged_bam <- args[3]
project_name <- unlist(strsplit(merged_bam, "/"))[2]
bam_dir <- paste0(paste(unlist(strsplit(merged_bam, "/"))[1:4], collapse = '/'), '/')
sample <- unlist(strsplit(bam_dir, '/'))[length(unlist(strsplit(bam_dir, '/')))]

print('Processing merged bam')
reads_list <- bamtolist(bamfolder = bam_dir, annotation = annotation,  refseq_sep = '|')
print('Filtering read lengths')
filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 10:45)
filtered_list <- duplicates_filter(data = filtered_list, extremity = "both")
print('Identifying P-sites')
psite_offset <- psite(filtered_list, flanking = 6)
write.table(psite_offset, paste0('riboWaltz/', project_name, '/offset.tsv'), row.names=FALSE, quote =FALSE, sep = '\t')
#write.table(psite_offset, paste0('riboWaltz/', project_name, '/', sample, '_offset.tsv'), row.names=FALSE, quote =FALSE, sep = '\t')
print('Successful operation')
