library(devtools)
library(data.table)
library(riboWaltz)

args <- commandArgs(trailingOnly = TRUE)

gtf_file <- args[1]
fasta_file <- args[2]
offset_tsv <- args[3]
sample_bam <- args[4]

project_name <- unlist(strsplit(sample_bam, "/"))[2]
sample <- unlist(strsplit(sample_bam, "/"))[4]

annotation <- create_annotation(gtf_file, dataSource = "gencode", organism = "Mus musculus")

reads_list <- bamtolist(bamfolder=paste0('bam/', project_name, '/temp/', sample), annotation = annotation, refseq_sep = '|')

filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 10:45)
merged_offsets <- data.frame(read.table(offset_tsv, sep = '\t', header = TRUE))
merged_offsets['sample'] <- sample

rownames(merged_offsets) <- NULL
psite_offset <- data.table(merged_offsets)

reads_psite_list <- psite_info(filtered_list, psite_offset, fastapath = fasta_file, gtfpath = gtf_file, dataSource = "gencode")

print(sample)
output_table <- reads_psite_list[[sample]]

periodicity_frames_stratified <- frame_psite_length(duplicates_filter(data = reads_psite_list, extremity = 'both'), sample = sample, region = "all", length_range = "all")
write.table(periodicity_frames_stratified$dt, paste0('riboWaltz/', project_name, '/', sample, "_periodicity_lengths_deduped.tsv"), sep = '\t', row.names = FALSE, quote = FALSE)
periodicity_frames_stratified <- frame_psite_length(reads_psite_list, sample = sample, region = "all", length_range = "all")
write.table(periodicity_frames_stratified$dt, paste0('riboWaltz/', project_name, '/', sample, "_periodicity_lengths_all.tsv"), sep = '\t', row.names = FALSE, quote = FALSE)

gzip_out <- gzfile(paste0('riboWaltz/', project_name, '/', sample, "_psites_sample.tsv.gz"), "w")
write.table(output_table, gzip_out, sep = '\t', row.names = FALSE, quote = FALSE)
close(gzip_out)


