library(BiocManager)
library(tidyverse)
library(EnsDb.Mmusculus.v79)

args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
this_sample <- args[2]

edbx <- EnsDb.Mmusculus.v79
pause_df <- read_tsv(paste0('pauses/', project_name, '/pauses.tsv')) 

print(this_sample)
pause_df <- as_tibble(pause_df)

pause_df <- pause_df %>% dplyr::filter(sample == this_sample)
pause_df <- tibble::rownames_to_column(pause_df)

rng_tx <- IRanges(start = pause_df$coordinate, width = 3, names = pause_df$transcript)
split_transcript_names <- strsplit(names(rng_tx), ".", fixed = TRUE)
ensembl_transcripts <- unlist(lapply(split_transcript_names, function(l) l[[1]]))
names(rng_tx) <- ensembl_transcripts
rng_gnm <- transcriptToGenome(rng_tx, edbx)
rng_gnm <- as_tibble(data.frame(seqname=seqnames(rng_gnm), start = start(rng_gnm), end = end(rng_gnm), strand = strand(rng_gnm)))
rng_gnm <- rng_gnm %>% mutate(seqname.value = as.character(seqname.value), strand.value = as.character(strand.value))

genome_plus_df <- rng_gnm %>% dplyr::filter(strand.value == '+') %>% group_by(seqname.group) %>%
    summarize(chrom = paste0('chr',dplyr::first(seqname.value)), spliced_start = min(start.value), spliced_end = max(end.value),  strand = dplyr::first(strand.value))
genome_minus_df <- rng_gnm %>% dplyr::filter(strand.value == '-') %>% group_by(seqname.group) %>% 
    summarize(chrom = paste0('chr', dplyr::first(seqname.value)),  spliced_start = min(end.value), spliced_end = max(start.value), strand = dplyr::first(strand.value))
genome_tib <- bind_rows(list(genome_plus_df, genome_minus_df)) %>% arrange(seqname.group)

genome_tib$seqname.group = as.character(genome_tib$seqname.group)
genome_pauses <- left_join(pause_df, genome_tib, by = c('rowname' = 'seqname.group'))
out_names <- names(genome_pauses)[2:length(names(genome_pauses))]
##sample transcript  coordinate  region  score   psite_density   BG_left BG_right    BG_transcript   P_codon A_codon relative_start  relative_end    chrom   spliced_start   spliced_end strand
##Hbs1l_embryo_rep1   ENSMUST00000000001.4    426 cds 25.4626748874   49  0   0   3.23098591549   GGG GAA 95  -260    chr3    108123560   108123558   -
good_pauses <- genome_pauses[genome_pauses$strand %in% c('-', '+'),]
good_pauses <- good_pauses[,out_names]

good_pauses <- good_pauses %>% arrange(chrom, relative_start)
write_tsv(good_pauses, path = paste0('pauses/', project_name, '/', this_sample, '_genome_pauses_collapsed.tsv.gz'))
#collapsed_genome_tib <- good_pauses[,out_names] %>% group_by(sample, chrom, spliced_start, spliced_end, strand) %>%
#    mutate(transcript = paste0(transcript, collapse= ';'), coordinate = paste0(coordinate, collapse= ';'), region = paste0(unique(region), collapse = ';'), score = paste0(score, collapse = ';'), psite_density = paste0(unique(psite_density), collapse = ';'), BG_left = paste0(unique(BG_left), collapse = ';'),  BG_right = paste0(unique(BG_right), collapse = ';'),  BG_transcript = paste0(BG_transcript, collapse = ';'), P_codon = paste0(unique(P_codon), collapse = ';'), A_codon = paste0(unique(A_codon), collapse = ';'), relative_start = paste0(relative_start, collapse = ';'), relative_end = paste0(relative_end, collapse = ';'), shortmer_fraction = paste0(unique(shortmer_fraction), collapse = ';')) 
#collapsed_genome_tib <- collapsed_genome_tib %>% arrange(chrom, relative_start)
#write_tsv(collapsed_genome_tib, path = paste0(base_dir, this_sample, '_0.5_genome_pauses_collapsed.tsv.gz'))
