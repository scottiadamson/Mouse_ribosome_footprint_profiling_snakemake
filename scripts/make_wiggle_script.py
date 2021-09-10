import csv, sys
from collections import defaultdict

project_name = sys.argv[1]
offset_file = sys.argv[2]
read_length_min = int(sys.argv[3])
read_length_max = int(sys.argv[4])

#length  total_percentage        start_percentage        around_start    offset_from_5   offset_from_3   corrected_offset_from_5 corrected_offset_from_3 sample
#27      8.7     11.4    T       11      15      11      15      July_2017_young_old_merged
a = open(offset_file, 'r')
reader = csv.DictReader(a, delimiter = '\t')
offset_dict = defaultdict(set)
for line in reader:
    if int(line['length']) >= read_length_min and int(line['length']) <= read_length_max:
        offset_dict[line['corrected_offset_from_5']].add(int(line['length']))
a.close()

b = open('scripts/make_wiggles_' + project_name + '.sh', 'w')
b.write('bam_dir="bam/' + project_name + '"\n')
b.write('chrom_sizes="references/mm10.chrom.sizes"\n')
b.write('base_dir="wiggles/' + project_name + '"\n')
b.write('sample=$1\n')
b.write('samtools view -H $bam_dir/$sample"_genome_sorted.bam" > $base_dir/$sample"_header.sam"\n')
for offset in offset_dict:
    awk_conditions = []
    for read_length in offset_dict[offset]:
        awk_conditions.append('length($10) == ' + str(read_length))
    b.write('samtools view $bam_dir/$sample"_genome_sorted.bam" | awk ' + "'(" + '||'.join(awk_conditions) + ")' - | cat $base_dir/$sample" + '"_header.sam" - | samtools view -bh - > $base_dir/$sample"_' + str(offset) + '.bam"\n')
    b.write('samtools index $base_dir/$sample"_' + str(offset) + '.bam"\n')
    b.write('bamCoverage --Offset ' + str(offset) + ' --bam $base_dir/$sample"_' + str(offset) + '.bam" --outFileName $base_dir/$sample"_' + str(offset) + '.bed" --outFileFormat bedgraph --binSize 1\n')
    b.write('read_total_' + str(offset) + '=$(samtools view $base_dir/$sample"_' + str(offset) + '.bam" -c )\n')
    b.write('rm $base_dir/$sample"_' + str(offset) + '.bam"\n')
    b.write('rm $base_dir/$sample"_' + str(offset) + '.bam.bai"\n')
b.write('let read_totals=' + '+'.join(['read_total_' + str(offset) for offset in list(offset_dict)]) + '\n')
b.write('bedtools unionbedg -i $base_dir/$sample*.bed | awk ' + "'{sum=$4+$5+$6; print $1,$2,$3,sum}' - | sort -k1,1 -k2,2n - | python scripts/normalize_wiggles.py $read_totals - > $base_dir/$sample" + '"_normed.bed"\n')
b.write('scripts/bedGraphToBigWig $base_dir/$sample"_normed.bed" $chrom_sizes wiggles/' + project_name + '/$sample".bw"\n')
b.write('rm $base_dir/$sample*.bed\n')
b.write('rm $base_dir/$sample"_header.sam"')
b.close()

