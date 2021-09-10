import csv, sys

unformatted_counts_file = sys.argv[1]
a = open(unformatted_counts_file, 'r')
reader =csv.reader(a , delimiter = '\t')
#ENSMUSG00000102693.1    chr1    3073253 3074322 +   1070    0   0   0   0   0   0   0   0   0   0   0   0   0   0
#Geneid Chr Start End Strand Length samples...
next(reader)
header = next(reader)
i=0
for item in header:
    header[i] = item.split('/')[-1].replace('_genome_sorted.bam', '')
    i+=1

samples = header[6:]

b = open(unformatted_counts_file.replace('.txt', '_formatted.txt'), 'w')

writer = csv.writer(b , delimiter = '\t')
writer.writerow(['gene_id'] + samples)
for line in reader:
    line = dict(zip(header, line))
    new_line = [line['Geneid']]
    for sample in samples:
        new_line.append(line[sample])
    writer.writerow(new_line)
a.close();b.close()





