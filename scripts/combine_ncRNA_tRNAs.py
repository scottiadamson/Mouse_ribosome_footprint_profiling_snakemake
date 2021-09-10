import gzip

tRNA_lines = ''
a = open('mm10-mature-tRNAs.fa', 'r')
for line in a:
    tRNA_lines += line
a.close()
a = gzip.open('GCF_000001635.27_GRCm39_rna.fna.gz', 'rt')
b = open('references/rRNA_snoRNA_tRNA.fa', 'w')
for line in a:
    if line[0] == '>':
        if 'NR' in line:
            if 'ribosomal RNA' in line and 'non-coding RNA' not in line:
                write = True
            elif 'small nucleolar RNA' in line and 'long non-coding RNA' not in line:
                write = True
            else:
                write = False
        else:
            write = False
    if write:
        b.write(line)
b.write(tRNA_lines)
a.close();b.close()

