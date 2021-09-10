import csv, gzip, sys

read_length_min = int(sys.argv[1])
read_length_max = int(sys.argv[2])
pc_fasta = sys.argv[3]
project_name = sys.argv[4]
samples = sys.argv[5:]

def parse_header(header):
    info = {}
    header = header.replace('>', '').split('|')
    info['transcript_id'] = header[0]
    info['gene_id'] = header[1]
    for field in header[7:]:
        if 'CDS' in field:
            info['CDS_bounds'] = [int(x) for x in field.replace('CDS:', '').split('-')]
    info['CDS_bounds'][0] -= 1
    if 'UTR3' not in str(header):
        info['CDS_bounds'][1] = None
    return info

def codon_dict(CDS_seq, offset):
    pos_dict = {};codon_list = []
    for i in range(0, int(len(CDS_seq)/3)):
        this_codon = CDS_seq[i*3:i*3+3]
        codon_list.append(this_codon)
        for i in range(i*3,i*3+3):
            pos_dict[i+offset] = this_codon
    return pos_dict, codon_list

codon_lists = {}
b=open(pc_fasta, 'r')
#>ENSMUST00000070533.4|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|AC157543.1-001|Xkr4|3634|UTR5:1-150|CDS:151-2094|UTR3:2095-3634
#AAGGAAAGAGGATAACACTTGAAATGTAAATAAAGAAAATACCTAATAAAAATAAATAAA
seq = ''; i=0;transcript_codons = {}
for line in b:
    line = line.replace('\n', '')
    if line[0] == '>' and i !=0:
        these_bounds = info['CDS_bounds']
        if these_bounds[1] == None:
            CDS = seq[info['CDS_bounds'][0]:]
            info['CDS_bounds'][1] = len(seq)
        else:
            CDS = seq[info['CDS_bounds'][0]: info['CDS_bounds'][1]]
        
        transcript_codons[transcript_id], codon_lists[transcript_id] = codon_dict(CDS, info['CDS_bounds'][0])
        seq = ''
        transcript_id = line.replace('>', '').split('|')[0]
        info = parse_header(line.replace('>', ''))
    elif i ==0:
        i=1
        transcript_id = line.replace('>', '').split('|')[0]
        info = parse_header(line.replace('>', ''))
    else:
        seq +=line
these_bounds = info['CDS_bounds']
if these_bounds[1] == None:
    CDS = seq[info['CDS_bounds'][0]:]
    info['CDS_bounds'][1] = len(seq)
else:
    CDS = seq[info['CDS_bounds'][0]: info['CDS_bounds'][1]]
transcript_codons[transcript_id], codon_lists[transcript_id]= codon_dict(CDS, info['CDS_bounds'][0])
b.close()

acceptable = set(range(read_length_min, read_length_max + 1))
codon_window = range(-15, 16)
for sample in samples:
    observed = {}; i=0;expected = {};these_reads = 0 
    b=gzip.open('riboWaltz/' + project_name + '/' + sample + '_psites_sample.tsv.gz', 'rt')
    reader =csv.DictReader(b, delimiter = '\t') 
    #transcript end5    psite   end3    length  cds_start   cds_stop    psite_from_start    psite_from_stop psite_region    p_site_codon    a_site_codon    e_site_codon    frame
    #ENSMUST00000000001.4    26  37  51  26  142 1206    -105    -1169   5utr    CGG AGA TGA 1
    for line in reader:
        if i == 0:
            prev_transcript = line['transcript']
            i=1
        read_length = int(line['length'])
        if line['psite'] != 'NA':
            psite = int(line['psite'])
            transcript_id = line['transcript']
            #expected
            if read_length in acceptable:
                if prev_transcript != transcript_id and prev_transcript in codon_lists:
                    this_seq = codon_lists[prev_transcript]
                    seq_len = len(this_seq)
                    read_increment = float(these_reads)/seq_len
                    j = 0
                    for codon in this_seq:
                        for relative_location in codon_window:
                            if relative_location + j >= 0 and relative_location + j < seq_len:
                                this_codon = this_seq[relative_location + j]
                                if this_codon not in expected:
                                    expected[this_codon] = {}
                                if relative_location not in expected[this_codon]:
                                    expected[this_codon][relative_location] = 0
                                expected[this_codon][relative_location] += read_increment
                        j+=1
                    these_reads = 0
                if line['psite_region'] == 'cds' and transcript_id in transcript_codons:
                    these_reads +=1
                    for relative_location in codon_window:
                        relative_psite = psite + relative_location*3
                        if relative_psite in transcript_codons[transcript_id]:
                            this_codon = transcript_codons[transcript_id][relative_psite]
                            if this_codon not in observed:
                                observed[this_codon] = {}
                            if relative_location not in observed[this_codon]:
                                observed[this_codon][relative_location] =0
                            observed[this_codon][relative_location] +=1
                prev_transcript = transcript_id
    b.close()
    if transcript_id in this_seq:
        this_seq = codon_lists[transcript_id]
        seq_len = len(this_seq)
        read_increment = float(these_reads)/seq_len
        j = 0
        for codon in this_seq:
            for relative_location in codon_window:
                if relative_location + j >= 0 and relative_location + j < seq_len:
                    this_codon = this_seq[relative_location + j]
                    if this_codon not in expected:
                        expected[this_codon] = {}
                    if relative_location not in expected[this_codon]:
                        expected[this_codon][relative_location] = 0
                    expected[this_codon][relative_location] += read_increment
            j+=1
    values = list(range(-17, -3)) + ['E', 'P', 'A'] + list(range(1,15))
    x_dict = dict(zip(list(range(-15,16)), values))
    c = open('obs_exp/' + project_name + '/' + sample + '_obs_exp_sample.tsv', 'w')
    writer = csv.writer(c, delimiter = '\t')
    writer.writerow(['codon', 'relative_location', 'observed', 'expected', 'ratio'])
    for codon in sorted(expected):
        for relative_loc in sorted(expected[codon]):
            if relative_loc not in observed[codon]:
                observed[codon][relative_loc] = 0
            writer.writerow([codon, x_dict[relative_loc], observed[codon][relative_loc], expected[codon][relative_loc], float(observed[codon][relative_loc])/expected[codon][relative_loc]])
    c.close()


