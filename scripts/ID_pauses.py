import csv, numpy,sys, gzip

#ID_pauses.py reads_percodon_threshold read_length_min read_length_max pc_fasta project_name samples... 
background_threshold = float(sys.argv[1]) 
read_length_min = int(sys.argv[2])
read_length_max = int(sys.argv[3])
pc_fasta = sys.argv[4]
project_name = sys.argv[5]
samples = sys.argv[6:]

print('Background threshold: ' + str(background_threshold))
print('Read Type: ' + str(project_name))
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
    pos_dict = {};codon_list = []; codon_centers = {}
    for i in range(0, int(len(CDS_seq)/3)):
        this_codon = CDS_seq[i*3:i*3+3]
        codon_list.append(this_codon)
        for j in range(i*3,i*3+3):
            pos_dict[j+offset] = this_codon
            codon_centers[j+offset] = i*3+offset
    return pos_dict, codon_list, codon_centers

def update_transcript_info(transcript_id, info, seq):
    these_bounds = info['CDS_bounds']
    if these_bounds[1] == None:
        CDS = seq[info['CDS_bounds'][0]:]
        info['CDS_bounds'][1] = len(seq)
    else:
        CDS = seq[info['CDS_bounds'][0]: info['CDS_bounds'][1]]
    transcript_codons[transcript_id], codon_lists[transcript_id], codon_centers[transcript_id] = codon_dict(CDS, info['CDS_bounds'][0])
    return info

codon_lists = {}
b=open(pc_fasta, 'r')
#>ENSMUST00000070533.4|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|AC157543.1-001|Xkr4|3634|UTR5:1-150|CDS:151-2094|UTR3:2095-3634
#AAGGAAAGAGGATAACACTTGAAATGTAAATAAAGAAAATACCTAATAAAAATAAATAAA
seq = ''; i=0;transcript_codons = {};codon_centers = {}
for line in b:
    line = line.replace('\n', '')
    if line[0] == '>' and i !=0:
        update_transcript_info(transcript_id, info, seq)
        seq = ''
        transcript_id = line.replace('>', '').split('|')[0]
        info = parse_header(line.replace('>', ''))
    elif i ==0:
        i=1
        transcript_id = line.replace('>', '').split('|')[0]
        info = parse_header(line.replace('>', ''))
    else:
        seq +=line
update_transcript_info(transcript_id, info, seq)
b.close()

def conditional_dict(dict_, key):
    if key in dict_:
        return dict_[key]
    else:
        return 0

def center_psites(psite_density, transcript_id, psite_info):
    codon_density = {}; codon_info ={}
    this_center_dict = codon_centers[transcript_id]
    for psite in psite_density:
        region = psite_info[psite]['region']
        if region == 'cds':
            try:
                codon_center = this_center_dict[psite]
                codon_info[codon_center] = psite_info[psite]
                if codon_center not in codon_density:
                    codon_density[codon_center] = 0 
                codon_density[codon_center] += psite_density[psite]
            except KeyError:#at stop codon causes issues
                codon_density[psite] = psite_density[psite]
                codon_info[psite] = psite_info[psite]
        else:
            codon_density[psite] = psite_density[psite]
            codon_info[psite] = psite_info[psite]
    return codon_density, codon_info

def identify_index_flanks(psite, transcript_id, region, codon_density):
    BG_left = [];BG_right = []
    if psite in codon_centers[transcript_id]:
        codon_center = codon_centers[transcript_id][psite]
    else:
        codon_center = psite
    if region == 'cds':
        offset = int(psite) - codon_center
    else:
        offset = 0
    for k in range(-9,0):
        this_loc = k*3 + offset + psite
        BG_left.append(conditional_dict(codon_density, this_loc))
    BG_right = []
    for k in range(1,10):
        this_loc = k*3 + offset +3 + psite
        BG_right.append(conditional_dict(codon_density, this_loc))
    BG_left = numpy.mean(BG_left)
    BG_right = numpy.mean(BG_right)
    return BG_left, BG_right

def psite_add_codons(line, transcript_id, psite):
    psite_info[psite] = dict(zip(['region', 'from_start', 'from_end'], [line[x] for x in ['psite_region', 'psite_from_start', 'psite_from_stop']]))
    if int(float(line['psite_from_start'])) >=0:
        if int(float(line['psite_from_stop'])) <= -3:
            psite_info[psite]['P_codon'] = transcript_codons[transcript_id][psite]
        else:
            psite_info[psite]['P_codon'] = 'NA'
        if int(float(line['psite_from_stop'])) <= -6:
            psite_info[psite]['A_codon'] = transcript_codons[transcript_id][psite+3]
        else:
            psite_info[psite]['A_codon'] = 'NA'
    else:
        psite_info[psite]['A_codon'] = 'NA'; psite_info[psite]['P_codon'] = 'NA'

def calculate_backgound_pauses(transcript):
    this_seq = codon_lists[transcript]
    pos_dict = transcript_codons[transcript]
    BG_transcript = 0
    for psite in psite_density:
        if psite in pos_dict:#only for P-sites in known CDS
            BG_transcript += psite_density[psite]
    BG_transcript = float(BG_transcript)/len(this_seq)#average footprints per codon
    if BG_transcript >= background_threshold:
        if transcript not in pass_threshold:
            pass_threshold[transcript] = set()
        pass_threshold[transcript].add(sample)
        codon_density, codon_info = center_psites(psite_density, transcript, psite_info)
        # for each codon observed, calculate the background on each flank
        for codon in sorted(codon_density):
            pause_count = codon_density[codon]
            BG_left, BG_right = identify_index_flanks(codon, transcript, codon_info[codon]['region'], codon_density) 
            BG = max(BG_left,BG_right,BG_transcript)
            pause_score = (pause_count - BG)/(BG**0.5)
            if pause_score >=10:
                if codon_info[codon]['region'] == 'cds':
                    writer.writerow([sample, transcript, codon, codon_info[codon]['region'], round(pause_score,5), pause_count, round(BG_left,5), round(BG_right, 5), round(BG_transcript,5), codon_info[codon]['P_codon'], codon_info[codon]['A_codon'], int((int(codon_info[codon]['from_start'])-1)/3), int(int(codon_info[codon]['from_end'])/3)])
                else:
                    writer.writerow([sample, transcript, codon, codon_info[codon]['region'], round(pause_score,5), pause_count, round(BG_left,5), round(BG_right,5), round(BG_transcript,5), codon_info[codon]['P_codon'], codon_info[codon]['A_codon'], codon_info[codon]['from_start'], codon_info[codon]['from_end']])

c=open('pauses/' + project_name + '/pauses.tsv', 'w')
writer = csv.writer(c, delimiter = '\t')
writer.writerow(['sample', 'transcript', 'coordinate', 'region', 'score', 'psite_density', 'BG_left', 'BG_right', 'BG_transcript', 'P_codon', 'A_codon', 'relative_start', 'relative_end'])
pauses = {};pass_threshold = {}
read_lengths = set(range(read_length_min, read_length_max + 1))
j=0
for sample in samples:
    print(sample)
    pauses[sample] = {}
    i=0; psite_density = {};psite_info = {}
    b=gzip.open('riboWaltz/' + project_name + '/' + sample + '_psites_sample.tsv.gz', 'rt')
    reader =csv.DictReader(b, delimiter = '\t') 
    #transcript end5    psite   end3    length  cds_start   cds_stop    psite_from_start    psite_from_stop psite_region    p_site_codon    a_site_codon    e_site_codon    frame
    #ENSMUST00000162897.1    6   18  35  30  0   0   0   0   NA  TCC AAC ATC 0
    for line in reader:
        if i == 0:
            prev_transcript = line['transcript']
            i=1
        read_length = int(line['length'])
        if line['psite'] != 'NA':
            psite = int(line['psite'])
            #check read length and CDS capacity and add psite density if so
            if read_length in read_lengths and line['transcript'] in codon_lists:
                transcript_id = line['transcript']
                #if the transcript is a coding transcript, calculate the background distribution
                if prev_transcript != transcript_id and prev_transcript in codon_lists:
                    j+=1
                    calculate_backgound_pauses(prev_transcript)
                    psite_density = {};psite_info = {}
                psite_density[psite] = conditional_dict(psite_density, psite)
                psite_density[psite] += 1
                #update psite info if the psite is new
                if psite not in psite_info:
                    psite_add_codons(line, transcript_id, psite)
                if psite_info[psite] == 'cds' and int(line['psite_from_stop']) < -2:
                    codon_center = codon_centers[transcript_id][psite]
                    psite_info[codon_center] = 'cds'
                prev_transcript = transcript_id
    if transcript_id in codon_lists:
        calculate_backgound_pauses(transcript_id)
    print(j)
    b.close()
c.close()

d=open('pauses/' + project_name + '/transcripts_pass_threshold.tsv', 'w')
writer = csv.writer(d, delimiter = '\t')
writer.writerow(['transcript_id', 'samples'])
for transcript in pass_threshold:
    writer.writerow([transcript, ','.join(list(pass_threshold[transcript]))])
d.close()

