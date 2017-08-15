# USAGE: python query_pondr.py prots.fa terse.tsv verbose.txt
# WORKS ON PONDR ON JULY 10 2017

import requests
import sys
from collections import OrderedDict

# load sequence fasta
aaseqs = OrderedDict()
cur_prot = None
growing_aaseq = ''
with open(sys.argv[1]) as fin:
    for line in fin:
        line = line.rstrip()
        if line[0] == '>':
            if cur_prot is not None:
                aaseqs[cur_prot] = growing_aaseq
            cur_prot = line[1:]
            growing_aaseq = ''
        else:
            growing_aaseq += line
    if cur_prot is not None:
        aaseqs[cur_prot] = growing_aaseq
del cur_prot
del growing_aaseq

# talk to PONDR
fout_terse = open(sys.argv[2], 'w')
fout_verbose = open(sys.argv[3], 'w')

pondr_url = 'http://www.pondr.com'
stats_division = '================================VLXT NNP STATISTICS================================'

cur_prot_num = 0
sys.stderr.write('Completed: 0 of {} (0%)'.format(len(aaseqs)))

for prot in aaseqs:
    try:
        cur_prot_num += 1
        payload = [
            ('VLXT', 'on'),
            ('CHStart', ''),
            ('CHEnd', ''),
            ('ProteinName', prot),
            ('AccessionCode', ''),
            ('Sequence', aaseqs[prot]),
            ('stats', 'on'),
            ('submit_result', 'Submit Query'),
            ('.cgifields', 'wcwraw'),
            ('.cgifields', 'stats'),
            ('.cgifields', 'graphic'),
            ('.cgifields', 'CH'),
            ('.cgifields', 'VSL2'),
            ('.cgifields', 'VLXT'),
            ('.cgifields', 'CDF'),
            ('.cgifields', 'XL1'),
            ('.cgifields', 'VL3'),
            ('.cgifields', 'seq'),
            ('.cgifields', 'CAN')]
        req = requests.post(pondr_url, data=payload)
        lines = req.text.split('\n')
        stats_start_lnnum = lines.index(stats_division)+1 #inclusive
        stats_end_lnnum = lines.index('</PRE>', stats_start_lnnum) #exclusive
        stats_lines = lines[stats_start_lnnum:stats_end_lnnum]
        stats_text = '\n'.join(stats_lines)
        fout_verbose.write(
            '====BEGIN====\nGene: {}\n{}\n====END====\n\n'.format(
                prot, stats_text))

        # extract terse info
        percent_disordered = None
        longseg_begin = 0 #inclusive 1-based
        longseg_end = 0 #inclusive!! 1-based
        avg_strength = 0
        for line in stats_lines:
            if 'Overall percent disordered' in line:
                percent_disordered = line[line.find(':')+2:line.find('\t')]
            elif 'Predicted disorder segment' in line:
                bracketed = line[line.find('['):line.find('\t')]
                brack_spl = bracketed.split('-')
                begin = int(brack_spl[0][1:-1])
                end = int(brack_spl[1][1:-1])
                strength = float(line[line.find('=')+2:])
                if (((end-begin+1) > (longseg_end-longseg_begin+1)) or
                    ((end-begin+1) == (longseg_end-longseg_begin+1) and
                        (strength > avg_strength))):
                    longseg_begin = begin
                    longseg_end = end
                    avg_strength = strength

        fout_terse.write('{}\t{}\t{}\t{}\t{}\n'.format(
            prot,
            percent_disordered,
            longseg_begin, longseg_end,
            avg_strength))
                    
        sys.stderr.write('\rCompleted: {} of {} ({}%)'.format(
            cur_prot_num,
            len(aaseqs),
            int((cur_prot_num/len(aaseqs))*100.0)))
    except:
        fout_verbose.write(
            '====BEGIN====\nGene: {}\n{}\n====END====\n\n'.format(
                prot, 'EXCEPTION CAUGHT! (MANUALLY RETRY?)'))
        fout_terse.write('{}{}\n'.format(prot, '\tERR'*4))
        # silently continue... (lol idek)

sys.stderr.write('\n')
fout_terse.close()
fout_verbose.close()
