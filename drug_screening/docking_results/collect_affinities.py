import sys
import os
from collections import defaultdict

def get_affinity(fpath):
    with open(fpath) as fin:
        text = fin.read()
        table_raw = text.split('-----+------------+----------+----------\n')[1]
        table_lines = table_raw.split('\n')
        first_fields = table_lines[0].split()
        try:
            return float(first_fields[1])
        except:
            return None # maybe vina failed or something??

for vf in os.listdir(sys.argv[1]):
    for bc in os.listdir('{}/{}'.format(sys.argv[1], vf)):
        for side in ['n', 'c']:
            fpath = '{}/{}/{}/{}side/log.txt'.format(sys.argv[1], vf, bc, side)
            affinity = None
            try:
                affinity = get_affinity(fpath)
            except:
                affinity = None
            if affinity is None:
                sys.stderr.write(
                    ('Failed to get affinity:\n'
                     '  VF  :  {}\n'
                     '  BC  :  {}\n'
                     '  Side:  {}\n'
                     '  Path:  {}\n'
                     'SKIPPING...\n').format(vf, bc, side, fpath))
                continue
            sys.stdout.write('{}\t{}\t{}\t{}\n'.format(vf, bc, side, affinity))
