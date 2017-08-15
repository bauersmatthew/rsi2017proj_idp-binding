# USAGE:
# python perform_docking.py order.txt search_spaces.tsv vf_models/ bc_models/ \
#        path_to_vina out_dir/
import sys
import os
import glob
import subprocess

def flatten(l):
    r = []
    for i in l:
        if type(i) in (list, tuple):
            r += flatten(i)
        else:
            r.append(i)
    return r

def load_tsv(path):
    with open(path) as fin:
        return [l.rstrip().split('\t') for l in fin if l.rstrip()]

def run_vina(vf, bc, ss, side):
    stdoe = open(
        '{}/{}/{}/{}/stdoe.txt'.format(sys.argv[6], vf, bc, side),
        'w')
    subprocess.call(
        (sys.argv[5],
         '--receptor', '{}/{}.pdbqt'.format(sys.argv[3], vf),
         '--ligand', bcpath,
         '--center_x', ss[0],
         '--center_y', ss[1],
         '--center_z', ss[2],
         '--size_x', ss[3],
         '--size_y', ss[4],
         '--size_z', ss[5],
         '--out', '{o}/{vf}/{bc}/{s}/{vf}_{bc}_{s}.pdbqt'.format(
             o=sys.argv[6], vf=vf, bc=bc, s=side),
         '--log', '{}/{}/{}/{}/log.txt'.format(sys.argv[6], vf, bc, side)),
        stdout=stdoe,
        stderr=stdoe)
    stdoe.close()

# load order
run_order = flatten(load_tsv(sys.argv[1]))

# load search spaces
ss_data = load_tsv(sys.argv[2])
search_spaces = {}
for rec in ss_data:
    search_spaces[rec[0]] = rec[1:]

# create main output directory
os.mkdir(sys.argv[6])

# go through each vf in order
for vf in run_order:
    # check that we have a search space for it
    if vf not in search_spaces:
        continue
    # alert user of progress
    sys.stderr.write('Processing VF {}...\n'.format(vf))
    # create subdir for this vf
    os.mkdir('{}/{}'.format(sys.argv[6], vf))
    # go through each bioactive compound
    bcpaths = glob.glob('{}/*'.format(sys.argv[4]))
    for i, bcpath in enumerate(bcpaths):
        # update status maybe
        if i%299 == 0: # ~3000 bcs as of right now, 300 ~ 10%
            sys.stderr.write('\t~{}%\n'.format(i//299))
        # create subsubdir for this bc
        bc = '.'.join(bcpath.split('/')[-1].split('.')[:-1])
        os.mkdir('{}/{}/{}'.format(sys.argv[6], vf, bc))
        ss_both = search_spaces[vf]
        # run vina on n-terminal side box
        os.mkdir('{}/{}/{}/{}'.format(sys.argv[6], vf, bc, 'nside'))
        run_vina(vf, bc, ss_both[0:6], 'nside')
        # run vina on c-terminal side box
        os.mkdir('{}/{}/{}/{}'.format(sys.argv[6], vf, bc, 'cside'))
        run_vina(vf, bc, ss_both[6:12], 'cside')
    # alert user of progress
    sys.stderr.write('\tDone!\n')

# done; alert me
sys.stderr.write('All finished!\n')
try:
    import cypy_email as email
    email.send_gmail(
        'bauer.s.matthew', 'Blue&Green77',
        '7753386715@mms.att.net',
        'Docking done!', '')
except:
    sys.stderr.write('SMS notification failed.\n')
