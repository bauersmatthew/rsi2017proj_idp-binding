import sys
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
import Bio.PDB.Polypeptide as PP
import traceback
import regex

seqdump = []
def maybe_do_seqdump():
    global seqdump
    if 'd' not in sys.argv[1] or not seqdump:
        return
    with open('ssspy_seqdump.txt', 'w') as fout:
        for seqs in seqdump:
            fout.write('{}\tNOT_IN\t{}\n'.format(seqs[0], seqs[1]))

def strifyl(iterable):
    return [str(x) for x in iterable]

def permissive_strsearch(query, body, permissivity=5):
    matches = regex.findall(
        '({}){{e<={}}}'.format(query, permissivity),
        body)
    if len(matches) > 1:
        raise RuntimeError('Multiple IDR matches!')
    elif len(matches) == 0:
        raise RuntimeError('No IDR matches!')
    else:
        return body.index(matches[0]), len(matches[0])

def get_sequence(structure):
    aaseq = ''
    for res in structure.get_residues():
        if res.get_resname() in PP.standard_aa_names:
            aaseq += (PP.three_to_one(res.get_resname()))
    return aaseq

def get_offset(structure, seq):
    pp_seq = get_sequence(structure)
    if len(pp_seq) < len(seq):
        seqdump.append((seq, pp_seq))
        raise RuntimeError('Structure smaller than IDR!')
    try:
        return permissive_strsearch(str(seq), str(pp_seq))
    except:
        seqdump.append((seq, pp_seq))
        raise

def flatten(l):
    r = []
    for i in l:
        if type(i) in [list, tuple]:
            r += flatten(i)
        else:
            r.append(i)
    return r

def get_box(structure, start, end):
    # assume one model one chain
    chain = structure[0].get_list()[0]
    residues = chain.get_list()[start:end]
    atoms = flatten([res.get_list() for res in residues])

    min_corner = [None, None, None] # ATOM COORDS
    max_corner = [None, None, None] # ATOM COORDS

    for atom in atoms:
        coord = atom.get_coord()
        for c in range(3):
            if min_corner[c] is None or coord[c] < min_corner[c]:
                min_corner[c] = coord[c]
            if max_corner[c] is None or coord[c] > max_corner[c]:
                max_corner[c] = coord[c]

    box_center = [(min_corner[i]+max_corner[i])/2.0 for i in range(3)]
    box_widths = [(max_corner[i]-min_corner[i]) for i in range(3)]

    return (*box_center, *box_widths)

if __name__ == '__main__':
    if '-h' in sys.argv or '--help' in sys.argv:
        sys.stderr.write(
            ('Usage: python select_search_spaces.py [tda0] idr_locs.tsv '
             'prot_seqs.fa 1.pdb 2.pdb ... > search_spaces.tsv\n'))
        sys.exit(1)

    # load file describing location of the disordered regions
    idr_location_info = {}
    with open(sys.argv[2]) as fin:
        for line in fin:
            fields = line.rstrip().split()
            if fields:
                idr_location_info[fields[0]] = (int(fields[1]), int(fields[2]))

    # load fasta of protein sequences
    prot_seqs = {}
    for rec in SeqIO.parse(sys.argv[3], 'fasta'):
        prot_seqs[rec.id] = rec.seq

    # process each input file
    pdbparser = PDBParser()
    for fpath in sys.argv[4:-1]:
        base = '.'.join(fpath.split('/')[-1].split('.')[:-1])

        idrloc = idr_location_info[base]
        idrseq = prot_seqs[base][idrloc[0]-1:idrloc[1]]

        structure = pdbparser.get_structure(base, fpath)

        try:
            pdb_idr_begin, pdb_idr_len = get_offset(structure, idrseq)
            pdb_idr_end = pdb_idr_begin+pdb_idr_len

            if (pdb_idr_begin-5 <= 0) or \
               (pdb_idr_end+5 > len(get_sequence(structure))):
                raise RuntimeError('Not enough space around the IDR!')

            nside = get_box(structure, pdb_idr_begin-5, pdb_idr_begin)
            cside = get_box(structure, pdb_idr_end, pdb_idr_end+5)
            sys.stdout.write(
                '{}\t{}\t{}\n'.format(
                    base,
                    '\t'.join(strifyl(nside)),
                    '\t'.join(strifyl(cside))))
        except Exception as exc:
            sys.stderr.write('On file {}...\n'.format(fpath))
            sys.stderr.write('{}\n'.format(exc))
            if 't' in sys.argv[1]:
                sys.stderr.write('Traceback:\n')
                traceback.print_tb(exc.__traceback__)
            if 'a' in sys.argv[1]:
                sys.stderr.write('Quitting...\n')
                maybe_do_seqdump()
                sys.exit(2)
            else:
                sys.stderr.write('Skipping...\n\n')
                continue

    maybe_do_seqdump()
    sys.exit(0)
