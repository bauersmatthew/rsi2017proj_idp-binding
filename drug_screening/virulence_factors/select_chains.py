import sys
from inspect import signature
import copy
from pdb_parser import *

def remove_unwanted_records(pdb, keep_chains):
    kill = []
    # find records/cont groups that we want to kill
    for i, cgrp in enumerate(pdb):
        t = cgrp[0]['type']

        if (t in ('DBREF', 'SEQADV', 'SEQRES', 'HET',
                  'ATOM', 'ANISOU', 'TER', 'HETATM',
                  'MODRES', 'DBREF1', 'DBREF2') and
            cgrp[0]['chainID'] not in keep_chains):
            kill.append(i)

        elif (t in ('COMPND',)):
            chain_text = None
            for rec in cgrp:
                if 'CHAIN' in rec['text']:
                    chain_text = rec['text']
                    break
            else:
                raise RuntimeError('COMPND block has no CHAIN entry!')
            compnd_chains = [x.strip()[0] for x in
                             chain_text.split(':')[1].split(',')]
            for i2, chid in enumerate(compnd_chains):
                if chid not in keep_chains:
                    del i2
            if compnd_chains:
                # only some chains were bad
                cgrp[2]['text'] = '{:<80s}'.format(
                    'COMPND   3 CHAIN: {};'.format(', '.join(compnd_chains)))
            else:
                # all chains were bad
                kill.append(i)

        elif (t in ('HELIX',) and
              (cgrp[0]['initChainID'] not in keep_chains or
               cgrp[0]['endChainID'] not in keep_chains)):
            kill.append(i)

        elif (t in ('SHEET',)):
            good = True
            for rec in cgrp:
                if (rec['initChainID'] not in keep_chains or
                    rec['endChainID'] not in keep_chains):
                    good = False
                    break
            if not good:
                kill.append(i)

        elif (t in ('LINK', 'SSBOND', 'CISPEP') and
              (cgrp[0]['chainID1'] not in keep_chains or
               cgrp[0]['chainID2'] not in keep_chains)):
            kill.append(i)

        elif (t in ('SITE') and
              (cgrp[0]['chainID1'] not in keep_chains or
               cgrp[0]['chainID2'] not in keep_chains or
               (cgrp[0]['chainID3'].strip() and
                cgrp[0]['chainID3'] not in keep_chains) or
               (cgrp[0]['chainID4'].strip() and
                cgrp[0]['chainID4'] not in keep_chains))):
            kill.append(i)

    # kill whatever's on our kill list
    for i, k in enumerate(kill):
        del pdb[k-i]

def clean_pdb(pdb):
    # remove records with references to records that don't exist
    # step 1: collect info
    hetIDs = []
    atom_serials = []
    for cgrp in pdb:
        t = cgrp[0]['type']
        if t == 'HET' and cgrp[0]['hetID'] not in hetIDs:
            hetIDs.append(cgrp[0]['hetID'])
        elif t in ('ATOM', 'HETATM'):
            atom_serials.append(cgrp[0]['serial'])
    # step 2: remove stuff that links nonexistent ids/serials
    kill = []
    for i, cgrp in enumerate(pdb):
        t = cgrp[0]['type']
        if (t in ('HETNAM', 'FORMUL', 'HETSYN') and
            cgrp[0]['hetID'] not in hetIDs):
            if not (t == 'FORMUL' and cgrp[0]['hetID'] == 'HOH'): # HOH is weird
                kill.append(i)
        elif (t in ('CONECT',)):
            if cgrp[0]['serial1'] not in atom_serials:
                kill.append(i)
            dests = [copy.copy(cgrp[0]['serial{}'.format(n)]) for n in range(2,6)]
            dests_valid = []
            for d_ser in dests:
                if d_ser and d_ser in atom_serials:
                    dests_valid.append(d_ser)
            for sernum in range(2,6): # clear all
                cgrp[0]['serial{}'.format(sernum)].val = -1
            for dnum in range(len(dests_valid)): # reinsert valids
                cgrp[0]['serial{}'.format(dnum+2)] = dests_valid[dnum]
    for i, k in enumerate(kill):
        del pdb[k-i]
                    
def transform_ids(pdb):
    consecutive_transformation = (
        lambda d: {b : a+1 for a, b in enumerate(d)})
    strid_to_num = (lambda s:
                    ((ord(s[0])-ord('A'))*26*9)+
                    ((ord(s[1])-ord('A'))*9)+
                    int(s[2]))
    num_to_strid = (lambda n:
                    chr(ord('A')+((n-1)//(26*9)))+
                    chr(ord('A')+(((n-1)//9)%26))+
                    str(((n-1)%9)+1))
    dom_atoms = []
    dom_helix_nums = []
    dom_helix_ids = []
    dom_sheet_ids = []
    for cgrp in pdb:
        t = cgrp[0]['type']
        if t in ('ATOM', 'HETATM', 'TER'):
            dom_atoms.append(cgrp[0]['serial'])
        elif t == 'HELIX':
            dom_helix_nums.append(cgrp[0]['serNum'])
            dom_helix_ids.append(cgrp[0]['helixID'])
        elif t == 'SHEET' and cgrp[0]['sheetID'] not in dom_sheet_ids:
            dom_sheet_ids.append(cgrp[0]['sheetID'])

    tform_atoms = consecutive_transformation(dom_atoms)
    tform_helix_nums = consecutive_transformation(dom_helix_nums)
    tform_helix_ids = consecutive_transformation(dom_helix_ids)
    tform_sheet_ids = consecutive_transformation(dom_sheet_ids)

    for cgrp in pdb:
        t = cgrp[0]['type']
        if t in ('ATOM', 'HETATM', 'TER'):
            cgrp[0]['serial'] = tform_atoms[cgrp[0]['serial']]
        elif t == 'HELIX':
            cgrp[0]['serNum'] = tform_helix_nums[cgrp[0]['serNum']]
            cgrp[0]['helixID'] = num_to_strid(
                tform_helix_ids[cgrp[0]['helixID']])
        elif t == 'SHEET':
            for rec in cgrp:
                rec['sheetID'] = num_to_strid(
                    tform_sheet_ids[rec['sheetID']])

def update_master(pdb):
    nremark = 0
    nhet = 0
    nhelix = 0
    nsheet = 0
    nsite = 0
    ntrans = 0
    natom = 0
    nter = 0
    nconect = 0
    nseqres = 0
    for cgrp in pdb:
        t = cgrp[0]['type']
        if t == 'REMARK':
            nremark += len(cgrp)
        elif t == 'HET':
            nhet += len(cgrp)
        elif t == 'HELIX':
            nhelix += len(cgrp)
        elif t == 'SHEET':
            nsheet += len(cgrp)
        elif t == 'SITE':
            nsite += len(cgrp)
        elif t in ('ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3',
                   'MTRIX1', 'MTRIX2', 'MTRIX3'):
            ntrans += len(cgrp)
        elif t in ('ATOM', 'HETATM'):
            natom += len(cgrp)
        elif t == 'TER':
            nter += len(cgrp)
        elif t == 'CONECT':
            nconect += len(cgrp)
        elif t == 'SEQRES':
            nseqres += len(cgrp)

        elif t == 'MASTER':
            # ASSUME MASTER COMES LAST!!!
            cgrp[0]['numRemark'] = nremark
            cgrp[0]['numHet'] = nhet
            cgrp[0]['numHelix'] = nhelix
            cgrp[0]['numSheet'] = nsheet
            cgrp[0]['numSite'] = nsite
            cgrp[0]['numXform'] = ntrans
            cgrp[0]['numCoord'] = natom
            cgrp[0]['numTer'] = nter
            cgrp[0]['numConect'] = nconect
            cgrp[0]['numSeq'] = nseqres

def select_chains(pdb, chains):
    # not incredibly efficient but whatever honestly
    remove_unwanted_records(pdb, chains)
    clean_pdb(pdb)
    transform_ids(pdb)
    update_master(pdb)
    return pdb

if __name__ == '__main__':
    if len(sys.argv) < 4 or '-h' in sys.argv or '--help' in sys.argv:
        sys.stderr.write(
            'Usage: python select_chains.py <wanted-chains.tsv> <file1.pdb> '
            '<file2.pdb> <....pdb> <outdir/>\n')
        sys.exit(1)

    path_chain_info = sys.argv[1]
    pathl_pdbinputs = sys.argv[2:-1]
    path_pdboutdir = sys.argv[-1]

    # load chain info
    wanted_chains = {}
    with open(path_chain_info) as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            wanted_chains[fields[0]] = fields[1].split(',')

    # process inputs
    for path_rawpdb in pathl_pdbinputs:
        base = '.'.join(path_rawpdb.split('/')[-1].split('.')[:-1])
        try:
            write_pdb(select_chains(load_pdb(path_rawpdb), wanted_chains[base]),
                      '{}/{}.pdb'.format(path_pdboutdir, base))
        except:
            sys.stderr.write('ERROR ON FILE: {}\n'.format(path_rawpdb))
            raise

    sys.exit(0)
