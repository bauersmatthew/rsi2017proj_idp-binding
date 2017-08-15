def R(w, d): # Real
    class X:
        _w=w
        _d=d
        __class__='Real'
        def __init__(self, s):
            self.val = float(s)
        def __str__(self):
            return '{{:{}.{}f}}'.format(self._w, self._d).format(self.val)
        def __repr__(self):
            return 'R({},{})'.format(self._w, self._d)
    return X

class Record(dict):
    def __str__(self):
        global pdb_format
        line = list(' '*80)
        line[0:6] = list('{:<6s}'.format(self['type']))
        fmt = pdb_format[self['type']]
        for f in fmt:
            templ = '{{:>{}s}}'.format(f[0][1]-f[0][0])
            ins = list(templ.format(str(self[f[1]])))
            line[f[0][0]:f[0][1]] = ins
        return ''.join(line)

# convert continuations (which can be blank) into ints
class Cont:
    def __init__(self, s):
        if not s.strip():
            self.val = 1
        else:
            self.val = int(s)
    def __str__(self):
        if self.val == 1:
            return ''
        else:
            return str(self.val)
    def __eq__(self, other):
        if type(other) == int:
            return self.val == other
        else:
            return self.val == other.val

# optional integer
class OInt:
    def __init__(self, s):
        if not s.strip():
            self.val = -1
        else:
            self.val = int(s)
    def __str__(self):
        if self.val == -1:
            return ''
        else:
            return str(self.val)
    def __eq__(self, other):
        if type(other) == int:
            return self.val == other
        else:
            return self.val == other.val
    def __bool__(self):
        return self.val != -1

keep_record = (((0,80), 'text'),)
pdb_format = {
    'HEADER' : keep_record,
    'TITLE' : keep_record,
    'COMPND' : keep_record,
    'SOURCE' : keep_record,
    'KEYWDS' : keep_record,
    'EXPDTA' : keep_record,
    'AUTHOR' : keep_record,
    'REVDAT' : keep_record,
    'JRNL' : keep_record,
    'REMARK' : (((7,10), 'num', int), ((11,79), 'text')),
    'DBREF' : (((7,11), 'idCode'), ((12,13), 'chainID'),
               ((14,18), 'seqBegin', int), ((18,19), 'insertBegin'),
               ((20,24), 'seqEnd', int), ((24,25), 'insertEnd'),
               ((26,32), 'database'), ((33,41), 'dbAccession'),
               ((42,54), 'dbIdCode'), ((55,60), 'dbseqBegin', int),
               ((60,61), 'dbinsBeg'), ((62,67), 'dbseqEnd', int),
               ((67,68), 'dbinsEnd')),
    'SEQRES' : (((7,10), 'serNum', int), ((11,12), 'chainID'),
                ((13,17), 'numRes', int), ((19,70), 'residues')),
    'HET' : (((7,10), 'hetID'), ((12,13), 'chainID'),
             ((13,17), 'seqNum', int), ((17,18), 'iCode'),
             ((20,25), 'numHetAtoms', int), ((30,70), 'text')),
    'HETNAM' : (((8,10), 'continuation', Cont), ((11,14), 'hetID'),
                ((15,70), 'text')),
    'FORMUL' : (((8,10), 'compNum', int), ((12,15), 'hetID'),
                ((16,18), 'continuation', Cont), ((18,19), 'asterisk'),
                ((19,70), 'text')),
    'HELIX' : (((7,10), 'serNum', int), ((11,14), 'helixID'),
               ((15,18), 'initResName'), ((19,20), 'initChainID'),
               ((21,25), 'initSeqNum', int), ((25,26), 'initICode'),
               ((27,30), 'endResName'), ((31,32), 'endChainID'),
               ((33,37), 'endSeqNum', int), ((37,38), 'endICode'),
               ((38,40), 'helixClass', int), ((40,70), 'comment'),
               ((71,76), 'length', int)),
    'SHEET' : (((7,10), 'strand', int), ((11,14), 'sheetID'),
               ((14,16), 'numStrands', int), ((17,20), 'initResName'),
               ((21,22), 'initChainID'), ((22,26), 'initSeqNum', int),
               ((26,27), 'initICode'), ((28,31), 'endResName'),
               ((32,33), 'endChainID'), ((33,37), 'endSeqNum', int),
               ((37,38), 'endICode'), ((38,40), 'sense', int),
               ((41,45), 'curAtom'), ((45,48), 'curResName'),
               ((49,50), 'curChainId'), ((50,54), 'curResSeq', OInt),
               ((54,55), 'curICode'), ((56,60), 'prevAtom'),
               ((60,63), 'prevResName'), ((64,65), 'prevChainId'),
               ((65,69), 'prevResSeq', OInt), ((69,70), 'prevICode')),
    'LINK' : (((12,16), 'name1'), ((16,17), 'altLoc1'),
              ((17,20), 'resName1'), ((21,22), 'chainID1'),
              ((22,26), 'resSeq1', int), ((26,27), 'iCode1'),
              ((42,46), 'name2'), ((46,47), 'altLoc2'),
              ((47,50), 'resName2'), ((51,52), 'chainID2'),
              ((52,56), 'resSeq2', int), ((56,57), 'iCode2'),
              ((59,65), 'sym1'), ((66,72), 'sym2'),
              ((73,78), 'length', R(5,2))),
    'SITE' : (((7,10), 'seqNum', int), ((11,14), 'siteID'),
              ((15,17), 'numRes', int), ((18,21), 'resName1'),
              ((22,23), 'chainID1'), ((23,27), 'seq1', int),
              ((27,28), 'iCode1'), ((29,32), 'resName2'),
              ((33,34), 'chainID2'), ((34,38), 'seq2', OInt),
              ((38,39), 'iCode2'), ((40,43), 'resName3'),
              ((44,45), 'chainID3'), ((45,49), 'seq3', OInt),
              ((49,50), 'iCode3'), ((51,54), 'resName4'),
              ((55,56), 'chainID4'), ((56,60), 'seq4', OInt),
              ((60,61), 'iCode4')),
    'CRYST1' : keep_record,
    'ORIGX1' : keep_record,
    'ORIGX2' : keep_record,
    'ORIGX3' : keep_record,
    'SCALE1' : keep_record,
    'SCALE2' : keep_record,
    'SCALE3' : keep_record,
    'ATOM' : (((6,11), 'serial', int), ((12,16), 'name'),
              ((16,17), 'altLoc'), ((17,20), 'resName'),
              ((21,22), 'chainID'), ((22,26), 'resSeq', int),
              ((26,27), 'iCode'), ((30,38), 'x', R(8,3)),
              ((38,46), 'y', R(8,3)), ((46,54), 'z', R(8,3)),
              ((54,60), 'occupancy', R(6,2)), ((60,66), 'tempFactor', R(6,2)),
              ((76,78), 'element'), ((78,80), 'charge')),
    'ANISOU' : (((6,11), 'serial', int), ((12,16), 'name'),
                ((16,17), 'altLoc'), ((17,20), 'resName'),
                ((21,22), 'chainID'), ((22,26), 'resSeq', int),
                ((26,27), 'iCode'), ((28,35), 'u00', int),
                ((35,42), 'u11', int), ((42,49), 'u22', int),
                ((49,56), 'u01', int), ((56,63), 'u02', int),
                ((63,70), 'u12', int), ((76,78), 'element'),
                ((78,80), 'charge')),
    'TER' : (((6,11), 'serial', int), ((17,22), 'resName'),
             ((21,22), 'chainID'), ((22,26), 'resSeq'),
             ((26,27), 'iCode')),
    'HETATM' : (((6,11), 'serial', int), ((12,16), 'name'),
                ((16,17), 'altLoc'), ((17,20), 'resName'),
                ((21,22), 'chainID'), ((22,26), 'resSeq', int),
                ((26,27), 'iCode'), ((30,38), 'x', R(8,3)),
                ((38,46), 'y', R(8,3)), ((46,54), 'z', R(8,3)),
                ((54,60),'occupancy',R(6,2)), ((60,66),'tempFactor',R(6,2)),
                ((76,78), 'element'), ((78,80), 'charge')),
    'CONECT' : (((6,11), 'serial1', OInt), ((11,16), 'serial2', OInt),
                ((16,21), 'serial3', OInt), ((21,26), 'serial4', OInt),
                ((26,31), 'serial5', OInt)),
    'MASTER' : (((10,15), 'numRemark', int), ((15,20), 'zeroes'),
                ((20,25), 'numHet', int), ((25,30), 'numHelix', int),
                ((30,35), 'numSheet', int), ((35,40), 'numTurn'),
                ((40,45), 'numSite', int), ((45,50), 'numXform', int),
                ((50,55), 'numCoord', int), ((55,60), 'numTer', int),
                ((60,65), 'numConect', int), ((65,70), 'numSeq', int)),
    'END' : keep_record,
    'SSBOND' : (((7,10), 'serNum', int), ((11,14), 'resName1'),
                ((15,16), 'chainID1'), ((17,21), 'seqNum1', int),
                ((21,22), 'iCode1'), ((25,28), 'resName2'),
                ((29,30), 'chainID2'), ((31,35), 'seqNum2', int),
                ((35,36), 'iCode2'), ((59,65), 'sym1'),
                ((66,72), 'sym2'), ((73,78), 'length', R(5,2))),
    'CISPEP' : (((7,10), 'serNum', int), ((11,14), 'pep1'),
                ((15,16), 'chainID1'), ((17,21), 'seqNum1', int),
                ((21,22), 'iCode1'), ((25,28), 'pep2'),
                ((29,30), 'chainID2'), ((31,35), 'seqNum2', int),
                ((35,36), 'iCode2'), ((43,46), 'modNum', int),
                ((53,59), 'measure', R(6,2))),
    'MODRES' : (((7,11), 'idCode'), ((12,15), 'resName'),
                ((16,17), 'chainID'), ((18,22), 'seqNum', int),
                ((22,23), 'iCode'), ((24,27), 'stdRes'),
                ((29,70), 'comment')),
    'HETSYN' : (((8,10), 'continuation', Cont), ((11,14), 'hetID'),
                ((15,70), 'hetSynonyms')),
    'DBREF1' : (((7,11), 'idCode'), ((12,13), 'chainID'),
                ((14,18), 'seqBegin', int), ((18,19), 'insertBegin'),
                ((20,24), 'seqEnd', int), ((24,25), 'insertEnd'),
                ((26,32), 'database'), ((47,67), 'dbIdCode')),
    'DBREF2' : (((7,11), 'idCode'), ((12,13), 'chainID'),
                ((18,40), 'dbAccession'), ((45,55), 'seqBegin', int),
                ((57,67), 'seqEnd', int)),
    'CAVEAT' : (((8,10), 'continuation', Cont), ((11,15), 'idCode'),
                ((19,79), 'comment')),
    'SEQADV' : (((7,11), 'idCode'), ((12,15), 'resName'),
                ((16,17), 'chainID'), ((18,22), 'seqNum', OInt),
                ((22,23), 'iCode'), ((24,28), 'database'),
                ((29,38), 'dbAccession'), ((39,42), 'dbRes'),
                ((43,48), 'dbSeq', OInt), ((49,70), 'conflict'))}

def load_record(line):
    if len(line) < 80:
        return None

    global pdb_format
    rec_type = line[:6].strip()
    if rec_type not in pdb_format:
        raise RuntimeError('Unrecognized record type: {}'.format(rec_type))

    record = Record()
    record['type'] = rec_type

    rec_fmt = pdb_format[rec_type]
    for field_fmt in rec_fmt:
        cast = field_fmt[2] if len(field_fmt) == 3 else str
        try:
            record[field_fmt[1]] = cast(
                line[field_fmt[0][0]:field_fmt[0][1]])
        except:
            sys.stderr.write('Cast failed on line:\n')
            sys.stderr.write('{}\n'.format(line))
            raise

    return record

def form_cont_groups(recs):
    grab = lambda t, r: (
        [r for r in recs if r['type'] in t] if type(t) is list else
        [r for r in recs if r['type'] == t])
    cgs = []
    # one time multiple lines
    otml = lambda t: cgs.append(grab(t, recs))
    # one time one line
    otol = otml
    # multiple times one line
    def add_sep_internal(recs, cgs):
        for rec in recs:
            cgs.append([rec])
    mtol = lambda t: add_sep_internal(grab(t, recs), cgs)
    # multiple times multiple lines
    def mt_ml_internal(recs, cgs, ng_key):
        group = []
        prev = None
        for rec in recs:
            do_new_group = None
            if len(signature(ng_key).parameters) == 2:
                do_new_group = ng_key(rec, prev)
            else:
                do_new_group = ng_key(rec)
            if do_new_group:
                if group:
                    cgs.append(group)
                group = [rec]
            else:
                group.append(rec)
            prev = rec
        if group:
            cgs.append(group)
    mtml = lambda t, k: mt_ml_internal(grab(t, recs), cgs, k)
    otol('HEADER')
    otol('TITLE')
    mtml('COMPND', lambda r: r['text'][10:16] == 'MOL_ID')
    mtml('SOURCE', lambda r: r['text'][10:16] == 'MOL_ID')
    otml('KEYWDS')
    otml('EXPDTA')
    otml('AUTHOR')
    otml('REVDAT')
    otml('JRNL')
    mtml('REMARK', lambda r: not r['text'].strip())
    mtol('DBREF')
    mtml(['DBREF1', 'DBREF2'], lambda r: r['type'] == 'DBREF1')
    mtol('SEQADV')
    mtml('SEQRES', lambda r: r['serNum'] == 1)
    mtol('MODRES')
    mtol('HET')
    mtml('HETNAM', lambda r: r['continuation'] == 1)
    mtml('HETSYN', lambda r: r['continuation'] == 1)
    mtml('FORMUL', lambda r: r['continuation'] == 1)
    mtol('HELIX')
    mtml('SHEET', lambda r, p: p is None or r['sheetID'] != p['sheetID'])
    mtol('SSBOND')
    mtol('LINK')
    mtol('CISPEP')
    mtml('SITE', lambda r, p: p is None or r['siteID'] != p['siteID'])
    otol('CRYST1')
    otol('ORIGX1')
    otol('ORIGX2')
    otol('ORIGX3')
    otol('SCALE1')
    otol('SCALE2')
    otol('SCALE3')
    mtol(['ATOM', 'ANISOU', 'HETATM', 'TER', 'CONECT'])
    otol('MASTER')
    otol('END')
    # remove empties
    for i, grp in enumerate(cgs):
        if not grp:
            del cgs[i]
    return cgs

# structure:
# pdb[cont_groups[indiv. records]]
def load_pdb(path):
    records = None
    with open(path) as fin:
        records = [load_record(line) for line in fin]
    for i, r in enumerate(records):
        if r is None:
            del records[i]
    return form_cont_groups(records)

def write_pdb(pdb, fpath):
    with open(fpath, 'w') as fout:
        for cgrp in pdb:
            for rec in cgrp:
                fout.write('{}\n'.format(str(rec)))
