import sys
import allosmod.util

def _enumerate_sequence(seq):
    """Pass through the sequence and yield alignment position, residue index,
       and the one-letter code"""
    ires = -1
    for alnpos, res in enumerate(seq.primary):
        if res != "-" and res != "/":
            ires += 1
        yield alnpos, ires, res

def _get_insertion_point(seq, lres):
    for alnpos, ires, res in _enumerate_sequence(seq):
        if ires == lres:
            return alnpos

def insert_gap(fname, sequence, from_res, to_res, fh=sys.stdout):
    """Insert a gap in the named alignment file.
       A gap is inserted in the PIR format alignment file `fname` such that
       residues `from_res` through `to_res` in the sequence with index
       `sequence` are now misaligned with all other sequences. The resulting
       alignment file is printed out."""
    p = allosmod.util.PIRFile()
    with open(fname) as f:
        seqs = [s for s in p.read(f)]
    insertion_point = _get_insertion_point(seqs[sequence], from_res)
    
    for iseq, seq in enumerate(seqs):
        insert = None
        for alnpos, ires, res in _enumerate_sequence(seq):
            if alnpos == insertion_point and iseq != sequence:
                insert = alnpos
            elif ires == to_res and iseq == sequence:
                insert = alnpos + 1
        seq.primary = seq.primary[:insert] + "-" * (to_res - from_res + 1) \
                      + seq.primary[insert:]
        p.write(fh, seq)
