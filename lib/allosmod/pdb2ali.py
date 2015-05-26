"""Convert a PDB file to a Modeller alignment file."""

from __future__ import print_function
import optparse
import sys
import collections

rescodes = {'HEM': 'h', 'PPI': 'f', 'RIB': 'r', 'ALA': 'A', 'ARG': 'R',
            'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q',
            'GLY': 'G', 'HIS': 'H', 'HID': 'H', 'HIE': 'H', 'HIP': 'H',
            'HSD': 'H', 'HSE': 'H', 'HSP': 'H', 'ILE': 'I', 'LEU': 'L',
            'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'ADE': 'a', 'CYT': 'c', 'GUA': 'g', 'THY': 's', 'URA': 'u',
            'A': 'a', 'C': 'c', 'G': 'g', 'T': 's', 'U': 'u',
            'DA': 'e', 'DC': 'j', 'DG': 'l', 'DT': 't', 'DU': 'v', 'PPP': '-'}

Residue = collections.namedtuple('Residue', ['resnam', 'chain', 'resnum'])

def get_residues(pdb_file):
    last = None
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resnam = line[17:20]
                chain = line[21]
                resnum = line[22:27]
                if (resnam, chain, resnum) != last:
                    last = (resnam, chain, resnum)
                    if resnam != 'HOH':
                        yield Residue(rescodes.get(resnam.strip(), '.'),
                                      chain, resnum)

class SequencePrinter(object):
    def __init__(self, fh):
        self.fh = fh
        self.num_printed = 0
    def __call__(self, c):
        print(c, end='', file=self.fh)
        self.num_printed += 1
        if self.num_printed % 75 == 0:
            print(file=self.fh)

def any_empty_chain_ids(residues):
    for r in residues:
        if r.chain == ' ':
            return True

def rewrite_chain_ids(pdb_file):
    """Replace any empty chain IDs with '@'"""
    with open(pdb_file) as fh:
        lines = fh.readlines()
    with open(pdb_file, 'w') as fh:
        for line in lines:
            if (line.startswith('ATOM') or line.startswith('HETATM')) \
               and line[21] == ' ':
                line = line[:21] + '@' + line[22:]
            fh.write(line)

def pdb2ali(pdb_file, fh=sys.stdout):
    residues = list(get_residues(pdb_file))
    if any_empty_chain_ids(residues):
        rewrite_chain_ids(pdb_file)
        residues = list(get_residues(pdb_file))
    print(">P1;" + pdb_file, file=fh)
    nres = len([r for r in residues if r.resnam != '-'])
    print("structureX:%s:%s:%s:+%d:%s:::-1.00:-1.00" %
          (pdb_file, residues[0].resnum, residues[0].chain, nres,
           residues[-1].chain), file=fh)
    last_chain = None
    p = SequencePrinter(fh)
    for r in residues:
        if r.chain != last_chain:
            if last_chain is not None:
                p('/')
            last_chain = r.chain
        p(r.resnam)
    p('*')
    print(file=fh)

def parse_args():
    usage = """%prog <PDB file>

Convert a PDB file to a Modeller alignment file, which is printed on
standard output.

All hetero atoms are output as BLK residues ("."), except for HOH/water
(which is ignored) and heme (which is written as "h").

If the PDB file contains any chains with no ID, the file is rewritten
with these chains given the '@' ID.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    return args[0]

def main():
    pdb_file = parse_args()
    pdb2ali(pdb_file)

if __name__ == '__main__':
    main()
