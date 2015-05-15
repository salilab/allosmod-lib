"""Set the chain ID on a PDB file."""

from __future__ import print_function
import optparse

def setchain(pdb_file, chain_id):
    """Replace any empty chain IDs with '@'"""
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                line = line[:21] + chain_id + line[22:]
            print(line, end='')

def parse_args():
    usage = """%prog <PDB file> <chainID>

Set the chain ID on a PDB file, and write it to standard output.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    return args[0], args[1][:1]

def main():
    pdb_file, chain_id = parse_args()
    setchain(pdb_file, chain_id)

if __name__ == '__main__':
    main()
