"""Set the chain ID on a PDB file."""

from __future__ import print_function
import optparse
import fileinput

def setchain(pdb_file, chain_id, inplace):
    """Replace any empty chain IDs with '@'"""
    for line in fileinput.input(pdb_file, inplace):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            line = line[:21] + chain_id + line[22:]
        print(line, end='')

def parse_args():
    usage = """%prog <PDB file> <chainID>

Set the chain ID on a PDB file, and write it to standard output (or modify
it inplace if --in-place is specified).
"""
    parser = optparse.OptionParser(usage)
    parser.add_option("--in-place", action="store_true", dest="inplace",
                      help="modify the file in place")

    opts, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    return args[0], args[1][:1], opts

def main():
    pdb_file, chain_id, opts = parse_args()
    setchain(pdb_file, chain_id, opts.inplace)

if __name__ == '__main__':
    main()
