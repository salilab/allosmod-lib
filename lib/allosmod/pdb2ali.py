"""Convert a PDB file to a Modeller alignment file."""

import optparse

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

if __name__ == '__main__':
    main()
