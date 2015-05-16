"""Get the center of mass of a PDB file"""

from __future__ import print_function
import optparse

def getcofm(pdb_file):
    nr = 0
    x = y = z = 0.
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x += float(line[30:38])
                y += float(line[38:46])
                z += float(line[46:54])
                nr += 1
    return x / nr, y / nr, z / nr

def parse_args():
    usage = """%prog <PDB file>

Get and print the center of mass of the given PDB file.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    return args[0]

def main():
    pdb_file = parse_args()
    print("%8.3f%8.3f%8.3f" % getcofm(pdb_file))

if __name__ == '__main__':
    main()
