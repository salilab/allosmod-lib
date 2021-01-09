"""Translate a PDB file"""

from __future__ import print_function
import optparse


def translatepdb(pdb_file, dx, dy, dz):
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                print(line[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + line[54:],
                      end='')


def parse_args():
    usage = """%prog <PDB file> <dx> <dy> <dz>

Translate all ATOM and HETATM records in a PDB file, and write them
to standard output.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 4:
        parser.error("incorrect number of arguments")
    return args[0], float(args[1]), float(args[2]), float(args[3])


def main():
    pdb_file, dx, dy, dz = parse_args()
    translatepdb(pdb_file, dx, dy, dz)


if __name__ == '__main__':
    main()
