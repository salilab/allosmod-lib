"""Get the radius of gyration of a PDB file"""

from __future__ import print_function
import optparse
import allosmod.util
from allosmod.getcofm import CenterOfMassPDBParser
import math

def getrofg(pdb_file):
    c = CenterOfMassPDBParser(allosmod.util.atom_hetatm_filter)
    cofm_x, cofm_y, cofm_z = c.get_cofm(open(pdb_file))
    nr = 0
    r = 0
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                r += math.sqrt((x-cofm_x)**2 + (y-cofm_y)**2 + (z-cofm_z)**2)
                nr += 1
    return r / nr

def parse_args():
    usage = """%prog <PDB file>

Get and print the radius of gyration of the given PDB file.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    return args[0]

def main():
    pdb_file = parse_args()
    print("%6.1f" % getrofg(pdb_file))

if __name__ == '__main__':
    main()
