"""Get the center of mass of a PDB file"""

from __future__ import print_function
import optparse
import allosmod.util

class CenterOfMassPDBParser(allosmod.util.PDBParser):
    def get_cofm(self, fh):
        x = y = z = 0.
        nr = 0
        for line in self.parse(fh):
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
    c = CenterOfMassPDBParser(allosmod.util.atom_hetatm_filter)
    print("%8.3f%8.3f%8.3f" % c.get_cofm(open(pdb_file)))

if __name__ == '__main__':
    main()
