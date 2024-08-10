"""Get the center of mass of a PDB file"""

import argparse
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
    parser = argparse.ArgumentParser(
        description="Get and print the center of mass of all ATOMs in the "
                    "given PDB file.")
    parser.add_argument("pdb", help="PDB file")

    args = parser.parse_args()
    return args.pdb


def main():
    pdb_file = parse_args()
    c = CenterOfMassPDBParser(allosmod.util.atom_filter)
    print("%8.3f %8.3f %8.3f" % c.get_cofm(open(pdb_file)))


if __name__ == '__main__':
    main()
