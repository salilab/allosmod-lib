"""Translate a PDB file"""

import argparse


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
    parser = argparse.ArgumentParser(
        description="Translate all ATOM and HETATM records in a PDB file, "
                    "and write them to standard output.")
    parser.add_argument("pdb", help="PDB file")
    parser.add_argument("dx", type=float, help="x offset")
    parser.add_argument("dy", type=float, help="y offset")
    parser.add_argument("dz", type=float, help="z offset")
    return parser.parse_args()


def main():
    args = parse_args()
    translatepdb(args.pdb, args.dx, args.dy, args.dz)


if __name__ == '__main__':
    main()
