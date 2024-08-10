"""Rotate a PDB file"""

import optparse
import math


def rotate_coordinate(x, y, z, dx, dy, dz):
    xa = x
    ya = math.cos(dx) * y - math.sin(dx) * z
    za = math.sin(dx) * y + math.cos(dx) * z

    xb = math.sin(dy) * za + math.cos(dy) * xa
    yb = ya
    zb = math.cos(dy) * za - math.sin(dy) * xa

    xc = math.cos(dz) * xb - math.sin(dz) * yb
    yc = math.sin(dz) * xb + math.cos(dz) * yb
    zc = zb
    return xc, yc, zc


def rotatepdb(pdb_file, dx, dy, dz):
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                x, y, z = rotate_coordinate(x, y, z, dx, dy, dz)
                print(line[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + line[54:],
                      end='')


def parse_args():
    usage = """%prog <PDB file> <dx> <dy> <dz>

Rotate all ATOM and HETATM records in a PDB file, and write them
to standard output. The coordinates are first rotated dx degrees about the
x axis, then dy degrees about the y axis, then dz degrees about the z axis.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 4:
        parser.error("incorrect number of arguments")
    deg2rad = math.pi / 180.
    return (args[0], float(args[1]) * deg2rad, float(args[2]) * deg2rad,
            float(args[3]) * deg2rad)


def main():
    pdb_file, dx, dy, dz = parse_args()
    rotatepdb(pdb_file, dx, dy, dz)


if __name__ == '__main__':
    main()
