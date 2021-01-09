"""Make a PDB file by averaging two others"""

from __future__ import print_function, absolute_import
import optparse
import subprocess
import allosmod.util
import sys


def getavgpdb(pdb1, pdb2, code1, code2):
    # No Python implementation yet - call the original script:
    script = allosmod.util.get_data_file('getavgpdb2.sh')
    sys.exit(subprocess.call([script, pdb1, pdb2, code1, code2,
                              sys.executable]))


def parse_args():
    usage = """%%prog [opts] <pdb1> <pdb2> <code1> <code2>

Make a PDB file by averaging two others (<pdb1> and <pdb2>, with sequences
in align.ali given by codes <code1> and <code2> respectively).
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 4:
        parser.error("incorrect number of arguments")
    return args


def main():
    pdb1, pdb2, code1, code2 = parse_args()
    getavgpdb(pdb1, pdb2, code1, code2)


if __name__ == '__main__':
    main()
