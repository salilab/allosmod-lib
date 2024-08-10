"""Make Modeller restraints between protein and sugars."""

import optparse
import subprocess
import allosmod.util
import sys


def get_glyc_restraint(pdb, allosmod_py):
    # No Python implementation yet - call the original script:
    script = allosmod.util.get_data_file('get_glyc_restraint.sh')
    sys.exit(subprocess.call([script, pdb, allosmod_py]))


def parse_args():
    usage = """%%prog [opts] <pdb> <allosmod.py>

Generate Modeller restraints between protein and sugars. Any protein-sugar
interactions on <pdb> defined by patches in <allosmod.py> are removed (the file
is modified in place) and roughly equivalent Modeller restraints are printed.
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    return args


def main():
    pdb, allosmod_py = parse_args()
    get_glyc_restraint(pdb, allosmod_py)


if __name__ == '__main__':
    main()
