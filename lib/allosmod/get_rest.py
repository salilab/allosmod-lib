"""Generate Modeller restraints for glycosylation."""

import argparse
import subprocess
import allosmod.util
import sys


def get_rest(pdb):
    # No Python implementation yet - call the original script:
    script = allosmod.util.get_data_file('get_rest.sh')
    sys.exit(subprocess.call([script, pdb]))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate Modeller restraints for glycosylation.")
    parser.add_argument("pdb", help="PDB file")

    args = parser.parse_args()
    return args.pdb


def main():
    pdb = parse_args()
    get_rest(pdb)


if __name__ == '__main__':
    main()
