"""Generate Modeller restraints for glycosylation."""

from __future__ import print_function, absolute_import
import optparse
import subprocess
import allosmod.util
import sys

def get_rest(pdb):
    # No Python implementation yet - call the original script:
    script = allosmod.util.get_data_file('get_rest.sh')
    sys.exit(subprocess.call([script, pdb]))

def parse_args():
    usage = """%%prog [opts] <pdb>

Generate Modeller restraints for glycosylation.
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    return args[0]

def main():
    pdb = parse_args()
    get_rest(pdb)

if __name__ == '__main__':
    main()
