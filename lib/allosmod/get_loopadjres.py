"""Strengthen residues adjacent to loops"""

from __future__ import print_function, absolute_import
import optparse
import subprocess
import allosmod.util


def get_loopadjres():
    # No Python implementation yet - call the original script:
    script = allosmod.util.get_data_file('get_loopadjres.sh')
    subprocess.check_call([script])


def parse_args():
    usage = """%%prog [opts]

Strengthen residues adjacent to loops.
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")
    return


def main():
    parse_args()
    get_loopadjres()


if __name__ == '__main__':
    main()
