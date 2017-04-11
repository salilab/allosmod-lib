"""Get RMSD between two structures."""

from __future__ import print_function, absolute_import
import optparse
import subprocess
import os
import re
import allosmod.util

def make_profit_file(prolist):
    return """multi %s
ATOMS CA
#ATOMS CA,N,C,O
#ATOMS ^N,CA,C,O
#ATOMS *
ALIGN
fit
#MWRITE
quit
""" % prolist

def min_rmsd(file1, file2):
    with allosmod.util.temporary_directory() as d:
        prolist_in = os.path.join(d, "prolist.in")
        with open(prolist_in, 'w') as fh:
            print(file1, file=fh)
            print(file2, file=fh)
        out = allosmod.util.check_output(["profit"],
                                         input=make_profit_file(prolist_in),
                                         universal_newlines=True)
    return re.findall('RMS: ([\d.-]+)', out)[-1]
    
def parse_args():
    usage = """%prog <PDB file 1> <PDB file 2>

Get RMSD between two structures using ProFit, and print it.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    return args

def main():
    file1, file2 = parse_args()
    rms = min_rmsd(file1, file2)
    print("%s %s %s" % (os.path.basename(file1), os.path.basename(file2), rms))

if __name__ == '__main__':
    main()
