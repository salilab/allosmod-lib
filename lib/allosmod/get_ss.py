"""Get secondary structure with DSSP."""

from __future__ import print_function, absolute_import
import optparse
import subprocess

def get_ss(pdb_file):
    p = subprocess.Popen(["mkdssp", pdb_file], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True)
    out, err = p.communicate()
    if p.returncode != 0:
        # DSSP considers a lack of amino acids to be an error, but we don't
        # (e.g. the system could be a piece of RNA/DNA). Return an empty
        # set of SS info in this case
        if 'empty protein, or no valid complete residues' in err:
            out = ''
        else:
            raise OSError("mkdssp %s exited with code %d, output %s, error %s"
                          % (pdb_file, p.returncode, out, err))
    start = False
    for line in out.split("\n"):
        if "RESIDUE AA STRUCTURE" in line:
            start = True
        elif start and len(line) > 16 and line[11] != ' ':
            yield line[16] if line[16] != ' ' else '-'

def parse_args():
    usage = """%prog <PDB file>

Get secondary structure of <PDB file> with DSSP, and print the secondary
structure type for each residue, one per line.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    return args[0]

def main():
    pdb_file = parse_args()
    for res in get_ss(pdb_file):
        print(res)

if __name__ == '__main__':
    main()
