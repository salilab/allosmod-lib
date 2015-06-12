"""Get secondary structure with DSSP."""

from __future__ import print_function, absolute_import
import optparse
import subprocess

def check_output(args, stderr=None, retcode=0, input=None, *other):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr,
                         stdin=subprocess.PIPE if input else None,
                         *other)
    stdout, stderr = p.communicate(input)
    if p.returncode != retcode:
        raise OSError("Process %s exited with code %d"
                      % (" ".join(args), p.returncode))
    return stdout

def get_ss(pdb_file):
    out = check_output(["dssp", pdb_file])
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
