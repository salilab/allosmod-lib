"""Modify PDB residues that Modeller cannot recognize."""

from __future__ import print_function
import optparse

def pdb_fix_res(pdb_file):
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if line[17:20] in ('HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP'):
                    line = 'ATOM  ' + line[6:17] + 'HIS' + line[20:]
                if line[17:20] == 'MSE':
                    line = 'ATOM  ' + line[6:17] + 'MET' + line[20:]
                    if line[12:14] == 'SE':
                        line = line[:12] + ' SD' + line[15:]
                print(line, end='')

def parse_args():
    usage = """%prog <PDB file>

Modify PDB residues that Modeller cannot recognize, and write the new PDB file
to standard output.

The following substitutions are made:
- Everything except ATOM and HETATM records is stripped from the file.
- MSE residues (either ATOM or HETATM) are mapped to HIS ATOM records,
  and any SE atoms in those residues are mapped to SD.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    return args[0]

def main():
    pdb_file = parse_args()
    pdb_fix_res(pdb_file)

if __name__ == '__main__':
    main()
